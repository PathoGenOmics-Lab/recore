//! Nucleotide Similarity Clustering (NSC) module.
//!
//! Detects and merges Panaroo clusters whose representative sequences are
//! highly similar at the nucleotide level but were split at the amino-acid level.
//! This eliminates inflated accessory genome size due to annotation-level redundancy.
//!
//! Strategy (improvements over PanQC):
//! - MinHash sketching instead of full k-mer sets → O(1) memory per gene
//! - Banded MinHash pre-filter → skip dissimilar pairs, reducing O(n²) to O(n×candidates)
//! - Union-Find for clustering → O(α(n)) per merge
//! - Rayon parallel sketch computation + comparison
//! - Integrated with Recore's PresenceAbsence for in-place matrix correction

use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use anyhow::{Context, Result};
use log::{debug, info};
use rayon::prelude::*;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;

use crate::types::ClusterIdx;

// ── MinHash sketching ──

const SKETCH_SIZE: usize = 128;
/// Number of bands for LSH pre-filter. More bands = fewer false negatives but more comparisons.
const NUM_BANDS: usize = 16;
const ROWS_PER_BAND: usize = SKETCH_SIZE / NUM_BANDS; // 8

/// FNV-1a inspired fast hash for k-mers. Deterministic, no random seed needed.
#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Canonical k-mer hash: min(hash(kmer), hash(revcomp(kmer))).
#[inline]
fn canonical_kmer_hash(kmer: &[u8]) -> u64 {
    let fwd = hash_kmer(kmer);
    let rc = hash_kmer_rc(kmer);
    fwd.min(rc)
}

/// Hash the reverse complement without allocating.
#[inline]
fn hash_kmer_rc(kmer: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in kmer.iter().rev() {
        let comp = match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        };
        h ^= comp as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Bottom-k MinHash sketch: keep the k smallest hash values.
fn minhash_sketch(seq: &[u8], ksize: usize) -> Vec<u64> {
    if seq.len() < ksize {
        return vec![u64::MAX; SKETCH_SIZE];
    }

    // Use a bounded max-heap via sorted vec (SKETCH_SIZE is small)
    let mut sketch = Vec::with_capacity(SKETCH_SIZE + 1);
    let mut max_in_sketch = u64::MAX;

    for i in 0..=(seq.len() - ksize) {
        let kmer = &seq[i..i + ksize];
        // Skip k-mers with ambiguous bases
        if kmer
            .iter()
            .any(|&b| !matches!(b, b'A' | b'a' | b'T' | b't' | b'C' | b'c' | b'G' | b'g'))
        {
            continue;
        }
        let h = canonical_kmer_hash(kmer);
        if sketch.len() < SKETCH_SIZE {
            sketch.push(h);
            if sketch.len() == SKETCH_SIZE {
                sketch.sort_unstable();
                sketch.dedup();
                // Refill if dedup shrank us
                if sketch.len() == SKETCH_SIZE {
                    max_in_sketch = sketch[SKETCH_SIZE - 1];
                }
            }
        } else if h < max_in_sketch {
            // Insert and keep only SKETCH_SIZE smallest
            sketch.push(h);
            sketch.sort_unstable();
            sketch.dedup();
            sketch.truncate(SKETCH_SIZE);
            max_in_sketch = sketch[SKETCH_SIZE - 1];
        }
    }

    sketch.sort_unstable();
    sketch.dedup();
    // Pad with MAX if we didn't get enough unique hashes
    sketch.resize(SKETCH_SIZE, u64::MAX);
    sketch
}

/// Estimate Jaccard similarity from two sorted MinHash sketches.
#[cfg(test)]
fn jaccard_from_sketches(a: &[u64], b: &[u64]) -> f64 {
    let mut i = 0;
    let mut j = 0;
    let mut intersection = 0u64;
    let mut union = 0u64;

    while i < a.len() && j < b.len() {
        if a[i] == u64::MAX && b[j] == u64::MAX {
            break;
        }
        if a[i] == b[j] {
            intersection += 1;
            union += 1;
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            union += 1;
            i += 1;
        } else {
            union += 1;
            j += 1;
        }
    }
    // Count remaining non-MAX elements
    while i < a.len() && a[i] != u64::MAX {
        union += 1;
        i += 1;
    }
    while j < b.len() && b[j] != u64::MAX {
        union += 1;
        j += 1;
    }

    if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    }
}

/// Estimate Jaccard *containment* of `a` in `b`: |A∩B| / |A|.
/// More appropriate than similarity when sequences differ in length.
fn containment_from_sketches(a: &[u64], b: &[u64]) -> f64 {
    let a_size = a.iter().filter(|&&v| v != u64::MAX).count();
    if a_size == 0 {
        return 0.0;
    }

    let mut i = 0;
    let mut j = 0;
    let mut intersection = 0u64;

    while i < a.len() && j < b.len() {
        if a[i] == u64::MAX {
            break;
        }
        if b[j] == u64::MAX {
            break;
        }
        if a[i] == b[j] {
            intersection += 1;
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }

    intersection as f64 / a_size as f64
}

/// Max containment: max(containment(a,b), containment(b,a)).
/// This is what PanQC uses and what makes biological sense for redundancy detection.
fn max_containment_from_sketches(a: &[u64], b: &[u64]) -> f64 {
    containment_from_sketches(a, b).max(containment_from_sketches(b, a))
}

// ── LSH banding for candidate pair pre-filter ──

/// Hash a band (slice of sketch values) to a single u64 for LSH bucketing.
fn band_hash(band: &[u64]) -> u64 {
    let mut h: u64 = 0;
    for &v in band {
        h = h.wrapping_mul(6364136223846793005).wrapping_add(v);
    }
    h
}

/// Find candidate pairs using LSH banding. Two sketches are candidates if they
/// share at least one identical band. This dramatically reduces comparisons.
fn find_candidate_pairs(sketches: &[Vec<u64>]) -> HashSet<(usize, usize)> {
    let n = sketches.len();
    let mut candidates = HashSet::new();

    for band_idx in 0..NUM_BANDS {
        let start = band_idx * ROWS_PER_BAND;
        let end = start + ROWS_PER_BAND;

        // Build bucket: band_hash → list of sketch indices
        let mut buckets: HashMap<u64, Vec<usize>> = HashMap::new();
        for (i, sketch) in sketches.iter().enumerate() {
            let bh = band_hash(&sketch[start..end]);
            buckets.entry(bh).or_default().push(i);
        }

        // All pairs within each bucket are candidates
        for members in buckets.values() {
            if members.len() > 1 {
                for (mi, &a) in members.iter().enumerate() {
                    for &b in &members[mi + 1..] {
                        let pair = if a < b { (a, b) } else { (b, a) };
                        candidates.insert(pair);
                    }
                }
            }
        }
    }

    debug!(
        "LSH banding: {} candidate pairs from {} sketches ({} bands)",
        candidates.len(),
        n,
        NUM_BANDS
    );
    candidates
}

// ── Union-Find ──

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]); // path compression
        }
        self.parent[x]
    }

    fn union(&mut self, x: usize, y: usize) -> bool {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return false;
        }
        // union by rank
        match self.rank[rx].cmp(&self.rank[ry]) {
            std::cmp::Ordering::Less => self.parent[rx] = ry,
            std::cmp::Ordering::Greater => self.parent[ry] = rx,
            std::cmp::Ordering::Equal => {
                self.parent[ry] = rx;
                self.rank[rx] += 1;
            }
        }
        true
    }

    fn components(&mut self) -> HashMap<usize, Vec<usize>> {
        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..self.parent.len() {
            let root = self.find(i);
            groups.entry(root).or_default().push(i);
        }
        groups
    }
}

// ── Public API ──

/// A merged cluster produced by nucleotide similarity clustering.
#[derive(Debug)]
pub struct NsCluster {
    /// Indices of original Panaroo clusters that were merged.
    pub members: Vec<ClusterIdx>,
    /// Representative name (name of the member with highest prevalence).
    pub representative: String,
}

/// Result of nucleotide similarity clustering.
pub struct NsClusterResult {
    /// Merged clusters (only groups with ≥2 members).
    pub merged_clusters: Vec<NsCluster>,
    /// Total number of original clusters that were merged.
    pub total_merged_genes: usize,
    /// Total number of NS clusters formed.
    pub n_clusters: usize,
    /// Map from original ClusterIdx → merged ClusterIdx (representative).
    /// Only contains entries for clusters that were merged.
    pub merge_map: HashMap<ClusterIdx, ClusterIdx>,
}

/// Run nucleotide similarity clustering on the pan-genome reference.
///
/// # Arguments
/// * `pan_ref` - Cluster name → list of allele sequences
/// * `cluster_names` - Ordered cluster names matching ClusterIdx
/// * `cluster_prevalences` - Prevalence (fraction) per cluster for picking representative
/// * `kmer_size` - k-mer size for sketching (default: 21)
/// * `threshold` - Max containment threshold for merging (default: 0.8)
/// * `threads` - Number of rayon threads
pub fn run_nscluster(
    pan_ref: &HashMap<String, Vec<Vec<u8>>>,
    cluster_names: &[String],
    cluster_prevalences: &[f64],
    kmer_size: usize,
    threshold: f64,
    threads: usize,
) -> NsClusterResult {
    // ── 1. Build concatenated sequences per cluster (all alleles joined) ──
    let cluster_seqs: Vec<(&str, Vec<u8>)> = cluster_names
        .iter()
        .filter_map(|name| {
            pan_ref.get(name).map(|alleles| {
                let mut concat = Vec::new();
                for (i, allele) in alleles.iter().enumerate() {
                    if i > 0 {
                        // Separator to avoid k-mers spanning allele boundaries
                        concat.extend_from_slice(&vec![b'N'; kmer_size]);
                    }
                    concat.extend_from_slice(allele);
                }
                (name.as_str(), concat)
            })
        })
        .collect();

    info!(
        "NSC: sketching {} clusters (k={}, sketch_size={})",
        cluster_seqs.len(),
        kmer_size,
        SKETCH_SIZE
    );

    // ── 2. Compute MinHash sketches in parallel ──
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap_or_else(|_| {
            rayon::ThreadPoolBuilder::new()
                .build()
                .expect("rayon fallback")
        });

    let sketches: Vec<Vec<u64>> = pool.install(|| {
        cluster_seqs
            .par_iter()
            .map(|(_, seq)| minhash_sketch(seq, kmer_size))
            .collect()
    });

    // ── 3. LSH banding to find candidate pairs ──
    let candidates = find_candidate_pairs(&sketches);
    info!("NSC: {} candidate pairs from LSH banding", candidates.len());

    // ── 4. Compute max containment for candidates, merge above threshold ──
    let n = cluster_seqs.len();
    let mut uf = UnionFind::new(n);
    let mut n_merged = 0usize;

    // Parallel containment computation
    let pairs_above: Vec<(usize, usize)> = pool.install(|| {
        candidates
            .into_par_iter()
            .filter(|&(a, b)| {
                max_containment_from_sketches(&sketches[a], &sketches[b]) >= threshold
            })
            .collect()
    });

    for (a, b) in &pairs_above {
        if uf.union(*a, *b) {
            n_merged += 1;
        }
    }

    info!(
        "NSC: {} pairs above threshold {:.2}, {} merges performed",
        pairs_above.len(),
        threshold,
        n_merged
    );

    // ── 5. Extract components with ≥2 members ──
    let components = uf.components();
    let mut merged_clusters = Vec::new();
    let mut merge_map: HashMap<ClusterIdx, ClusterIdx> = HashMap::new();
    let mut total_merged_genes = 0usize;

    // Build name→index map
    let name_to_idx: HashMap<&str, usize> = cluster_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    for members in components.values() {
        if members.len() < 2 {
            continue;
        }

        // Map back to ClusterIdx via cluster_seqs names
        let member_cis: Vec<ClusterIdx> = members
            .iter()
            .filter_map(|&seq_idx| {
                let name = cluster_seqs[seq_idx].0;
                name_to_idx.get(name).map(|&ci| ClusterIdx(ci))
            })
            .collect();

        if member_cis.len() < 2 {
            continue;
        }

        // Pick representative: highest prevalence
        let rep_ci = *member_cis
            .iter()
            .max_by(|a, b| {
                cluster_prevalences[a.val()]
                    .partial_cmp(&cluster_prevalences[b.val()])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("member_cis guaranteed non-empty");

        let rep_name = cluster_names[rep_ci.val()].clone();

        for &ci in &member_cis {
            if ci != rep_ci {
                merge_map.insert(ci, rep_ci);
            }
        }

        total_merged_genes += member_cis.len();

        merged_clusters.push(NsCluster {
            members: member_cis,
            representative: rep_name,
        });
    }

    let n_clusters = merged_clusters.len();
    info!(
        "NSC: {} NS clusters formed, {} total genes merged",
        n_clusters, total_merged_genes
    );

    NsClusterResult {
        merged_clusters,
        total_merged_genes,
        n_clusters,
        merge_map,
    }
}

/// Apply NS cluster merges to a presence/absence matrix.
/// For each merged group, the representative cluster gets the union of all member presences.
/// Returns: (clusters_to_remove, presence_updates) where:
/// - clusters_to_remove: ClusterIdx values to drop from the output
/// - presence_updates: (representative_ci, gi) pairs that become present
pub fn apply_ns_merges(
    merge_map: &HashMap<ClusterIdx, ClusterIdx>,
    is_present: impl Fn(ClusterIdx, crate::types::GenomeIdx) -> bool,
    n_genomes: usize,
) -> (
    HashSet<ClusterIdx>,
    Vec<(ClusterIdx, crate::types::GenomeIdx)>,
) {
    let clusters_to_remove: HashSet<ClusterIdx> = merge_map.keys().copied().collect();
    let mut presence_updates = Vec::new();

    // Group by representative
    let mut rep_to_members: HashMap<ClusterIdx, Vec<ClusterIdx>> = HashMap::new();
    for (&member, &rep) in merge_map {
        rep_to_members.entry(rep).or_default().push(member);
    }

    for (rep, members) in &rep_to_members {
        for gi_val in 0..n_genomes {
            let gi = crate::types::GenomeIdx(gi_val);
            if is_present(*rep, gi) {
                continue; // already present
            }
            // Check if ANY member was present → union
            if members.iter().any(|&m| is_present(m, gi)) {
                presence_updates.push((*rep, gi));
            }
        }
    }

    (clusters_to_remove, presence_updates)
}

// ── MMseqs2 linclust engine ──

/// Run nucleotide similarity clustering using MMseqs2 linclust.
/// Requires `mmseqs` binary in PATH.
///
/// MMseqs2 linclust runs in linear time and uses alignment-based identity,
/// giving exact results (not approximate like MinHash).
pub fn run_nscluster_mmseqs2(
    pan_ref: &HashMap<String, Vec<Vec<u8>>>,
    cluster_names: &[String],
    cluster_prevalences: &[f64],
    min_seq_id: f64,
    coverage: f64,
    threads: usize,
) -> Result<NsClusterResult> {
    // Verify mmseqs2 is available
    let version_check = Command::new("mmseqs").arg("version").output();
    if version_check.is_err() {
        anyhow::bail!(
            "mmseqs2 not found in PATH. Install with: conda install -c bioconda mmseqs2\n\
             Or use --nsc-engine minhash (built-in, no dependencies)"
        );
    }

    let tmpdir = tempfile::tempdir().context("Cannot create temp dir for mmseqs2")?;
    let input_fasta = tmpdir.path().join("input.fasta");
    let db_path = tmpdir.path().join("db");
    let cluster_db = tmpdir.path().join("cluster");
    let mmseqs_tmp = tmpdir.path().join("tmp");
    let tsv_path = tmpdir.path().join("cluster.tsv");

    // ── 1. Write concatenated FASTA (longest allele per cluster) ──
    {
        let file =
            std::fs::File::create(&input_fasta).context("Cannot create input FASTA for mmseqs2")?;
        let mut w = BufWriter::new(file);
        for name in cluster_names {
            if let Some(alleles) = pan_ref.get(name) {
                // Use longest allele as representative
                if let Some(longest) = alleles.iter().max_by_key(|a| a.len()) {
                    writeln!(w, ">{}", name)?;
                    w.write_all(longest)?;
                    writeln!(w)?;
                }
            }
        }
        w.flush()?;
    }

    info!("NSC-mmseqs2: creating sequence database...");

    // ── 2. createdb ──
    let output = Command::new("mmseqs")
        .args([
            "createdb",
            &input_fasta.to_string_lossy(),
            &db_path.to_string_lossy(),
        ])
        .output()
        .context("Failed to run mmseqs createdb")?;
    if !output.status.success() {
        anyhow::bail!(
            "mmseqs createdb failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    info!(
        "NSC-mmseqs2: running linclust (min_id={:.0}%, cov={:.0}%)...",
        min_seq_id * 100.0,
        coverage * 100.0
    );

    // ── 3. linclust ──
    std::fs::create_dir_all(&mmseqs_tmp)?;
    let output = Command::new("mmseqs")
        .args([
            "linclust",
            &db_path.to_string_lossy(),
            &cluster_db.to_string_lossy(),
            &mmseqs_tmp.to_string_lossy(),
            "--min-seq-id",
            &min_seq_id.to_string(),
            "-c",
            &coverage.to_string(),
            "--cov-mode",
            "0", // bidirectional coverage
            "--cluster-mode",
            "2", // connected component clustering
            "--threads",
            &threads.to_string(),
            "--seq-id-mode",
            "0", // alignment length based
        ])
        .output()
        .context("Failed to run mmseqs linclust")?;
    if !output.status.success() {
        anyhow::bail!(
            "mmseqs linclust failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    // ── 4. createtsv ──
    let output = Command::new("mmseqs")
        .args([
            "createtsv",
            &db_path.to_string_lossy(),
            &db_path.to_string_lossy(),
            &cluster_db.to_string_lossy(),
            &tsv_path.to_string_lossy(),
        ])
        .output()
        .context("Failed to run mmseqs createtsv")?;
    if !output.status.success() {
        anyhow::bail!(
            "mmseqs createtsv failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    // ── 5. Parse TSV: representative\tmember ──
    let name_to_idx: HashMap<&str, usize> = cluster_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    let mut uf = UnionFind::new(cluster_names.len());
    let file = std::fs::File::open(&tsv_path).context("Cannot open mmseqs2 cluster TSV")?;
    let reader = BufReader::new(file);

    for line_result in reader.lines() {
        let line = line_result?;
        let mut fields = line.split('\t');
        let Some(rep) = fields.next() else { continue };
        let Some(member) = fields.next() else {
            continue;
        };
        if rep == member {
            continue;
        } // self-hit

        if let (Some(&ri), Some(&mi)) = (name_to_idx.get(rep), name_to_idx.get(member)) {
            uf.union(ri, mi);
        }
    }

    // ── 6. Extract components ≥ 2 members ──
    let components = uf.components();
    let mut merged_clusters = Vec::new();
    let mut merge_map: HashMap<ClusterIdx, ClusterIdx> = HashMap::new();
    let mut total_merged_genes = 0usize;

    for members in components.values() {
        if members.len() < 2 {
            continue;
        }

        let member_cis: Vec<ClusterIdx> = members.iter().map(|&i| ClusterIdx(i)).collect();

        // Pick representative: highest prevalence
        let rep_ci = *member_cis
            .iter()
            .max_by(|a, b| {
                cluster_prevalences[a.val()]
                    .partial_cmp(&cluster_prevalences[b.val()])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("member_cis guaranteed non-empty");

        let rep_name = cluster_names[rep_ci.val()].clone();

        for &ci in &member_cis {
            if ci != rep_ci {
                merge_map.insert(ci, rep_ci);
            }
        }

        total_merged_genes += member_cis.len();

        merged_clusters.push(NsCluster {
            members: member_cis,
            representative: rep_name,
        });
    }

    let n_clusters = merged_clusters.len();
    info!(
        "NSC-mmseqs2: {} NS clusters formed, {} total genes merged",
        n_clusters, total_merged_genes
    );

    // Cleanup is automatic (tmpdir drops)
    Ok(NsClusterResult {
        merged_clusters,
        total_merged_genes,
        n_clusters,
        merge_map,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonical_kmer_hash_is_strand_independent() {
        // ATCG and its RC CGAT should give same canonical hash
        let fwd = b"ATCG";
        let rc = b"CGAT";
        assert_eq!(canonical_kmer_hash(fwd), canonical_kmer_hash(rc));
    }

    #[test]
    fn test_minhash_sketch_size() {
        let seq = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        let sketch = minhash_sketch(seq, 11);
        assert_eq!(sketch.len(), SKETCH_SIZE);
    }

    #[test]
    fn test_identical_sequences_perfect_similarity() {
        let seq = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        let s1 = minhash_sketch(seq, 11);
        let s2 = minhash_sketch(seq, 11);
        let j = jaccard_from_sketches(&s1, &s2);
        assert!(
            (j - 1.0).abs() < 0.001,
            "identical seqs should have J≈1.0, got {}",
            j
        );
    }

    #[test]
    fn test_very_different_sequences_low_similarity() {
        let seq1 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let seq2 = b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let s1 = minhash_sketch(seq1, 11);
        let s2 = minhash_sketch(seq2, 11);
        let j = jaccard_from_sketches(&s1, &s2);
        assert!(j < 0.1, "very different seqs should have low J, got {}", j);
    }

    #[test]
    fn test_containment_asymmetric() {
        // Short sequence contained in long one → high containment(short, long)
        let short = b"ATCGATCGATCGATCGATCG";
        let long = b"ATCGATCGATCGATCGATCGTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGG";
        let s_short = minhash_sketch(short, 7);
        let s_long = minhash_sketch(long, 7);
        let c_short_in_long = containment_from_sketches(&s_short, &s_long);
        let c_long_in_short = containment_from_sketches(&s_long, &s_short);
        assert!(
            c_short_in_long > c_long_in_short,
            "short-in-long ({:.2}) should be > long-in-short ({:.2})",
            c_short_in_long,
            c_long_in_short
        );
    }

    #[test]
    fn test_max_containment_symmetric() {
        let a = b"ATCGATCGATCGATCGATCG";
        let b_seq = b"ATCGATCGATCGATCGATCGTTTTTTTTTTTTTTTTTTTT";
        let sa = minhash_sketch(a, 7);
        let sb = minhash_sketch(b_seq, 7);
        assert_eq!(
            max_containment_from_sketches(&sa, &sb),
            max_containment_from_sketches(&sb, &sa)
        );
    }

    #[test]
    fn test_union_find_basic() {
        let mut uf = UnionFind::new(5);
        uf.union(0, 1);
        uf.union(2, 3);
        uf.union(1, 3);
        assert_eq!(uf.find(0), uf.find(3));
        assert_ne!(uf.find(0), uf.find(4));
        let comps = uf.components();
        assert_eq!(comps.len(), 2);
    }

    #[test]
    fn test_band_hash_deterministic() {
        let band = [1u64, 2, 3, 4, 5, 6, 7, 8];
        assert_eq!(band_hash(&band), band_hash(&band));
    }

    #[test]
    fn test_short_sequence_sketch() {
        // Sequence shorter than k → all MAX
        let seq = b"ATCG";
        let sketch = minhash_sketch(seq, 11);
        assert!(sketch.iter().all(|&v| v == u64::MAX));
    }

    #[test]
    fn test_run_nscluster_merges_identical() {
        let mut pan_ref: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
        let seq = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_vec();
        // Two "different" clusters with identical sequences
        pan_ref.insert("clA".into(), vec![seq.clone()]);
        pan_ref.insert("clB".into(), vec![seq.clone()]);
        // One truly different cluster
        pan_ref.insert(
            "clC".into(),
            vec![b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec()],
        );

        let names = vec!["clA".into(), "clB".into(), "clC".into()];
        let prevs = vec![0.8, 0.5, 0.9];

        let result = run_nscluster(&pan_ref, &names, &prevs, 11, 0.8, 1);
        assert_eq!(result.n_clusters, 1, "clA and clB should merge");
        assert_eq!(result.merged_clusters[0].members.len(), 2);
        // Representative should be clA (higher prevalence)
        assert_eq!(result.merged_clusters[0].representative, "clA");
        // clC should NOT be merged
        assert!(!result.merge_map.contains_key(&ClusterIdx(2)));
    }

    #[test]
    fn test_apply_ns_merges() {
        use crate::types::GenomeIdx;

        // clB (idx 1) merges into clA (idx 0)
        let mut merge_map: HashMap<ClusterIdx, ClusterIdx> = HashMap::new();
        merge_map.insert(ClusterIdx(1), ClusterIdx(0));

        // clA present in genome 0, clB present in genome 1
        let is_present = |ci: ClusterIdx, gi: GenomeIdx| -> bool {
            matches!((ci.val(), gi.val()), (0, 0) | (1, 1))
        };

        let (to_remove, updates) = apply_ns_merges(&merge_map, is_present, 2);
        assert!(to_remove.contains(&ClusterIdx(1)));
        // clA should gain presence in genome 1 (from clB)
        assert!(updates.contains(&(ClusterIdx(0), GenomeIdx(1))));
    }
}
