use crate::types::ClusterIdx;
use ahash::{HashMap, HashMapExt};

/// A CDS as read from a GFF: the owned record an index is built from.
#[derive(Debug, Clone)]
pub struct GffRecord {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub gene_id: String,
}

/// One CDS of an index: a view into the index's own storage, cheap to copy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GffGene<'a> {
    pub contig: &'a str,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub gene_id: &'a str,
    /// The gene's position in its index.
    pub(crate) index: u32,
}

/// One CDS, 32 bytes. Its names (the gene ID, then its aliases) are the
/// `n_names` names from position `first_name` on.
#[derive(Debug, Clone, Copy)]
struct Record {
    start: u64,
    end: u64,
    first_name: u32,
    /// The record whose aliases answer for this one: the last record with the
    /// same gene ID, which is whose aliases a map keyed by gene ID kept.
    alias_owner: u32,
    contig: u32,
    n_names: u16,
    strand: u8,
}

/// The CDS annotation of one genome, sorted by contig and start.
///
/// Every name of every CDS (its gene ID, then its aliases) is stored once,
/// back to back in one buffer, and found through one sorted table of hashes.
/// The index used to hold each name as its own String in up to six places
/// (the gene, an ID map, an alias map and a label map, keys and values), plus
/// every product and note: on PGAP-style annotations that was about 4.7 MB per
/// genome of 4,000 genes, which is what made memory grow with every genome a
/// run added. Products and notes are only read to recognise mobile elements,
/// so that is worked out once, here, and only the few mobile genes keep their
/// labels.
pub struct GffIndex {
    text: String,
    /// Where each name ends in `text`; each starts where the one before ends.
    name_ends: Vec<u32>,
    records: Vec<Record>,
    contig_names: Vec<Box<str>>,
    contig_ids: HashMap<Box<str>, u32>,
    /// Records of each contig, as a range, by contig id.
    contig_ranges: Vec<(u32, u32)>,
    /// Largest gene end seen so far within each contig, aligned with the
    /// records. Genes are sorted by start, and their ends are not monotone once
    /// one CDS contains another, so a binary search on `end` itself is wrong:
    /// it skipped a CDS that contained a nested one, and the annotated locus
    /// read as having no CDS. The running maximum is monotone.
    prefix_max_end: Vec<u64>,
    /// Longest gene on each contig, by contig id, which bounds how far back an
    /// upstream search has to look.
    max_gene_len: Vec<u64>,
    /// (hash of a name, position of the name that answers for it), one entry
    /// per distinct name, sorted.
    name_index: Vec<(u32, u32)>,
    /// Mobile-element genes by alias owner: the kind and the joined labels.
    mobile: Vec<(u32, &'static str, Box<str>)>,
    hasher: ahash::RandomState,
}

/// What a gene's names and descriptions say about it being a mobile element:
/// the kind, and the labels joined as they were read.
fn mobile_kind_from_labels(labels: &[&str]) -> Option<(&'static str, String)> {
    let joined = labels.join("|");
    let normalized = joined
        .to_ascii_lowercase()
        .replace("%20", " ")
        .replace('_', " ");
    let kind = if normalized.contains("is6110") {
        "IS6110"
    } else if normalized.contains("is1081") {
        "IS1081"
    } else if normalized.contains("insertion sequence")
        || normalized.contains("transposase")
        || normalized.contains("is element")
    {
        "insertion_sequence"
    } else if normalized.contains("integrase") || normalized.contains("phage") {
        "phage_or_integrase"
    } else if normalized.contains("recombinase") {
        "recombinase"
    } else {
        return None;
    };
    Some((kind, joined))
}

fn name_hasher() -> ahash::RandomState {
    ahash::RandomState::with_seeds(
        0x243f_6a88_85a3_08d3,
        0x1319_8a2e_0370_7344,
        0xa409_3822_299f_31d0,
        0x082e_fa98_ec4e_6c89,
    )
}

/// The name at `pos` of a name buffer.
fn name_at<'t>(text: &'t str, name_ends: &[u32], pos: u32) -> &'t str {
    let start = match pos {
        0 => 0,
        _ => name_ends[pos as usize - 1] as usize,
    };
    &text[start..name_ends[pos as usize] as usize]
}

impl GffIndex {
    #[allow(dead_code)]
    pub fn build(records: Vec<GffRecord>) -> Self {
        Self::build_with_aliases(
            records
                .into_iter()
                .map(|record| (record, Vec::new()))
                .collect(),
        )
    }

    pub fn build_with_aliases(records: Vec<(GffRecord, Vec<String>)>) -> Self {
        Self::build_with_aliases_and_labels(
            records
                .into_iter()
                .map(|(record, aliases)| {
                    let labels = aliases.clone();
                    (record, aliases, labels)
                })
                .collect(),
        )
    }

    pub fn build_with_aliases_and_labels(
        mut records: Vec<(GffRecord, Vec<String>, Vec<String>)>,
    ) -> Self {
        // Stable, by contig name then start: the order every lookup and every
        // output that walks a contig relies on.
        records.sort_by(|a, b| a.0.contig.cmp(&b.0.contig).then(a.0.start.cmp(&b.0.start)));

        let mut text = String::new();
        let mut name_ends: Vec<u32> = Vec::new();
        let mut contig_names: Vec<Box<str>> = Vec::new();
        let mut contig_ranges: Vec<(u32, u32)> = Vec::new();
        let mut max_gene_len: Vec<u64> = Vec::new();
        let mut out: Vec<Record> = Vec::with_capacity(records.len());
        let mut prefix_max_end: Vec<u64> = Vec::with_capacity(records.len());
        // Labels of every record whose gene ID is not empty, until the mobile
        // genes have been picked out.
        let mut labels_by_record: Vec<Option<Vec<String>>> = Vec::with_capacity(records.len());

        for (i, (gene, aliases, labels)) in records.into_iter().enumerate() {
            // The genes of one contig are contiguous after the sort, so a
            // contig is new exactly when it is not the one just seen.
            let contig = match contig_names.last() {
                Some(last) if **last == *gene.contig => contig_names.len() as u32 - 1,
                _ => {
                    contig_names.push(gene.contig.into_boxed_str());
                    contig_ranges.push((i as u32, i as u32));
                    max_gene_len.push(0);
                    contig_names.len() as u32 - 1
                }
            };
            contig_ranges[contig as usize].1 = i as u32 + 1;
            let new_contig = contig_ranges[contig as usize].0 == i as u32;
            let running = if new_contig {
                gene.end
            } else {
                prefix_max_end[i - 1].max(gene.end)
            };
            prefix_max_end.push(running);
            let len = gene.end.saturating_sub(gene.start);
            max_gene_len[contig as usize] = max_gene_len[contig as usize].max(len);

            // Names: the gene ID, then each alias once, as the alias map held
            // them. A record with no gene ID has none.
            let first_name = name_ends.len() as u32;
            let mut n_names = 0u16;
            let mut stored_labels = None;
            if !gene.gene_id.is_empty() {
                let mut all_aliases: Vec<&str> = vec![gene.gene_id.as_str()];
                for alias in &aliases {
                    if !alias.is_empty() && !all_aliases.contains(&alias.as_str()) {
                        all_aliases.push(alias);
                    }
                }
                let mut all_labels: Vec<String> =
                    all_aliases.iter().map(|s| s.to_string()).collect();
                for label in labels {
                    if !label.is_empty() && !all_labels.contains(&label) {
                        all_labels.push(label);
                    }
                }
                stored_labels = Some(all_labels);
                all_aliases.truncate(u16::MAX as usize);
                for name in &all_aliases {
                    text.push_str(name);
                    name_ends.push(u32::try_from(text.len()).expect("annotation names under 4 GB"));
                }
                n_names = all_aliases.len() as u16;
            }
            labels_by_record.push(stored_labels);
            out.push(Record {
                start: gene.start,
                end: gene.end,
                first_name,
                alias_owner: i as u32,
                contig,
                n_names,
                strand: u8::try_from(gene.strand).unwrap_or(b'+'),
            });
        }

        // Which name answers for each spelling, exactly as the maps did: a gene
        // ID claims its name (a later record with the same ID takes it over),
        // an alias only takes a name nobody has claimed yet. And whose aliases
        // and labels a record answers with: the last record with its gene ID.
        let hasher = name_hasher();
        let mut claims: HashMap<&str, u32> = HashMap::new();
        let mut owner_of_id: HashMap<&str, u32> = HashMap::new();
        for (i, record) in out.iter().enumerate() {
            if record.n_names == 0 {
                continue;
            }
            let gene_id = name_at(&text, &name_ends, record.first_name);
            claims.insert(gene_id, record.first_name);
            owner_of_id.insert(gene_id, i as u32);
            for pos in record.first_name + 1..record.first_name + u32::from(record.n_names) {
                claims.entry(name_at(&text, &name_ends, pos)).or_insert(pos);
            }
        }
        for record in out.iter_mut() {
            if record.n_names > 0 {
                record.alias_owner = owner_of_id[name_at(&text, &name_ends, record.first_name)];
            }
        }
        let mut name_index: Vec<(u32, u32)> = claims
            .into_iter()
            .map(|(name, pos)| (hasher.hash_one(name) as u32, pos))
            .collect();
        name_index.sort_unstable();
        drop(owner_of_id);

        // The mobile kind of each alias owner, from the labels it kept.
        let mut mobile: Vec<(u32, &'static str, Box<str>)> = Vec::new();
        for (i, labels) in labels_by_record.into_iter().enumerate() {
            let Some(labels) = labels else {
                continue;
            };
            if out[i].alias_owner != i as u32 {
                continue;
            }
            let refs: Vec<&str> = labels.iter().map(String::as_str).collect();
            if let Some((kind, joined)) = mobile_kind_from_labels(&refs) {
                mobile.push((i as u32, kind, joined.into_boxed_str()));
            }
        }

        let contig_ids = contig_names
            .iter()
            .enumerate()
            .map(|(id, name)| (name.clone(), id as u32))
            .collect();
        text.shrink_to_fit();
        name_ends.shrink_to_fit();
        GffIndex {
            text,
            name_ends,
            records: out,
            contig_names,
            contig_ids,
            contig_ranges,
            prefix_max_end,
            max_gene_len,
            name_index,
            mobile,
            hasher,
        }
    }

    fn name(&self, pos: u32) -> &str {
        name_at(&self.text, &self.name_ends, pos)
    }

    fn view(&self, i: usize) -> GffGene<'_> {
        let record = &self.records[i];
        GffGene {
            contig: &self.contig_names[record.contig as usize],
            start: record.start,
            end: record.end,
            strand: record.strand as char,
            gene_id: if record.n_names == 0 {
                ""
            } else {
                self.name(record.first_name)
            },
            index: i as u32,
        }
    }

    /// Every CDS in contig and start order.
    pub fn genes(&self) -> impl Iterator<Item = GffGene<'_>> + '_ {
        (0..self.records.len()).map(|i| self.view(i))
    }

    /// The records of one contig, as a range of positions.
    fn contig_range(&self, contig: &str) -> Option<(usize, usize)> {
        let id = *self.contig_ids.get(contig)?;
        let (start, end) = self.contig_ranges[id as usize];
        Some((start as usize, end as usize))
    }

    /// Position of the name that answers for `name`, if the annotation uses it.
    fn name_position(&self, name: &str) -> Option<u32> {
        let hash = self.hasher.hash_one(name) as u32;
        let from = self.name_index.partition_point(|entry| entry.0 < hash);
        self.name_index[from..]
            .iter()
            .take_while(|entry| entry.0 == hash)
            .map(|entry| entry.1)
            .find(|&pos| self.name(pos) == name)
    }

    /// The record a name belongs to. Records hold their names in order, so it
    /// is the last record whose names start at or before it.
    fn record_of_name(&self, pos: u32) -> usize {
        self.records
            .partition_point(|record| record.first_name <= pos)
            - 1
    }

    pub fn get_by_id(&self, gene_id: &str) -> Option<GffGene<'_>> {
        self.name_position(gene_id)
            .map(|pos| self.view(self.record_of_name(pos)))
    }

    /// The record behind a gene this index handed out. The lookups below take
    /// the gene's position, not its ID, so a gene from another index would
    /// read the wrong record.
    fn record_of(&self, gene: &GffGene<'_>) -> &Record {
        let record = &self.records[gene.index as usize];
        debug_assert!(
            record.start == gene.start && record.end == gene.end,
            "gene {} is not from this index",
            gene.gene_id
        );
        record
    }

    /// The names that answer for a gene: its alias owner's, gene ID first.
    fn alias_positions(&self, record: &Record) -> std::ops::Range<u32> {
        let owner = &self.records[record.alias_owner as usize];
        owner.first_name..owner.first_name + u32::from(owner.n_names)
    }

    /// The gene ID and every alias the gene answers to, gene ID first.
    pub fn aliases_for_gene<'a>(&'a self, gene: &GffGene<'a>) -> Vec<&'a str> {
        let record = self.record_of(gene);
        if record.n_names == 0 {
            return vec![gene.gene_id];
        }
        self.alias_positions(record)
            .map(|pos| self.name(pos))
            .collect()
    }

    /// What the gene's names and descriptions say about it being a mobile
    /// element: the kind and the joined labels.
    pub fn mobile_annotation(&self, gene: &GffGene<'_>) -> Option<(&'static str, &str)> {
        let record = self.record_of(gene);
        if record.n_names == 0 {
            return None;
        }
        let owner = record.alias_owner;
        let at = self.mobile.partition_point(|entry| entry.0 < owner);
        self.mobile
            .get(at)
            .filter(|entry| entry.0 == owner)
            .map(|entry| (entry.1, &*entry.2))
    }

    pub fn overlapping(&self, contig: &str, query_start: u64, query_end: u64) -> Vec<GffGene<'_>> {
        let Some((range_start, range_end)) = self.contig_range(contig) else {
            return Vec::new();
        };
        let running_end = &self.prefix_max_end[range_start..range_end];
        let scan_start = range_start + running_end.partition_point(|&end| end <= query_start);

        let mut results = Vec::new();
        for i in scan_start..range_end {
            let record = &self.records[i];
            if record.start >= query_end {
                break;
            }
            if record.end > query_start {
                results.push(self.view(i));
            }
        }
        results
    }

    pub fn flanking(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> (Option<GffGene<'_>>, Option<GffGene<'_>>) {
        let Some((range_start, range_end)) = self.contig_range(contig) else {
            return (None, None);
        };
        let slice = &self.records[range_start..range_end];

        // The first gene starting at or after `end`: starts are sorted.
        let down = slice.partition_point(|g| g.start < end);
        let downstream = (down < slice.len()).then(|| self.view(range_start + down));

        // The gene ending closest before `start`. Ends are not sorted when
        // genes nest, so look back from the last gene starting before `start`,
        // stopping once no earlier gene can end later than the best so far.
        let contig_id = self.records[range_start].contig as usize;
        let longest = self.max_gene_len[contig_id];
        let mut upstream: Option<usize> = None;
        let before = slice.partition_point(|g| g.start < start);
        for k in (0..before).rev() {
            let g = &slice[k];
            if let Some(best) = upstream {
                if g.start.saturating_add(longest) < slice[best].end {
                    break;
                }
            }
            if g.end <= start && upstream.is_none_or(|best| g.end > slice[best].end) {
                upstream = Some(k);
            }
        }

        (upstream.map(|k| self.view(range_start + k)), downstream)
    }

    pub fn context_window(
        &self,
        contig: &str,
        start: u64,
        end: u64,
        center_gene_id: Option<&str>,
        n_neighbors: usize,
    ) -> (Vec<GffGene<'_>>, Vec<GffGene<'_>>) {
        let Some((range_start, range_end)) = self.contig_range(contig) else {
            return (Vec::new(), Vec::new());
        };
        let slice = &self.records[range_start..range_end];
        if slice.is_empty() {
            return (Vec::new(), Vec::new());
        }

        let center_pos = center_gene_id
            .and_then(|gene_id| self.name_position(gene_id))
            .map(|pos| self.record_of_name(pos))
            .and_then(|idx| {
                if idx >= range_start && idx < range_end {
                    Some(idx - range_start)
                } else {
                    None
                }
            })
            .or_else(|| {
                slice
                    .iter()
                    .enumerate()
                    .filter_map(|(i, gene)| {
                        let overlap = gene.end.min(end).saturating_sub(gene.start.max(start));
                        let query_len = end.saturating_sub(start);
                        let gene_len = gene.end.saturating_sub(gene.start);
                        let substantial = overlap > 0
                            && ((query_len > 0 && overlap as f64 / query_len as f64 >= 0.5)
                                || (gene_len > 0 && overlap as f64 / gene_len as f64 >= 0.5));
                        substantial.then_some((i, overlap))
                    })
                    .max_by_key(|(_, overlap)| *overlap)
                    .map(|(i, _)| i)
            });

        let views = |positions: &mut dyn Iterator<Item = usize>| -> Vec<GffGene<'_>> {
            positions.map(|k| self.view(range_start + k)).collect()
        };
        if let Some(center_pos) = center_pos {
            let upstream = views(&mut (0..center_pos).rev().take(n_neighbors));
            let downstream = views(&mut (center_pos + 1..slice.len()).take(n_neighbors));
            return (upstream, downstream);
        }

        let midpoint = start.saturating_add(end.saturating_sub(start) / 2);
        let pos = slice.partition_point(|g| {
            let gene_mid = g.start.saturating_add(g.end.saturating_sub(g.start) / 2);
            gene_mid < midpoint
        });
        let upstream = views(&mut (0..pos).rev().take(n_neighbors));
        let downstream = views(&mut (pos..slice.len()).take(n_neighbors));
        (upstream, downstream)
    }

    /// Bytes this index holds on the heap, roughly.
    #[cfg(test)]
    pub(crate) fn heap_bytes(&self) -> usize {
        self.text.capacity()
            + self.name_ends.capacity() * 4
            + self.records.capacity() * std::mem::size_of::<Record>()
            + self.prefix_max_end.capacity() * 8
            + self.name_index.capacity() * 8
            + self.contig_ranges.capacity() * 8
            + self.max_gene_len.capacity() * 8
            + self
                .contig_names
                .iter()
                .map(|name| 2 * (name.len() + 16) + 4)
                .sum::<usize>()
            + self
                .mobile
                .iter()
                .map(|m| m.2.len() + std::mem::size_of::<(u32, &str, Box<str>)>())
                .sum::<usize>()
    }
}

const NO_CLUSTER: u32 = u32::MAX;

fn pack_cluster(ci: ClusterIdx) -> u32 {
    match u32::try_from(ci.val()) {
        Ok(packed) if packed != NO_CLUSTER => packed,
        _ => panic!("cluster index {} does not fit in 32 bits", ci.val()),
    }
}

fn unpack_cluster(packed: u32) -> Option<ClusterIdx> {
    (packed != NO_CLUSTER).then_some(ClusterIdx(packed as usize))
}

/// The cluster each name of one genome's annotation belongs to.
///
/// This was a map from every name (matrix cell IDs and their variants, gene
/// IDs, aliases) to its cluster, a String and a hash slot each: about 1.1 MB
/// per genome on PGAP-style annotations, held twice once the corrected copy
/// was made. Every lookup a run makes is for a name the annotation uses, so one
/// cluster for each name of the index answers it as the map did, in 4 bytes.
pub struct GeneClusters<'a> {
    gff: Option<&'a GffIndex>,
    /// Aligned with the index's names.
    by_name: Vec<u32>,
    /// What the map held for the empty name, which a gene with no ID asks for.
    empty_name: u32,
}

impl<'a> GeneClusters<'a> {
    /// For a genome with no annotation, where nothing has a cluster.
    pub fn none() -> Self {
        GeneClusters {
            gff: None,
            by_name: Vec::new(),
            empty_name: NO_CLUSTER,
        }
    }

    /// What `map` says about every name of `gff`.
    pub fn from_map(gff: &'a GffIndex, map: &HashMap<String, ClusterIdx>) -> Self {
        let pack = |ci: Option<&ClusterIdx>| ci.map_or(NO_CLUSTER, |ci| pack_cluster(*ci));
        GeneClusters {
            gff: Some(gff),
            by_name: (0..gff.name_ends.len() as u32)
                .map(|pos| pack(map.get(gff.name(pos))))
                .collect(),
            empty_name: pack(map.get("")),
        }
    }

    fn at(&self, pos: u32) -> Option<ClusterIdx> {
        unpack_cluster(self.by_name[pos as usize])
    }

    /// The cluster of a name the annotation uses: what the map held for it.
    /// A name the annotation does not use has none.
    pub fn get(&self, name: &str) -> Option<ClusterIdx> {
        if name.is_empty() {
            return unpack_cluster(self.empty_name);
        }
        let gff = self.gff?;
        gff.name_position(name).and_then(|pos| self.at(pos))
    }

    /// The gene's cluster: its ID's, or else that of the first of its names
    /// that has one.
    pub fn cluster_for_gene(&self, gene: &GffGene<'_>) -> Option<ClusterIdx> {
        let gff = self.gff?;
        let record = gff.record_of(gene);
        if record.n_names == 0 {
            return unpack_cluster(self.empty_name);
        }
        self.at(record.first_name)
            .or_else(|| gff.alias_positions(record).find_map(|pos| self.at(pos)))
    }

    /// Gives each correction's cluster to every name spelled like its gene ID
    /// that has none yet, as adding the gene ID to the map did when it was not
    /// already there. The first correction of a name wins.
    pub fn fill(&mut self, corrections: &[(&str, ClusterIdx)]) {
        let Some(gff) = self.gff else {
            return;
        };
        let mut first: HashMap<&str, u32> = HashMap::new();
        for &(name, ci) in corrections {
            first.entry(name).or_insert(pack_cluster(ci));
        }
        if let Some(&ci) = first.get("") {
            if self.empty_name == NO_CLUSTER {
                self.empty_name = ci;
            }
        }
        for (pos, slot) in self.by_name.iter_mut().enumerate() {
            if *slot == NO_CLUSTER {
                if let Some(&ci) = first.get(gff.name(pos as u32)) {
                    *slot = ci;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;
