use clap::{Parser, ValueEnum};
use std::path::PathBuf;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum MatrixPolicy {
    /// Paper-grade: correct only strict same-locus or cluster-alias evidence.
    Strict,
    /// Default: correct strict/alias/full-sequence evidence; keep partial/disrupted evidence multistate-only.
    Balanced,
    /// Exploratory: correct all accepted rescues, including partial/disrupted evidence.
    Sensitive,
}

/// Post-hoc rescue of falsely absent genes in pangenome matrices.
#[derive(Parser, Debug)]
#[command(name = "recore", version, about)]
pub struct Args {
    /// Panaroo gene_presence_absence.csv
    #[arg(short = 'p', long)]
    pub presence_absence: PathBuf,

    /// Panaroo pan_genome_reference.fa (representative sequences per cluster)
    #[arg(short = 'r', long)]
    pub pan_reference: PathBuf,

    /// Directory of genome FASTA files
    #[arg(short = 'g', long)]
    pub genomes: PathBuf,

    /// Directory of GFF files for annotation-based classification (optional)
    #[arg(short = 'a', long)]
    pub annotations: Option<PathBuf>,

    /// Output directory
    #[arg(short = 'o', long, default_value = "recore_output")]
    pub output: PathBuf,

    /// Minimum % identity for a hit
    #[arg(long, default_value_t = 90.0)]
    pub min_identity: f64,

    /// Minimum % query coverage for a hit (single hit or merged fragments)
    #[arg(long, default_value_t = 80.0)]
    pub min_coverage: f64,

    /// Number of threads for parallel genome processing
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,

    /// minimap2 preset: asm5 (≤5% divergence, faster), asm20 (≤20%, safer)
    #[arg(long, default_value = "asm5")]
    pub preset: String,

    /// Minimum cluster prevalence (fraction of genomes) to attempt rescue.
    /// Clusters present in fewer genomes than this are likely truly accessory.
    /// Set to 0 to try rescuing everything.
    #[arg(long, default_value_t = 0.1)]
    pub min_prevalence: f64,

    /// Use all alleles per cluster as queries (slower but more sensitive).
    /// By default, only the longest representative is used.
    #[arg(long, default_value_t = false)]
    pub multi_allele: bool,

    /// Check flanking gene synteny to validate rescues (requires -a).
    /// Rescued genes whose flanking context doesn't match are downgraded.
    #[arg(long, default_value_t = false)]
    pub synteny: bool,

    /// Write per-genome GFF3 files with Recore rescue features appended.
    /// Original GFF records are preserved when --annotations is provided.
    #[arg(long, default_value_t = false)]
    pub write_corrected_gff: bool,

    /// Policy for inserting accepted rescues into the binary corrected matrix.
    #[arg(long, value_enum, default_value_t = MatrixPolicy::Balanced)]
    pub matrix_policy: MatrixPolicy,

    /// Genome column to audit as the reference/background genome.
    /// Accepts exact matrix names or embedded AS identifiers such as AS00146.
    #[arg(long)]
    pub reference_genome: Option<String>,

    // ── Nucleotide Similarity Clustering ──
    /// Enable nucleotide similarity clustering to merge redundant Panaroo clusters.
    /// Detects clusters split at the AA level but identical at the DNA level.
    #[arg(long, default_value_t = false)]
    pub nscluster: bool,

    /// k-mer size for nucleotide similarity sketching.
    #[arg(long, default_value_t = 21)]
    pub nsc_kmer_size: usize,

    /// Max containment threshold for merging clusters (0.0-1.0).
    /// Higher = stricter (fewer merges). PanQC default: 0.8.
    #[arg(long, default_value_t = 0.8)]
    pub nsc_threshold: f64,

    /// NSC engine: "minhash" (built-in, no deps) or "mmseqs2" (exact, requires mmseqs2 in PATH).
    /// mmseqs2 linclust runs in linear time with alignment-based identity.
    #[arg(long, default_value = "minhash")]
    pub nsc_engine: String,

    /// Verbose debug logging: show sample matches, GFF parsing details, synteny decisions.
    #[arg(short = 'v', long, default_value_t = false)]
    pub verbose: bool,

    /// Minimum sequence identity for mmseqs2 linclust clustering (0.0-1.0).
    #[arg(long, default_value_t = 0.9)]
    pub nsc_min_seq_id: f64,

    /// Minimum coverage for mmseqs2 linclust (0.0-1.0).
    /// Coverage mode: bidirectional (shorter sequence must be covered).
    #[arg(long, default_value_t = 0.8)]
    pub nsc_coverage: f64,
}
