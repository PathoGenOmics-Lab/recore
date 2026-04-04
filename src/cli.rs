use clap::Parser;
use std::path::PathBuf;

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
}
