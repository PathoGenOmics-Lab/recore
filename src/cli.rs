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

    /// Directory of GFF files for validation (optional)
    #[arg(short = 'a', long)]
    pub annotations: Option<PathBuf>,

    /// Output directory
    #[arg(short = 'o', long, default_value = "recore_output")]
    pub output: PathBuf,

    /// Minimum % identity for a hit
    #[arg(long, default_value_t = 90.0)]
    pub min_identity: f64,

    /// Minimum % query coverage for a hit
    #[arg(long, default_value_t = 80.0)]
    pub min_coverage: f64,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,
}
