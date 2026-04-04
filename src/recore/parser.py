"""Parse Panaroo/Roary output files."""

import csv
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

log = logging.getLogger("recore")


def parse_presence_absence(csv_path: Path) -> pd.DataFrame:
    """Parse gene_presence_absence.csv (Panaroo/Roary format).

    Returns a DataFrame where:
    - Index: gene cluster name
    - Columns: genome names
    - Values: gene ID if present, NaN if absent
    """
    # Detect format: Panaroo uses comma, Roary uses comma but with extra metadata cols
    df = pd.read_csv(csv_path, low_memory=False)

    # Panaroo format: Gene, Non-unique Gene name, Annotation, No. isolates,
    #                 No. sequences, Avg sequences per isolate, Genome Fragment,
    #                 Order within Fragment, Accessory Fragment, Accessory Order,
    #                 QC, Min group size nuc, Max group size nuc, Avg group size nuc,
    #                 then genome columns
    # Roary format: similar but slightly different column names

    # Find genome columns: everything after the metadata columns
    metadata_cols = {
        "Gene",
        "Non-unique Gene name",
        "Annotation",
        "No. isolates",
        "No. sequences",
        "Avg sequences per isolate",
        "Genome Fragment",
        "Order within Fragment",
        "Accessory Fragment",
        "Accessory Order with Fragment",  # Roary
        "Accessory Order",                # Panaroo
        "QC",
        "Min group size nuc",
        "Max group size nuc",
        "Avg group size nuc",
    }

    genome_cols = [c for c in df.columns if c not in metadata_cols]

    if not genome_cols:
        raise ValueError(
            f"No genome columns found in {csv_path}. "
            f"Expected Panaroo or Roary gene_presence_absence.csv format."
        )

    log.info(
        f"Parsed presence/absence matrix: {len(df)} clusters × {len(genome_cols)} genomes"
    )

    # Build presence/absence matrix
    pa = df.set_index("Gene")[genome_cols]

    return pa


def find_gaps(pa: pd.DataFrame) -> list[tuple[str, str]]:
    """Find all (gene, genome) pairs where the gene is absent.

    Returns list of (gene_cluster, genome_name) tuples.
    """
    gaps = []
    for gene in pa.index:
        for genome in pa.columns:
            val = pa.at[gene, genome]
            if pd.isna(val) or (isinstance(val, str) and val.strip() == ""):
                gaps.append((gene, genome))

    log.info(f"Found {len(gaps)} gaps in presence/absence matrix")
    return gaps


def parse_pan_reference(fasta_path: Path) -> dict[str, str]:
    """Parse pan_genome_reference.fa → {cluster_name: sequence}.

    Panaroo format: >cluster_name~~~gene_id description
    """
    refs = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        # Panaroo names: cluster_name or cluster_name~~~gene_id
        name = record.id.split("~~~")[0] if "~~~" in record.id else record.id
        # Keep the longest sequence per cluster
        seq = str(record.seq)
        if name not in refs or len(seq) > len(refs[name]):
            refs[name] = seq

    log.info(f"Parsed {len(refs)} representative sequences from {fasta_path}")
    return refs


def find_genome_files(genomes_dir: Path) -> dict[str, Path]:
    """Map genome names to FASTA file paths.

    Tries common extensions: .fasta, .fa, .fna, .fasta.gz, .fa.gz, .fna.gz
    """
    extensions = [".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz"]
    genome_files = {}

    for ext in extensions:
        for f in genomes_dir.glob(f"*{ext}"):
            name = f.name
            for e in extensions:
                if name.endswith(e):
                    name = name[: -len(e)]
                    break
            if name not in genome_files:
                genome_files[name] = f

    log.info(f"Found {len(genome_files)} genome files in {genomes_dir}")
    return genome_files


def parse_gff_genes(gff_path: Path) -> list[dict]:
    """Parse CDS features from a GFF file.

    Returns list of dicts: {seqid, start, end, strand, gene_id, product}
    """
    genes = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "\t" not in line:
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "CDS":
                continue

            attrs = {}
            for attr in fields[8].split(";"):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attrs[key] = val

            genes.append({
                "seqid": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "gene_id": attrs.get("ID", attrs.get("locus_tag", "")),
                "product": attrs.get("product", ""),
            })

    return genes
