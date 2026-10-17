#!/usr/bin/env python3
"""Generate the small dataset the documentation tutorial walks through.

A tutorial that cannot be run is a description. This writes a complete, tiny
input set (matrix, pan-reference FASTA, assemblies, GFF3 annotations and a
Newick tree) with known problems planted in it, so every claim the tutorial
makes about the output is reproducible from the repository alone:

  false absence      a full-length gene sits in the assembly, absent in the
                     matrix. Recore should rescue it (DROPPED).
  missing CDS        the gene is in the assembly but the GFF has no CDS there.
  pseudogene         an internal stop codon truncates the reading frame.
  paralog trap       a diverged copy in the one genome that lacks the real
                     gene, in the wrong place on the contig. It must not become
                     a presence call.
  true absence       the gene really is not there, and must stay absent.
  clade acquisition  a gene only the L2 genomes carry.

Genome names follow the MTBC convention the report understands (L4_AS00101), so
the viewer can colour tips by lineage and root on the L8 outgroup.

Usage:

    python3 benchmarks/make_tutorial_dataset.py --out tutorial_data
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

METADATA_COLS = [
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Order",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
]

# Lineage -> genomes. L8 is the outgroup the report roots on by default.
LINEAGES = {
    "L8": ["L8_AS00801"],
    "L1": ["L1_AS00101", "L1_AS00102"],
    "L2": ["L2_AS00201", "L2_AS00202", "L2_AS00203"],
    "L3": ["L3_AS00301"],
    "L4": ["L4_AS00401", "L4_AS00402", "L4_AS00403", "L4_AS00404"],
    "L5": ["L5_AS00501"],
}
GENOMES = [g for group in LINEAGES.values() for g in group]
L2 = set(LINEAGES["L2"])
L4 = set(LINEAGES["L4"])

# One entry per cluster, in genome order along the contig so the flanking-gene
# context of a cluster is stable and synteny has something to validate against.
#   name        Rv-style locus name used as the cluster id
#   product     annotation text
#   carriers    genomes that genuinely carry the gene
#   note        which planted problem, if any, this cluster demonstrates
CLUSTERS = [
    ("Rv0001", "chromosomal replication initiator DnaA", set(GENOMES), None),
    ("Rv0002", "DNA polymerase III subunit beta", set(GENOMES), None),
    ("Rv0003", "DNA recombination protein RecF", set(GENOMES), "false_absence"),
    ("Rv0004", "hypothetical protein", set(GENOMES), None),
    ("Rv0005", "DNA gyrase subunit B", set(GENOMES), None),
    ("Rv0006", "DNA gyrase subunit A", set(GENOMES), "missing_cds"),
    ("Rv0007", "conserved membrane protein", set(GENOMES), None),
    ("Rv0008c", "conserved membrane protein", set(GENOMES) - {"L5_AS00501"}, None),
    ("Rv0009", "peptidyl-prolyl cis-trans isomerase A", set(GENOMES), "pseudogene"),
    ("Rv0010c", "conserved hypothetical protein", set(GENOMES), None),
    ("Rv0011c", "septation inhibitor protein", set(GENOMES), None),
    ("Rv0012", "conserved hypothetical protein", set(GENOMES), "false_absence"),
    ("Rv0013", "transcription elongation factor GreA", set(GENOMES), None),
    ("Rv0014c", "transmembrane serine/threonine kinase PknB", set(GENOMES), None),
    ("Rv0015c", "transmembrane serine/threonine kinase PknA", set(GENOMES), None),
    ("Rv0016c", "peptidoglycan synthase PbpA", set(GENOMES), None),
    ("Rv0017c", "cell division protein RodA", set(GENOMES), None),
    ("Rv0018c", "protein phosphatase PstP", set(GENOMES) - {"L4_AS00404"}, "paralog_trap"),
    ("Rv0019c", "FHA-domain protein", set(GENOMES), None),
    ("Rv0020c", "FHA-domain protein FhaA", set(GENOMES), None),
    ("Rv1180", "polyketide synthase Pks3", L2, "clade_acquisition"),
    ("Rv1181", "polyketide synthase Pks4", L2, None),
    ("Rv2930", "fatty-acid-AMP ligase FadD26", L4, None),
    ("Rv3810", "complement-binding protein PE-PGRS", {"L8_AS00801"}, "true_absence"),
]

CODONS = [
    "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "AAT", "AAC", "GAT", "GAC",
    "TGT", "TGC", "CAA", "CAG", "GAA", "GAG", "GGT", "GGC", "CAT", "CAC",
    "ATT", "ATC", "ATA", "TTA", "TTG", "CTT", "CTC", "AAA", "AAG", "ATG",
    "TTT", "TTC", "CCT", "CCC", "TCT", "TCC", "ACT", "ACC", "TGG", "TAT",
    "TAC", "GTT", "GTC", "GTA", "GTG",
]


def coding_sequence(rng: random.Random, codons: int) -> str:
    """An ORF with a start, no internal stops, and one terminal stop."""
    body = "".join(rng.choice(CODONS) for _ in range(codons))
    return "ATG" + body + "TAA"


STOPS = {"TAA", "TAG", "TGA"}


def diverge(rng: random.Random, seq: str, identity: float) -> str:
    """Mutate a sequence down to roughly `identity`, keeping the frame open.

    Mutations that would create an in-frame stop are rejected: a stop that
    arrives by accident turns a planted false absence into a planted
    pseudogene, and the tutorial would then be explaining the wrong thing.
    """
    bases = list(seq)
    changes = int(len(bases) * (1.0 - identity))
    for _ in range(changes):
        i = rng.randrange(3, len(bases) - 3)
        codon_at = (i // 3) * 3
        original = bases[i]
        for base in rng.sample([b for b in "ACGT" if b != original], 3):
            bases[i] = base
            if "".join(bases[codon_at:codon_at + 3]) not in STOPS:
                break
            bases[i] = original
    return "".join(bases)


def with_internal_stop(seq: str) -> str:
    """Plant an in-frame stop about 40% of the way in."""
    at = (len(seq) // 3) // 5 * 2 * 3
    return seq[:at] + "TAA" + seq[at + 3:]


def spacer(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def wrap(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def newick() -> str:
    """A tree whose shape matches the lineage grouping, L8 as the outgroup."""

    def tips(lineage: str, length: float) -> str:
        names = LINEAGES[lineage]
        if len(names) == 1:
            return f"{names[0]}:{length:.4f}"
        inner = ",".join(f"{n}:{length:.4f}" for n in names)
        return f"({inner}):{length * 2:.4f}"

    modern = f"({tips('L1', 0.0020)},{tips('L2', 0.0018)},{tips('L3', 0.0022)},{tips('L4', 0.0016)}):0.0040"
    animal = tips("L5", 0.0035)
    return f"(({modern},{animal}):0.0060,{tips('L8', 0.0090)});\n"


def build(out: Path, seed: int) -> None:
    rng = random.Random(seed)
    genomes_dir = out / "genomes"
    gff_dir = out / "annotations"
    for directory in (genomes_dir, gff_dir):
        directory.mkdir(parents=True, exist_ok=True)

    # Representative sequence per cluster. Lengths vary so the report's gene
    # arrows are not all the same width.
    reference = {
        name: coding_sequence(rng, rng.randrange(180, 460))
        for name, _product, _carriers, _note in CLUSTERS
    }

    with (out / "pan_genome_reference.fa").open("w") as fh:
        for name, seq in reference.items():
            fh.write(f">{name}\n{wrap(seq)}\n")

    # Lay every genome out as one contig, genes in cluster order, so the
    # flanking context of a cluster is the same everywhere it occurs.
    locus_of: dict[tuple[str, str], str] = {}
    matrix_absent: set[tuple[str, str]] = set()

    for genome in GENOMES:
        contig: list[str] = [spacer(rng, 400)]
        features: list[tuple[int, int, str, str, str]] = []
        index = 0

        def place(seq: str, name: str, product: str, annotate: bool = True) -> str:
            """Append one gene to the contig and return its locus tag."""
            nonlocal index
            index += 1
            locus = f"{genome}_{index:05d}"
            start = sum(len(chunk) for chunk in contig) + 1
            contig.append(seq)
            if annotate:
                features.append((start, start + len(seq) - 1, locus, name, product))
            contig.append(spacer(rng, rng.randrange(120, 320)))
            return locus

        for name, product, carriers, note in CLUSTERS:
            # The paralog trap lives in the one genome that does not carry the
            # real gene: a diverged copy is only a trap when there is a gap for
            # Recore to be tempted to fill with it.
            #
            # 85%, well under the default --min-identity of 90. It was 96%, and
            # it only stayed out of the matrix because identity used to be
            # chain-approximate and read it in the 80s. With exact identity a
            # 96% full-length annotated copy at the orthologous locus is,
            # by Recore's own rules, a DROPPED rescue, which is not a trap.
            if note == "paralog_trap" and genome not in carriers:
                contig.append(spacer(rng, 900))
                place(diverge(rng, reference[name], 0.85), f"{name}_paralog",
                      f"{product} (paralogous copy)")
                continue
            if genome not in carriers:
                continue

            seq = diverge(rng, reference[name], 0.995)
            annotate = True
            if note == "pseudogene" and genome in L2:
                seq = with_internal_stop(seq)
            if note == "missing_cds" and genome in {"L1_AS00101", "L3_AS00301"}:
                # In the assembly, invisible to the annotation.
                annotate = False

            locus = place(seq, name, product, annotate)
            if annotate:
                locus_of[(name, genome)] = locus

        contig.append(spacer(rng, 400))
        sequence = "".join(contig)

        with (genomes_dir / f"{genome}.fasta").open("w") as fh:
            fh.write(f">{genome}_contig1\n{wrap(sequence)}\n")

        with (gff_dir / f"{genome}.gff").open("w") as fh:
            fh.write("##gff-version 3\n")
            fh.write(f"##sequence-region {genome}_contig1 1 {len(sequence)}\n")
            for start, end, locus, name, product in features:
                strand = "-" if name.rstrip("_paralog").endswith("c") else "+"
                attrs = (
                    f"ID={locus};locus_tag={locus};gene={name};"
                    f"product={product.replace(';', ',')}"
                )
                fh.write(
                    f"{genome}_contig1\ttutorial\tCDS\t{start}\t{end}\t.\t"
                    f"{strand}\t0\t{attrs}\n"
                )

    # The matrix: what the upstream pangenome pipeline thought it saw. The
    # planted false absences are the cells it got wrong.
    for name, _product, carriers, note in CLUSTERS:
        if note == "false_absence":
            matrix_absent.update((name, g) for g in sorted(carriers)[:3])
        if note == "missing_cds":
            matrix_absent.update((name, g) for g in ("L1_AS00101", "L3_AS00301"))
        if note == "pseudogene":
            matrix_absent.update((name, g) for g in sorted(L2))

    header = METADATA_COLS + GENOMES
    with (out / "gene_presence_absence.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for name, product, carriers, _note in CLUSTERS:
            cells = []
            for genome in GENOMES:
                in_matrix = (
                    genome in carriers
                    and (name, genome) not in matrix_absent
                    and (name, genome) in locus_of
                )
                cells.append(locus_of[(name, genome)] if in_matrix else "")
            present = sum(1 for cell in cells if cell)
            length = len(reference[name])
            meta = [
                name, name, product, str(present), str(present), "1.0",
                "0", "0", "0", "", str(length), str(length), str(length),
            ]
            writer.writerow(meta + cells)

    (out / "species.tree").write_text(newick())

    planted = sum(
        1 for name, _p, carriers, _n in CLUSTERS
        for g in carriers if (name, g) in matrix_absent
    )
    print(f"genomes            : {len(GENOMES)}")
    print(f"clusters           : {len(CLUSTERS)}")
    print(f"planted absences   : {planted}")
    print(f"written to         : {out}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default=str(ROOT / "tutorial_data"),
        help="Directory to write the dataset into.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=11,
        help="Random seed, so the dataset is the same every time.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build(Path(args.out), args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
