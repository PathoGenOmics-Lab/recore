/// Integration test: mini dataset through the full pipeline.
///
/// Creates a tiny pangenome (3 clusters, 3 genomes) with 2 gaps,
/// runs the pipeline, and verifies all output files are generated correctly.
use std::fs;
use std::io::Write;
use std::path::Path;
use tempfile::tempdir;

fn write_file(dir: &Path, name: &str, content: &str) {
    let path = dir.join(name);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).unwrap();
    }
    let mut f = fs::File::create(&path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

fn minimap2_available() -> bool {
    std::process::Command::new("minimap2")
        .arg("--version")
        .output()
        .is_ok()
}

#[test]
fn test_e2e_mini_pipeline() {
    if !minimap2_available() {
        // Skipping keeps the suite usable on a machine without minimap2, but a
        // silent skip in CI means this test has never actually run there. CI sets
        // RECORE_REQUIRE_MINIMAP2 so a missing dependency fails loudly instead.
        assert!(
            std::env::var_os("RECORE_REQUIRE_MINIMAP2").is_none(),
            "RECORE_REQUIRE_MINIMAP2 is set but minimap2 is not in PATH"
        );
        eprintln!("SKIPPED: minimap2 not in PATH");
        return;
    }
    let tmp = tempdir().unwrap();
    let base = tmp.path();

    // ── Presence/absence CSV ──
    // 3 clusters (gA, gB, gC), 3 genomes (s1, s2, s3)
    // gA: present in all (core)
    // gB: present in s1 only → gaps in s2, s3
    // gC: present in s1, s2 → gap in s3
    write_file(base, "pa.csv", "\
Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,Accessory Order,QC,Min group size nuc,Max group size nuc,Avg group size nuc,s1,s2,s3
gA,gA,product_A,3,3,1.0,0,0,0,,100,100,100,gA_s1,gA_s2,gA_s3
gB,gB,product_B,1,1,1.0,0,0,0,,90,90,90,gB_s1,,
gC,gC,product_C,2,2,1.0,0,0,0,,80,80,80,gC_s1,gC_s2,
");

    // ── Pan genome reference ──
    // gB sequence: 100bp gene
    let gene_b = "ATGAAAGCGTTTGCGAAAGCGTTTGCGAAAGCGTTTGCGAAAGCGTTTGCGAAAGCGTTTGCGAAAGCGTTTGCGAAAGCGTTTGCGATCGTAA";
    // What s2 and s3 actually carry of gB: its first 78 of 94 bp, 83%, between
    // --min-coverage (80) and the adaptive floor for a short cluster (85). Both
    // used to carry all of gB and the test still expected it held back at
    // "80.9% coverage", which was chain-approximate alignment misreading a
    // repetitive sequence; aligned exactly, an intact gB is 100%.
    let gene_b_partial = &gene_b[..78];
    // gC sequence: 80bp gene
    let gene_c =
        "ATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTGA";

    write_file(base, "pan_ref.fa", &format!(
        ">gA\nATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAA\n\
         >gB\n{}\n\
         >gC\n{}\n", gene_b, gene_c));
    write_file(base, "tree.nwk", "((s1,s2),s3);\n");

    // ── Genome FASTAs ──
    // s2 and s3 contain the gene sequences (they should be rescued)
    let genomes_dir = base.join("genomes");
    fs::create_dir_all(&genomes_dir).unwrap();

    // s1: has all genes
    write_file(&genomes_dir, "s1.fasta", &format!(
        ">contig1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA{}TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT{}CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC{}GGGGGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "ATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAA", // gA
        gene_b, gene_c));

    // s2: has gA and gC, but also has gB hidden in the assembly
    write_file(&genomes_dir, "s2.fasta", &format!(
        ">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA{}TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT{}CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC{}GGGGGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "ATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAA",
        gene_b_partial,  // part of gB present in assembly, marked absent in CSV
        gene_c));

    // s3: has gA, and also has gB and gC hidden
    write_file(&genomes_dir, "s3.fasta", &format!(
        ">scaffold_1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA{}TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT{}CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC{}GGGGGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "ATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAA",
        gene_b_partial,  // part of gB present
        gene_c)); // gC present

    let output_dir = base.join("output");

    // ── Run Recore ──
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_recore"))
        .args([
            "-p",
            base.join("pa.csv").to_str().unwrap(),
            "-r",
            base.join("pan_ref.fa").to_str().unwrap(),
            "-g",
            genomes_dir.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-t",
            "1",
            "--min-prevalence",
            "0.0",
            "--write-corrected-gff",
            "--write-html-report",
            "--tree",
            base.join("tree.nwk").to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run recore binary");

    assert!(status.success(), "recore exited with non-zero status");

    // ── Verify outputs exist ──
    assert!(output_dir.join("gene_presence_absence_recore.csv").exists());
    assert!(output_dir
        .join("gene_presence_absence_recore_states.csv")
        .exists());
    assert!(output_dir.join("rescue_report.tsv").exists());
    assert!(output_dir.join("cluster_relationships.tsv").exists());
    assert!(output_dir.join("evolutionary_events.tsv").exists());
    assert!(output_dir.join("summary.txt").exists());
    assert!(output_dir.join("run_manifest.json").exists());
    assert!(output_dir.join("recore_report.html").exists());
    // gC in s3 is the only evidence accepted into the binary matrix, so it is the
    // only genome that gets a corrected GFF. The partial gB clears --min-coverage
    // at 83% but not the adaptive floor for a short, rare cluster, so it stays
    // audit-only.
    assert!(output_dir
        .join("corrected_gff")
        .join("s3.recore.gff3")
        .exists());
    assert!(
        !output_dir
            .join("corrected_gff")
            .join("s2.recore.gff3")
            .exists(),
        "s2 has only audit-only evidence, so it must not get a corrected GFF"
    );

    // ── Verify rescue report ──
    let report = fs::read_to_string(output_dir.join("rescue_report.tsv")).unwrap();
    let lines: Vec<&str> = report.lines().collect();
    assert!(
        lines.len() >= 2,
        "rescue_report should have header + at least 1 rescue"
    );
    assert!(
        lines[0].contains("gene_cluster"),
        "header should contain gene_cluster"
    );

    // gB should be rescued in s2 and s3
    let gb_rescues: Vec<&&str> = lines.iter().filter(|l| l.starts_with("gB\t")).collect();
    assert!(
        gb_rescues.len() >= 2,
        "gB should be rescued in s2 and s3, got {}",
        gb_rescues.len()
    );

    // ── Verify corrected CSV ──
    let csv_content =
        fs::read_to_string(output_dir.join("gene_presence_absence_recore.csv")).unwrap();
    let csv_lines: Vec<&str> = csv_content.lines().collect();
    assert_eq!(csv_lines.len(), 4, "header + 3 clusters");

    // The whole point of the design: accepted evidence reaches the binary matrix,
    // guarded evidence does not, and both stay visible in the multistate matrix.
    let gc_row = csv_lines.iter().find(|l| l.starts_with("gC")).unwrap();
    assert!(
        gc_row.contains("RESCUED:"),
        "gC/s3 clears the adaptive floor and should be corrected: {}",
        gc_row
    );
    let gb_row = csv_lines.iter().find(|l| l.starts_with("gB")).unwrap();
    assert!(
        !gb_row.contains("RESCUED:"),
        "gB is below the adaptive coverage floor and must stay audit-only: {}",
        gb_row
    );

    let state_content =
        fs::read_to_string(output_dir.join("gene_presence_absence_recore_states.csv")).unwrap();
    assert!(state_content.contains("PRESENT_INPUT"));
    assert!(state_content.contains("SEQUENCE_ONLY"));

    // ... and the audit table has to say why gB was held back.
    let evidence = fs::read_to_string(output_dir.join("orthology_evidence.tsv")).unwrap();
    // gene_cluster is not the first column here, so match the field rather than
    // the line prefix.
    let gb_evidence: Vec<&str> = evidence
        .lines()
        .skip(1)
        .filter(|line| line.split('\t').any(|field| field == "gB"))
        .collect();
    assert_eq!(gb_evidence.len(), 2, "gB has evidence rows for s2 and s3");
    for row in &gb_evidence {
        assert!(
            row.contains("audit_only_quality"),
            "gB evidence should be audit-only: {row}"
        );
        assert!(
            row.contains("adaptive_sequence_below_floor"),
            "gB should carry the adaptive-floor risk flag: {row}"
        );
    }

    // ── Verify summary ──
    let summary = fs::read_to_string(output_dir.join("summary.txt")).unwrap();
    assert!(summary.contains("Rescued:"));
    assert!(summary.contains("Core genome"));
    let manifest = fs::read_to_string(output_dir.join("run_manifest.json")).unwrap();
    assert!(manifest.contains("\"schema_version\": \"recore.outputs.v7\""));
    assert!(manifest.contains("\"created_utc\": \"20"));
    assert!(manifest.contains("\"command_line\": "));
    assert!(manifest.contains("\"matrix_policy\": \"Balanced\""));
    assert!(manifest.contains("\"write_html_report\": true"));
    assert!(manifest.contains("\"matrix_corrected_total\""));

    let html_report = fs::read_to_string(output_dir.join("recore_report.html")).unwrap();
    assert!(html_report.contains("<title>Recore run report</title>"));
    assert!(html_report.contains("Orthology actions"));
    assert!(html_report.contains("id='runDetails'"));
    assert!(html_report.contains("Unique syntenies"));
    assert!(html_report.contains("Top candidate clusters"));
    assert!(html_report.contains("Corrected matrix"));
    assert!(html_report.contains("window.recoreExplorerData"));
    assert!(html_report.contains("columnHeaderCanvas"));
    assert!(html_report.contains("tipColorMode"));
    assert!(html_report.contains("rowHeight"));
    assert!(html_report.contains("matrixColorMode"));
    assert!(html_report.contains("processLegend"));
    assert!(html_report.contains("lineageLegend"));
    assert!(html_report.contains("data-mechanism"));
    assert!(html_report.contains("targetReferenceGene"));
    assert!(html_report.contains("value='lineage'"));
    assert!(html_report.contains("phylogramMode"));
    assert!(html_report.contains("referenceOnly"));
    assert!(html_report.contains("visibleColumnRange"));
    // Find-anything palette, shortcut sheet and the per-pane menus that replaced
    // the shared control strip.
    assert!(html_report.contains("id='palette'"));
    assert!(html_report.contains("id='clusterOmnibox'"));
    assert!(html_report.contains("id='shortcutSheet'"));
    assert!(html_report.contains("id='treeMenu'"));
    assert!(html_report.contains("id='matrixMenu'"));
    assert!(html_report.contains("id='genomeDetails'"));
    assert!(html_report.contains("id='overviewCanvas'"));
    assert!(html_report.contains("id='treeScaleBar'"));
    // The mark ships web-ready from tools/make_logo_variants.py: cropped to the
    // artwork, and free of the export's fixed size and preserveAspectRatio="none",
    // which would let any CSS box stretch it out of shape.
    assert!(html_report.contains("viewBox=\"98 120 1857 854\""));
    assert!(!html_report.contains("preserveAspectRatio=\"none\""));
    assert!(!html_report.contains("sodipodi"));
    // Nothing in the mark is painted white any more, so it needs no plate and
    // reads as holes rather than blobs on a dark panel.
    assert!(!html_report.contains("fill=\"rgb(255,255,255)\""));
    // The report has to stay openable from an offline results directory. Links to
    // the project are fine; loading a resource is not.
    for tag in ["<link", "<script", "<img"] {
        let mut rest = html_report.as_str();
        while let Some(at) = rest.find(tag) {
            rest = &rest[at + tag.len()..];
            let end = rest.find('>').unwrap_or(rest.len());
            assert!(
                !rest[..end].contains("http://") && !rest[..end].contains("https://"),
                "report loads a remote resource in a {tag} element"
            );
        }
    }

    // ── Verify corrected GFF ──
    let gff_content =
        fs::read_to_string(output_dir.join("corrected_gff").join("s3.recore.gff3")).unwrap();
    assert!(gff_content.contains("##gff-version 3"));
    assert!(gff_content.contains("Recore\tCDS"));
    assert!(gff_content.contains("gene_cluster=gC"));
    assert!(
        !gff_content.contains("gene_cluster=gB"),
        "audit-only evidence must not be appended to a corrected GFF"
    );

    // ── Verify report is sorted (reproducibility) ──
    let data_lines: Vec<&str> = lines[1..].to_vec();
    if data_lines.len() >= 2 {
        for pair in data_lines.windows(2) {
            let a_cluster: &str = pair[0].split('\t').next().unwrap();
            let b_cluster: &str = pair[1].split('\t').next().unwrap();
            assert!(
                a_cluster <= b_cluster,
                "rescue_report should be sorted by cluster: {} > {}",
                a_cluster,
                b_cluster
            );
        }
    }
}

/// Tiny deterministic generator, so the scale test needs no extra dependency.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() % 4) as usize]
    }
    fn seq(&mut self, n: usize) -> String {
        (0..n).map(|_| self.base() as char).collect()
    }
    /// An open reading frame with no internal stop codon.
    fn orf(&mut self, codons: usize) -> String {
        let mut s = String::from("ATG");
        while s.len() < 3 * (codons + 1) {
            let c = self.seq(3);
            if !matches!(c.as_str(), "TAA" | "TAG" | "TGA") {
                s.push_str(&c);
            }
        }
        s + "TAA"
    }
    /// Change one base in each of `n` codons, never creating a stop codon.
    fn mutate(&mut self, s: &str, n: usize) -> String {
        let mut b = s.as_bytes().to_vec();
        let codons = b.len() / 3;
        let mut done = 0;
        while done < n {
            let c = 1 + (self.next() as usize % (codons - 2));
            let i = 3 * c + 2;
            let old = b[i];
            let new = self.base();
            if new == old {
                continue;
            }
            b[i] = new;
            if matches!(&b[3 * c..3 * c + 3], b"TAA" | b"TAG" | b"TGA") {
                b[i] = old;
                continue;
            }
            done += 1;
        }
        String::from_utf8(b).unwrap()
    }
}

/// The search used to put every genome into one minimap2 target, under one
/// `-N 50` secondary budget, and without base-level alignment. Once more than
/// about fifty genomes carried a gene, their copies used up the budget and the
/// gap genomes' hits were dropped: on 120 genomes almost none of the planted
/// false absences came back, with exit 0. And identity was chain-approximate,
/// reading these 99.5% copies in the low 90s.
#[test]
fn test_search_rescues_past_fifty_carriers_with_exact_identity() {
    if !minimap2_available() {
        assert!(
            std::env::var_os("RECORE_REQUIRE_MINIMAP2").is_none(),
            "RECORE_REQUIRE_MINIMAP2 is set but minimap2 is not in PATH"
        );
        eprintln!("SKIPPED: minimap2 not in PATH");
        return;
    }
    let tmp = tempdir().unwrap();
    let base = tmp.path();
    let mut rng = Lcg(7);

    let n = 120;
    let genomes: Vec<String> = (0..n).map(|i| format!("S{i:03}")).collect();
    let blank: Vec<usize> = (0..10).map(|k| k * n / 10 + 1).collect();
    // accA in even genomes and accB in odd ones, so every genome has a true gap
    // and lands in the search, which is what made the shared target crowded.
    let clusters = ["core1", "tgt", "accA", "accB"];
    let refs: Vec<String> = clusters.iter().map(|_| rng.orf(200)).collect();
    let carries = |c: &str, i: usize| match c {
        "accA" => i % 2 == 0,
        "accB" => i % 2 == 1,
        _ => true,
    };

    let mut fasta = String::new();
    for (name, seq) in clusters.iter().zip(&refs) {
        fasta.push_str(&format!(">{name}\n{seq}\n"));
    }
    write_file(base, "pan_ref.fa", &fasta);

    let mut csv = String::from(
        "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,\
         Avg sequences per isolate,Genome Fragment,Order within Fragment,\
         Accessory Order,QC,Min group size nuc,Max group size nuc,Avg group size nuc",
    );
    for g in &genomes {
        csv.push(',');
        csv.push_str(g);
    }
    csv.push('\n');
    for (ci, name) in clusters.iter().enumerate() {
        let cells: Vec<String> = (0..n)
            .map(|i| {
                let present = carries(name, i) && !(*name == "tgt" && blank.contains(&i));
                if present {
                    format!("{}_{ci}", genomes[i])
                } else {
                    String::new()
                }
            })
            .collect();
        let count = cells.iter().filter(|c| !c.is_empty()).count();
        csv.push_str(&format!(
            "{name},,{name} protein,{count},{count},1,0,0,0,,603,603,603,{}\n",
            cells.join(",")
        ));
    }
    write_file(base, "pa.csv", &csv);

    for (i, g) in genomes.iter().enumerate() {
        let mut contig = rng.seq(300);
        for (name, seq) in clusters.iter().zip(&refs) {
            if carries(name, i) {
                // 3 changes in 603 bp: a 99.5% copy, in every genome that
                // carries the gene, blanked in the matrix or not.
                contig.push_str(&rng.mutate(seq, 3));
                contig.push_str(&rng.seq(200));
            }
        }
        write_file(
            &base.join("genomes"),
            &format!("{g}.fasta"),
            &format!(">{g}_c1\n{contig}\n"),
        );
    }

    let output_dir = base.join("output");
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_recore"))
        .args([
            "-p",
            base.join("pa.csv").to_str().unwrap(),
            "-r",
            base.join("pan_ref.fa").to_str().unwrap(),
            "-g",
            base.join("genomes").to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-t",
            "4",
        ])
        .status()
        .expect("Failed to run recore binary");
    assert!(status.success(), "recore exited with non-zero status");

    let report = fs::read_to_string(output_dir.join("rescue_report.tsv")).unwrap();
    let header: Vec<&str> = report.lines().next().unwrap().split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let tgt: Vec<Vec<&str>> = report
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| f[col("gene_cluster")] == "tgt")
        .collect();
    assert_eq!(
        tgt.len(),
        10,
        "every planted false absence should be rescued, got {}",
        tgt.len()
    );
    for row in &tgt {
        let identity: f64 = row[col("identity")].parse().unwrap();
        assert!(
            identity >= 99.0,
            "a 99.5% copy must read as such, not chain-approximate: {identity}"
        );
    }
}

/// A gap the search never looked at used to be written as ABSENT, the word for
/// a gap that was searched and came back empty. Each of the three ways a gap is
/// left out of the search has to say so, and a searched empty gap still reads
/// ABSENT.
#[test]
fn test_unsearched_gaps_say_why_instead_of_absent() {
    if !minimap2_available() {
        assert!(
            std::env::var_os("RECORE_REQUIRE_MINIMAP2").is_none(),
            "RECORE_REQUIRE_MINIMAP2 is set but minimap2 is not in PATH"
        );
        eprintln!("SKIPPED: minimap2 not in PATH");
        return;
    }
    let tmp = tempdir().unwrap();
    let base = tmp.path();
    let mut rng = Lcg(11);

    // At --min-prevalence 0.5: gB (1 of 4) is below the floor, gD (2 of 4)
    // passes it but has no sequence to search with, and gE's only gap is in s4,
    // which has no FASTA. gC's gap in s3 is searched and s3 lacks the gene.
    let genomes = ["s1", "s2", "s3", "s4"];
    let clusters = ["gA", "gB", "gC", "gD", "gE"];
    let present = |c: &str, g: &str| match c {
        "gB" => g == "s1",
        "gC" => g != "s3",
        "gD" => g == "s1" || g == "s2",
        "gE" => g != "s4",
        _ => true,
    };
    let seqs: Vec<String> = clusters.iter().map(|_| rng.orf(150)).collect();
    let mut fasta = String::new();
    for (name, seq) in clusters.iter().zip(&seqs) {
        if *name != "gD" {
            fasta.push_str(&format!(">{name}\n{seq}\n"));
        }
    }
    write_file(base, "pan_ref.fa", &fasta);

    let mut csv = String::from(
        "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,\
         Avg sequences per isolate,Genome Fragment,Order within Fragment,\
         Accessory Order,QC,Min group size nuc,Max group size nuc,Avg group size nuc,s1,s2,s3,s4\n",
    );
    for name in clusters {
        let cells: Vec<String> = genomes
            .iter()
            .map(|g| {
                if present(name, g) {
                    format!("{name}_{g}")
                } else {
                    String::new()
                }
            })
            .collect();
        let n = cells.iter().filter(|c| !c.is_empty()).count();
        csv.push_str(&format!(
            "{name},{name},product,{n},{n},1.0,0,0,0,,450,450,450,{}\n",
            cells.join(",")
        ));
    }
    write_file(base, "pa.csv", &csv);
    write_file(base, "tree.nwk", "(((s1,s2),s3),s4);\n");

    // s2 carries gB although the matrix says it does not, so had gB been
    // searched it would have been rescued: the state below is not a search
    // that found nothing.
    let genomes_dir = base.join("genomes");
    for g in ["s1", "s2", "s3"] {
        let mut contig = rng.seq(300);
        for (name, seq) in clusters.iter().zip(&seqs) {
            if present(name, g) || (*name == "gB" && g == "s2") {
                contig.push_str(seq);
                contig.push_str(&rng.seq(300));
            }
        }
        write_file(
            &genomes_dir,
            &format!("{g}.fasta"),
            &format!(">{g}_1\n{contig}\n"),
        );
    }

    let output_dir = base.join("output");
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_recore"))
        .args([
            "-p",
            base.join("pa.csv").to_str().unwrap(),
            "-r",
            base.join("pan_ref.fa").to_str().unwrap(),
            "-g",
            genomes_dir.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-t",
            "1",
            "--min-prevalence",
            "0.5",
            "--write-html-report",
            "--tree",
            base.join("tree.nwk").to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run recore binary");
    assert!(status.success(), "recore exited with non-zero status");

    let states =
        fs::read_to_string(output_dir.join("gene_presence_absence_recore_states.csv")).unwrap();
    let mut lines = states.lines();
    let header: Vec<&str> = lines.next().unwrap().split(',').collect();
    let state = |cluster: &str, genome: &str| -> String {
        let col = header.iter().position(|h| *h == genome).unwrap();
        let row = states
            .lines()
            .find(|l| l.split(',').next() == Some(cluster))
            .unwrap();
        row.split(',').nth(col).unwrap().to_string()
    };
    for g in ["s2", "s3", "s4"] {
        assert_eq!(state("gB", g), "NOT_SEARCHED_LOW_PREVALENCE", "gB/{g}");
    }
    for g in ["s3", "s4"] {
        assert_eq!(state("gD", g), "NOT_SEARCHED_NO_CLUSTER_SEQUENCE", "gD/{g}");
    }
    assert_eq!(state("gE", "s4"), "NOT_SEARCHED_NO_FASTA");
    assert_eq!(
        state("gC", "s3"),
        "ABSENT",
        "a searched empty gap stays ABSENT"
    );
    assert_eq!(state("gC", "s1"), "PRESENT_INPUT");

    // The viewer gets the same reasons: one code per cluster and the genomes
    // without a FASTA, and the run details say how many gaps were not searched.
    let html = fs::read_to_string(output_dir.join("recore_report.html")).unwrap();
    assert!(
        html.contains(
            "\"notSearched\":{\"minPrevalence\":0.5,\"clusters\":\"01020\",\"genomesWithoutFasta\":[3]}"
        ),
        "report payload is missing the not-searched cells"
    );
    // Per genome: present, corrected, NSC, kept out, absent, not searched.
    assert!(
        html.contains("\"genomeStats\":[[5,0,0,0,0,0],[4,0,0,0,0,1],[2,0,0,0,1,2],[2,0,0,0,0,3]]"),
        "genome totals must count not-searched cells apart from absent ones"
    );
    assert!(html.contains("6 gaps were never searched (3 in clusters below the 50% prevalence floor, 2 in clusters with no sequence in the pan-genome reference, 1 in genomes with no FASTA file)"));

    // The evolutionary audit counts the absences it cannot vouch for, says why,
    // and does not present a call made of them as an audited one.
    let events = fs::read_to_string(output_dir.join("evolutionary_events.tsv")).unwrap();
    let mut rows = events.lines();
    let header: Vec<&str> = rows.next().unwrap().split('\t').collect();
    let column = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let event = |cluster: &str| -> Vec<String> {
        events
            .lines()
            .find(|l| l.split('\t').next() == Some(cluster))
            .unwrap()
            .split('\t')
            .map(str::to_string)
            .collect()
    };
    let (count, reason) = (column("absent_not_searched"), column("not_searched_reason"));
    for (cluster, n, why) in [
        ("gA", "0", ""),
        ("gB", "3", "low_prevalence"),
        ("gC", "0", ""),
        ("gD", "2", "no_cluster_sequence"),
        ("gE", "1", "no_fasta"),
    ] {
        let row = event(cluster);
        assert_eq!(
            (row[count].as_str(), row[reason].as_str()),
            (n, why),
            "{cluster}"
        );
    }
    let gb = event("gB");
    assert!(gb[column("mechanism_detail")].ends_with(";absences_not_searched"));
    assert!(gb[column("recommended_wording")]
        .contains("None of the 3 absences behind this call were searched"));
    let gc = event("gC");
    assert!(!gc[column("mechanism_detail")].contains("absences_not_searched"));
    assert!(!gc[column("recommended_wording")].contains("never searched"));
}
