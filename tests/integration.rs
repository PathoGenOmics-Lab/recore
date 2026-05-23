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
        gene_b,  // gB present in assembly but marked absent in CSV
        gene_c));

    // s3: has gA, and also has gB and gC hidden
    write_file(&genomes_dir, "s3.fasta", &format!(
        ">scaffold_1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA{}TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT{}CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC{}GGGGGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "ATGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAA",
        gene_b,  // gB present
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
    assert!(output_dir
        .join("corrected_gff")
        .join("s2.recore.gff3")
        .exists());

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

    // Check that rescued cells contain "RESCUED:"
    let gb_row = csv_lines.iter().find(|l| l.starts_with("gB")).unwrap();
    assert!(
        gb_row.contains("RESCUED:"),
        "gB row should have RESCUED cells: {}",
        gb_row
    );

    let state_content =
        fs::read_to_string(output_dir.join("gene_presence_absence_recore_states.csv")).unwrap();
    assert!(state_content.contains("PRESENT_INPUT"));
    assert!(state_content.contains("SEQUENCE_ONLY"));

    // ── Verify summary ──
    let summary = fs::read_to_string(output_dir.join("summary.txt")).unwrap();
    assert!(summary.contains("Rescued:"));
    assert!(summary.contains("Core genome"));
    let manifest = fs::read_to_string(output_dir.join("run_manifest.json")).unwrap();
    assert!(manifest.contains("\"schema_version\": \"recore.outputs.v6\""));
    assert!(manifest.contains("\"matrix_policy\": \"Balanced\""));
    assert!(manifest.contains("\"write_html_report\": true"));
    assert!(manifest.contains("\"matrix_corrected_total\""));

    let html_report = fs::read_to_string(output_dir.join("recore_report.html")).unwrap();
    assert!(html_report.contains("<title>Recore run report</title>"));
    assert!(html_report.contains("Orthology actions"));
    assert!(html_report.contains("Top candidate clusters"));
    assert!(html_report.contains("Phylogeny and corrected matrix explorer"));
    assert!(html_report.contains("window.recoreExplorerData"));
    assert!(html_report.contains("columnHeaderCanvas"));
    assert!(html_report.contains("tipColorMode"));
    assert!(html_report.contains("rowHeight"));
    assert!(html_report.contains("matrixColorMode"));
    assert!(html_report.contains("processLegend"));
    assert!(html_report.contains("lineageLegend"));
    assert!(html_report.contains("data-mechanism"));
    assert!(html_report.contains("targetReferenceGene"));
    assert!(html_report.contains("Color tips by lineage"));
    assert!(html_report.contains("phylogramMode"));
    assert!(html_report.contains("referenceOnly"));
    assert!(html_report.contains("rebuildBaseMatrix"));

    // ── Verify corrected GFF ──
    let gff_content =
        fs::read_to_string(output_dir.join("corrected_gff").join("s2.recore.gff3")).unwrap();
    assert!(gff_content.contains("Recore\tCDS"));
    assert!(gff_content.contains("gene_cluster=gB"));

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
