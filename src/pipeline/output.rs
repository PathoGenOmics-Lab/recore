use super::cellcodec;
use super::evolution::{PhyloNode, PhyloTree};
use super::reference_audit::resolve_reference_genome;
use super::*;

/// What the corrected matrix actually holds, counted as it is written.
///
/// The summary, the manifest and the report used to compute these separately,
/// and drifted: rescues in rows that NSC merged away were counted as matrix
/// corrections although those rows are not written, and "Remaining gaps" was
/// searched gaps minus corrections, which is not what is left in the file.
#[derive(Debug, Default, Clone, Copy)]
pub(super) struct CorrectedMatrixCounts {
    pub(super) rescued_cells: usize,
    pub(super) nsc_cells: usize,
    pub(super) absent_cells: usize,
}

pub(super) fn write_corrected_csv(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    output_dir: &std::path::Path,
    nsc_result: &Option<nscluster::NsClusterResult>,
) -> Result<CorrectedMatrixCounts> {
    let mut counts = CorrectedMatrixCounts::default();
    // Build set of clusters to skip (merged into another by NSC)
    let clusters_to_skip: ahash::HashSet<usize> = nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();

    // Build NSC presence updates: (rep_ci, gi) pairs to mark present
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let path = output_dir.join("gene_presence_absence_recore.csv");
    let file = fs::File::create(&path)?;
    let writer = BufWriter::with_capacity(256 * 1024, file);
    let mut csv_writer = csv::Writer::from_writer(writer);

    csv_writer.write_record(pa.header())?;

    let genome_idx_set: ahash::HashSet<usize> = pa.genome_col_indices().iter().copied().collect();
    let meta_col_indices: Vec<usize> = (0..pa.header().len())
        .filter(|i| !genome_idx_set.contains(i))
        .collect();

    let isolates_col = pa.header().iter().position(|h| h == "No. isolates");
    let sequences_col = pa.header().iter().position(|h| h == "No. sequences");
    let avg_seq_col = pa
        .header()
        .iter()
        .position(|h| h == "Avg sequences per isolate");

    let n_cols = pa.header().len();
    let mut row: Vec<String> = vec![String::new(); n_cols];
    let mut rescue_buf = String::new();

    for ci_val in 0..pa.n_clusters() {
        // Skip clusters merged by NSC (their presence is merged into representative)
        if clusters_to_skip.contains(&ci_val) {
            continue;
        }

        let ci = ClusterIdx(ci_val);

        for (mi, &col_idx) in meta_col_indices.iter().enumerate() {
            let meta = pa.metadata(ci);
            if mi < meta.len() {
                row[col_idx].clone_from(&meta[mi]);
            } else {
                row[col_idx].clear();
            }
        }

        let mut isolate_count = 0usize;
        let mut sequence_count = 0usize;
        for (gi_val, &col_idx) in pa.genome_col_indices().iter().enumerate() {
            let gi = GenomeIdx(gi_val);
            if let Some(r) = rescue_map.get(&(ci, gi)) {
                rescue_buf.clear();
                use std::fmt::Write as FmtWrite;
                // 1-based coords for the corrected CSV (biologist-friendly)
                write!(rescue_buf, "RESCUED:{}:{}-{}", r.contig, r.start + 1, r.end)
                    .expect("write to String is infallible");
                row[col_idx].clone_from(&rescue_buf);
                isolate_count += 1;
                sequence_count += 1;
                counts.rescued_cells += 1;
            } else if nsc_presence.contains(&(ci_val, gi_val)) {
                // Gene absent in this cluster but present in a merged NSC member
                row[col_idx].clear();
                row[col_idx].push_str("NSC_MERGED");
                isolate_count += 1;
                sequence_count += 1;
                counts.nsc_cells += 1;
            } else {
                row[col_idx].clear();
                let cell = pa.cell(ci, gi);
                row[col_idx].push_str(cell);
                if !cell.trim().is_empty() {
                    isolate_count += 1;
                    sequence_count += cell_sequence_count(cell);
                } else {
                    counts.absent_cells += 1;
                }
            }
        }

        use std::fmt::Write as FmtWrite;
        if let Some(col) = isolates_col {
            row[col].clear();
            write!(row[col], "{}", isolate_count).expect("write to String is infallible");
        }
        if let Some(col) = sequences_col {
            row[col].clear();
            write!(row[col], "{}", sequence_count).expect("write to String is infallible");
        }
        if let Some(col) = avg_seq_col {
            row[col].clear();
            let avg = if isolate_count > 0 {
                sequence_count as f64 / isolate_count as f64
            } else {
                0.0
            };
            write!(row[col], "{avg:.3}").expect("write to String is infallible");
        }

        csv_writer.write_record(row.iter().map(|s| s.as_str()))?;
    }
    csv_writer.flush()?;
    Ok(counts)
}

pub(super) fn clear_optional_outputs(output_dir: &std::path::Path) -> Result<()> {
    for filename in [
        "synteny_rejected_report.tsv",
        "paralog_rejected_report.tsv",
        "reference_genome_audit.tsv",
        "reference_cluster_events.tsv",
        "nscluster_report.tsv",
        "cluster_locus_context.tsv",
        "cluster_locus_summary.tsv",
        "cluster_synteny_blocks.tsv",
        "recore_report.html",
        "recore_report.md",
    ] {
        let path = output_dir.join(filename);
        match fs::remove_file(&path) {
            Ok(()) => {}
            Err(err) if err.kind() == std::io::ErrorKind::NotFound => {}
            Err(err) => {
                return Err(err).with_context(|| format!("Cannot remove {}", path.display()))
            }
        }
    }

    let corrected_gff = output_dir.join("corrected_gff");
    match fs::remove_dir_all(&corrected_gff) {
        Ok(()) => {}
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {}
        Err(err) => {
            return Err(err).with_context(|| format!("Cannot remove {}", corrected_gff.display()))
        }
    }
    Ok(())
}

fn cell_sequence_count(cell: &str) -> usize {
    cell.split(';')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .count()
        .max(usize::from(!cell.trim().is_empty()))
}

pub(super) struct MultistateMatrixInputs<'a> {
    pub(super) pa: &'a crate::parser::PresenceAbsence,
    pub(super) accepted_rescue_map: &'a HashMap<(ClusterIdx, GenomeIdx), &'a RescueEvent>,
    pub(super) nsc_rescue_map: &'a HashMap<(ClusterIdx, GenomeIdx), &'a RescueEvent>,
    pub(super) synteny_rejected: &'a [RescueEvent],
    pub(super) paralog_rejected: &'a [RescueEvent],
    pub(super) gene_to_cluster_maps: &'a HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    pub(super) output_dir: &'a std::path::Path,
    pub(super) nsc_result: &'a Option<nscluster::NsClusterResult>,
}

pub(super) fn write_multistate_matrix(inputs: MultistateMatrixInputs<'_>) -> Result<()> {
    let pa = inputs.pa;
    let clusters_to_skip: ahash::HashSet<usize> = inputs
        .nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();

    let nsc_presence = nsc_presence_updates(pa, inputs.nsc_result, inputs.nsc_rescue_map);

    let synteny_rejected_map = build_best_rescue_map(inputs.synteny_rejected.iter());
    let paralog_rejected_map = build_best_rescue_map(inputs.paralog_rejected.iter());

    let path = inputs
        .output_dir
        .join("gene_presence_absence_recore_states.csv");
    let file = fs::File::create(&path)?;
    let writer = BufWriter::with_capacity(256 * 1024, file);
    let mut csv_writer = csv::Writer::from_writer(writer);
    csv_writer.write_record(pa.header())?;

    let genome_idx_set: ahash::HashSet<usize> = pa.genome_col_indices().iter().copied().collect();
    let meta_col_indices: Vec<usize> = (0..pa.header().len())
        .filter(|i| !genome_idx_set.contains(i))
        .collect();

    let n_cols = pa.header().len();
    let mut row: Vec<String> = vec![String::new(); n_cols];

    for ci_val in 0..pa.n_clusters() {
        if clusters_to_skip.contains(&ci_val) {
            continue;
        }

        let ci = ClusterIdx(ci_val);
        for (mi, &col_idx) in meta_col_indices.iter().enumerate() {
            let meta = pa.metadata(ci);
            if mi < meta.len() {
                row[col_idx].clone_from(&meta[mi]);
            } else {
                row[col_idx].clear();
            }
        }

        for (gi_val, &col_idx) in pa.genome_col_indices().iter().enumerate() {
            let gi = GenomeIdx(gi_val);
            row[col_idx].clear();
            if pa.is_present(ci, gi) {
                row[col_idx].push_str("PRESENT_INPUT");
            } else if let Some(rescue) = inputs.accepted_rescue_map.get(&(ci, gi)) {
                row[col_idx].push_str(accepted_rescue_evidence_state(
                    rescue,
                    pa,
                    inputs.gene_to_cluster_maps,
                ));
            } else if let Some(rescue) = synteny_rejected_map.get(&(ci, gi)) {
                let query_cluster = pa.cluster_name(rescue.cluster_idx);
                let annotation_cluster =
                    annotation_cluster_name(rescue, pa, inputs.gene_to_cluster_maps);
                row[col_idx].push_str(evidence_state(
                    rescue,
                    false,
                    query_cluster,
                    &annotation_cluster,
                    "synteny_fail",
                ));
            } else if let Some(rescue) = paralog_rejected_map.get(&(ci, gi)) {
                let query_cluster = pa.cluster_name(rescue.cluster_idx);
                let annotation_cluster =
                    annotation_cluster_name(rescue, pa, inputs.gene_to_cluster_maps);
                row[col_idx].push_str(evidence_state(
                    rescue,
                    false,
                    query_cluster,
                    &annotation_cluster,
                    "paralog_without_synteny",
                ));
            } else if nsc_presence.contains(&(ci_val, gi_val)) {
                row[col_idx].push_str("ALIAS_CLUSTER");
            } else {
                row[col_idx].push_str("ABSENT");
            }
        }

        csv_writer.write_record(row.iter().map(|s| s.as_str()))?;
    }
    csv_writer.flush()?;
    Ok(())
}

#[cfg(test)]
pub(super) fn copy_text_maybe_gzip_with_header(
    path: &std::path::Path,
    writer: &mut impl Write,
) -> Result<()> {
    // Test-only, and it goes through the shared opener too. A test helper that
    // reads gzip differently from the code under test is a helper that can hide
    // the very defect it is meant to exercise.
    let mut reader: Box<dyn BufRead> = crate::parser::open_maybe_gzip(path)?;

    let mut prefix = Vec::with_capacity(256);
    let mut line = Vec::with_capacity(4096);
    let mut has_gff_header = false;
    let mut ended_newline = true;

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        ended_newline = line.ends_with(b"\n");

        let text = String::from_utf8_lossy(&line);
        let trimmed = text.trim_start();
        if trimmed.is_empty() {
            prefix.extend_from_slice(&line);
            continue;
        }
        has_gff_header = trimmed.starts_with("##gff-version");

        prefix.extend_from_slice(&line);
        break;
    }

    if !has_gff_header {
        writeln!(writer, "##gff-version 3")?;
    }
    writer.write_all(&prefix)?;

    let mut buf = Vec::with_capacity(64 * 1024);
    loop {
        buf.clear();
        let n = reader.read_until(b'\n', &mut buf)?;
        if n == 0 {
            break;
        }
        ended_newline = buf.ends_with(b"\n");
        writer.write_all(&buf)?;
    }

    if !ended_newline {
        writeln!(writer)?;
    }
    Ok(())
}

pub(super) fn copy_gff_records_before_fasta_with_header(
    path: &std::path::Path,
    writer: &mut impl Write,
) -> Result<Option<Vec<u8>>> {
    // Through the shared opener rather than a second copy of the same decision.
    // This function held its own extension test and its own `GzDecoder`, so a
    // bgzipped annotation was truncated here as well, and the corrected GFF3
    // this writes kept the records of the first member and dropped the rest.
    let mut reader: Box<dyn BufRead> = crate::parser::open_maybe_gzip(path)?;

    let mut prefix = Vec::with_capacity(256);
    let mut line = Vec::with_capacity(4096);
    let mut has_gff_header = false;
    let mut ended_newline = true;
    let mut fasta_tail: Option<Vec<u8>> = None;

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        ended_newline = line.ends_with(b"\n");

        let text = String::from_utf8_lossy(&line);
        let trimmed = text.trim_start();
        if trimmed.starts_with("##FASTA") {
            let mut tail = Vec::new();
            tail.extend_from_slice(&line);
            reader.read_to_end(&mut tail)?;
            fasta_tail = Some(tail);
            break;
        }
        if trimmed.is_empty() {
            prefix.extend_from_slice(&line);
            continue;
        }
        has_gff_header = trimmed.starts_with("##gff-version");

        prefix.extend_from_slice(&line);
        break;
    }

    if !has_gff_header {
        writeln!(writer, "##gff-version 3")?;
    }
    writer.write_all(&prefix)?;

    let mut buf = Vec::with_capacity(64 * 1024);
    while fasta_tail.is_none() {
        buf.clear();
        let n = reader.read_until(b'\n', &mut buf)?;
        if n == 0 {
            break;
        }
        ended_newline = buf.ends_with(b"\n");
        let text = String::from_utf8_lossy(&buf);
        if text.trim_start().starts_with("##FASTA") {
            let mut tail = Vec::new();
            tail.extend_from_slice(&buf);
            reader.read_to_end(&mut tail)?;
            fasta_tail = Some(tail);
            break;
        }
        writer.write_all(&buf)?;
    }

    if !ended_newline {
        writeln!(writer)?;
    }
    Ok(fasta_tail)
}

pub(super) fn gff_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '%' => out.push_str("%25"),
            ';' => out.push_str("%3B"),
            '=' => out.push_str("%3D"),
            ',' => out.push_str("%2C"),
            '\t' => out.push_str("%09"),
            '\n' | '\r' => out.push(' '),
            ' ' => out.push_str("%20"),
            _ => out.push(ch),
        }
    }
    out
}

pub(super) fn safe_id_component(value: &str) -> String {
    value
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | ':') {
                ch
            } else {
                '_'
            }
        })
        .collect()
}

pub(super) fn write_corrected_gffs(
    rescues: &[RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    gff_paths: &HashMap<GenomeIdx, PathBuf>,
    output_dir: &std::path::Path,
) -> Result<()> {
    let gff_dir = output_dir.join("corrected_gff");
    match fs::remove_dir_all(&gff_dir) {
        Ok(()) => {}
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {}
        Err(err) => {
            return Err(err).with_context(|| format!("Cannot remove {}", gff_dir.display()))
        }
    }
    fs::create_dir_all(&gff_dir)?;

    let mut by_genome: HashMap<GenomeIdx, Vec<&RescueEvent>> = HashMap::new();
    for rescue in rescues {
        by_genome.entry(rescue.genome_idx).or_default().push(rescue);
    }

    let mut genome_indices: Vec<GenomeIdx> = by_genome.keys().copied().collect();
    genome_indices.sort();

    for gi in genome_indices {
        let genome_name = pa.genome_name(gi);
        let path = gff_dir.join(format!("{}.recore.gff3", genome_name));
        let file = fs::File::create(&path)?;
        let mut w = BufWriter::new(file);

        let fasta_tail = if let Some(source_gff) = gff_paths.get(&gi) {
            copy_gff_records_before_fasta_with_header(source_gff, &mut w)?
        } else {
            writeln!(w, "##gff-version 3")?;
            None
        };

        writeln!(w, "###")?;
        writeln!(w, "# Recore rescue features appended below")?;

        let mut rows = by_genome.remove(&gi).unwrap_or_default();
        rows.sort_by(|a, b| {
            a.contig
                .cmp(&b.contig)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
                .then(a.cluster_idx.cmp(&b.cluster_idx))
        });

        for rescue in rows {
            let cluster = pa.cluster_name(rescue.cluster_idx);
            let feature_id = format!(
                "recore_{}_{}_{}_{}",
                safe_id_component(cluster),
                safe_id_component(&rescue.contig),
                rescue.start + 1,
                rescue.end
            );
            let synteny = match rescue.synteny_ok {
                Some(true) => "PASS",
                Some(false) => "FAIL",
                None => "NA",
            };
            let mut attrs = vec![
                format!("ID={}", gff_escape(&feature_id)),
                format!("gene_cluster={}", gff_escape(cluster)),
                format!(
                    "rescue_type={}",
                    gff_escape(&rescue.rescue_type.to_string())
                ),
                format!("identity={:.1}", rescue.identity),
                format!("coverage={:.1}", rescue.coverage),
                format!("confidence={:.3}", rescue.confidence),
                format!("synteny={}", synteny),
            ];
            if let Some(annotation_match) = rescue.annotation_match.as_deref() {
                attrs.push(format!("annotation_match={}", gff_escape(annotation_match)));
            }
            if rescue.fragment_evidence.n_fragments > 0 {
                attrs.push(format!(
                    "fragment_n={}",
                    rescue.fragment_evidence.n_fragments
                ));
                attrs.push(format!(
                    "fragment_query_gap_bp={}",
                    rescue.fragment_evidence.total_query_gap_bp
                ));
                attrs.push(format!(
                    "fragment_target_gap_bp={}",
                    rescue.fragment_evidence.total_target_gap_bp
                ));
                if let Some(note) = rescue.fragment_evidence.note.as_deref() {
                    attrs.push(format!("fragment_note={}", gff_escape(note)));
                }
            }
            attrs.push("Note=Recore_rescue_feature".to_string());

            writeln!(
                w,
                "{}\tRecore\tCDS\t{}\t{}\t.\t{}\t0\t{}",
                rescue.contig,
                rescue.start + 1,
                rescue.end,
                rescue.strand,
                attrs.join(";")
            )?;
        }
        if let Some(tail) = fasta_tail {
            writeln!(w, "###")?;
            w.write_all(&tail)?;
            if !tail.ends_with(b"\n") {
                writeln!(w)?;
            }
        }
        w.flush()?;
    }

    info!("Corrected GFF3 files written to {}", gff_dir.display());
    Ok(())
}

pub(super) fn count_core_after(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
) -> usize {
    let n = pa.n_genomes();
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let clusters_to_skip: ahash::HashSet<usize> = nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();

    (0..pa.n_clusters())
        .filter(|&ci_val| {
            if clusters_to_skip.contains(&ci_val) {
                return false;
            }
            let ci = ClusterIdx(ci_val);
            let present = pa.present_count(ci);
            if present == n {
                return true;
            }
            let extra = (0..n)
                .filter(|&gi_val| {
                    let gi = GenomeIdx(gi_val);
                    !pa.is_present(ci, gi)
                        && (rescue_map.contains_key(&(ci, gi))
                            || nsc_presence.contains(&(ci_val, gi_val)))
                })
                .count();
            present + extra == n
        })
        .count()
}

pub(super) fn count_soft_core_after(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
) -> usize {
    let n = pa.n_genomes();
    let threshold = (n as f64 * 0.95).ceil() as usize;
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let clusters_to_skip: ahash::HashSet<usize> = nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();

    (0..pa.n_clusters())
        .filter(|&ci_val| {
            if clusters_to_skip.contains(&ci_val) {
                return false;
            }
            let ci = ClusterIdx(ci_val);
            let present = pa.present_count(ci);
            if present >= threshold {
                return true;
            }
            let extra = (0..n)
                .filter(|&gi_val| {
                    let gi = GenomeIdx(gi_val);
                    !pa.is_present(ci, gi)
                        && (rescue_map.contains_key(&(ci, gi))
                            || nsc_presence.contains(&(ci_val, gi_val)))
                })
                .count();
            present + extra >= threshold
        })
        .count()
}

pub(super) struct SummaryStats {
    pub(super) core_before: usize,
    pub(super) total_gaps: usize,
    pub(super) filtered_gaps: usize,
    pub(super) skipped_by_prevalence: usize,
    pub(super) searchable_gaps: usize,
    pub(super) missing_fasta_gaps: usize,
    pub(super) dropped: usize,
    pub(super) missing_cds: usize,
    pub(super) split: usize,
    pub(super) fragmented: usize,
    pub(super) alignment_split: usize,
    pub(super) start_variant: usize,
    pub(super) paralog: usize,
    pub(super) pseudogene: usize,
    pub(super) synteny_rejected: usize,
    pub(super) paralog_rejected: usize,
    pub(super) matrix_policy: MatrixPolicy,
    pub(super) matrix_corrected: usize,
    pub(super) nsc_matrix_corrected: usize,
    pub(super) synteny_enabled: bool,
    pub(super) nsc_clusters: usize,
    pub(super) nsc_merged_genes: usize,
    pub(super) skipped_no_reference: usize,
    pub(super) remaining_gaps: usize,
}

pub(super) fn write_nsc_report(
    nsc: &nscluster::NsClusterResult,
    cluster_names: &[String],
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("nscluster_report.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "ns_cluster\trepresentative\tn_members\tmembers")?;
    for (i, nsc_cl) in nsc.merged_clusters.iter().enumerate() {
        let member_names: Vec<&str> = nsc_cl
            .members
            .iter()
            .map(|ci| cluster_names[ci.val()].as_str())
            .collect();
        write_tsv_fields(
            &mut w,
            &[
                format!("NS_{}", i),
                nsc_cl.representative.clone(),
                nsc_cl.members.len().to_string(),
                member_names.join(";"),
            ],
        )?;
    }
    w.flush()?;
    Ok(())
}

pub(super) fn write_summary(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    rescues: &[RescueEvent],
    stats: &SummaryStats,
    nsc_result: &Option<nscluster::NsClusterResult>,
    output_dir: &std::path::Path,
) -> Result<()> {
    let core_after = count_core_after(pa, rescue_map, nsc_result);
    let soft_core_before = pa.soft_core_size();
    let soft_core_after = count_soft_core_after(pa, rescue_map, nsc_result);

    let path = output_dir.join("summary.txt");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "Recore Summary")?;
    writeln!(w, "{}", "=".repeat(50))?;
    writeln!(w, "Clusters:              {}", pa.n_clusters())?;
    writeln!(w, "Genomes:               {}", pa.n_genomes())?;
    writeln!(w, "Total gaps:            {}", stats.total_gaps)?;
    writeln!(w, "Gaps passing filter:   {}", stats.filtered_gaps)?;
    writeln!(w, "Gaps searched:         {}", stats.searchable_gaps)?;
    if stats.missing_fasta_gaps > 0 {
        writeln!(w, "Gaps skipped (FASTA):  {}", stats.missing_fasta_gaps)?;
    }
    if stats.skipped_no_reference > 0 {
        writeln!(w, "Gaps skipped (no seq): {}", stats.skipped_no_reference)?;
    }
    writeln!(w, "Skipped (prevalence):  {}", stats.skipped_by_prevalence)?;
    writeln!(w, "Rescued:               {}", rescues.len())?;
    writeln!(w, "Matrix policy:         {:?}", stats.matrix_policy)?;
    let total_matrix_corrected = stats.matrix_corrected + stats.nsc_matrix_corrected;
    writeln!(w, "Matrix-corrected:      {}", total_matrix_corrected)?;
    writeln!(w, "  rescue cells:        {}", stats.matrix_corrected)?;
    if stats.nsc_matrix_corrected > 0 {
        writeln!(w, "  NSC merged cells:    {}", stats.nsc_matrix_corrected)?;
    }
    if stats.synteny_rejected > 0 {
        writeln!(w, "Rejected (synteny):    {}", stats.synteny_rejected)?;
    }
    if stats.paralog_rejected > 0 {
        writeln!(w, "Rejected (paralog):    {}", stats.paralog_rejected)?;
    }
    // Absent cells left in the corrected matrix as written, so that
    // Total gaps - Matrix-corrected = Remaining gaps when nothing was merged.
    writeln!(w, "Remaining gaps:        {}", stats.remaining_gaps)?;
    writeln!(w)?;
    let core_diff = core_after as i64 - stats.core_before as i64;
    let soft_diff = soft_core_after as i64 - soft_core_before as i64;
    writeln!(
        w,
        "Core genome (strict):  {} → {} ({:+})",
        stats.core_before, core_after, core_diff
    )?;
    writeln!(
        w,
        "Soft-core (≥95%):      {} → {} ({:+})",
        soft_core_before, soft_core_after, soft_diff
    )?;
    writeln!(w)?;
    writeln!(w, "Rescue types:")?;
    writeln!(w, "  DROPPED:       {}", stats.dropped)?;
    writeln!(w, "  MISSING_CDS:   {}", stats.missing_cds)?;
    writeln!(w, "  SPLIT:         {}", stats.split)?;
    writeln!(w, "  FRAGMENTED:    {}", stats.fragmented)?;
    writeln!(w, "  ALIGNMENT_SPLIT: {}", stats.alignment_split)?;
    writeln!(w, "  START_VARIANT: {}", stats.start_variant)?;
    writeln!(w, "  PARALOG:       {}", stats.paralog)?;
    writeln!(w, "  PSEUDOGENE:    {}", stats.pseudogene)?;
    if stats.synteny_enabled {
        let syn_pass = rescues
            .iter()
            .filter(|r| r.synteny_ok == Some(true))
            .count();
        let syn_fail = rescues
            .iter()
            .filter(|r| r.synteny_ok == Some(false))
            .count();
        writeln!(w)?;
        writeln!(w, "Synteny validation:")?;
        writeln!(w, "  PASS:  {}", syn_pass)?;
        writeln!(w, "  FAIL:  {}", syn_fail)?;
        if stats.synteny_rejected > 0 {
            writeln!(w, "  REJECTED_FAIL: {}", stats.synteny_rejected)?;
        }
    }
    if stats.paralog_rejected > 0 {
        writeln!(w)?;
        writeln!(w, "Orthology rejection:")?;
        writeln!(w, "  PARALOG_WITHOUT_SYNTENY: {}", stats.paralog_rejected)?;
    }
    if stats.nsc_clusters > 0 {
        writeln!(w)?;
        writeln!(w, "Nucleotide similarity clustering:")?;
        writeln!(w, "  NS clusters:     {}", stats.nsc_clusters)?;
        writeln!(w, "  Genes merged:    {}", stats.nsc_merged_genes)?;
        writeln!(
            w,
            "  Clusters removed: {}",
            stats.nsc_merged_genes - stats.nsc_clusters
        )?;
    }
    w.flush()?;
    Ok(())
}

pub(super) fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\0' => out.push_str("\\u0000"),
            c if c.is_control() => {
                use std::fmt::Write as _;
                write!(out, "\\u{:04x}", c as u32).expect("write to String is infallible");
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

pub(super) fn json_path(path: &std::path::Path) -> String {
    json_escape(&path.display().to_string())
}

/// Version of the output set, shared by run_manifest.json and the report so the
/// two cannot name different schemas.
pub(super) const OUTPUT_SCHEMA: &str = "recore.outputs.v7";

/// The current time as an ISO 8601 UTC timestamp, without a date crate.
pub(super) fn utc_timestamp() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let (days, rem) = ((secs / 86_400) as i64, secs % 86_400);
    // Civil date from days since 1970-01-01 (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        rem / 3_600,
        rem % 3_600 / 60,
        rem % 60
    )
}

/// The command line this process was started with, quoted where needed so it
/// can be pasted back into a shell.
pub(super) fn command_line() -> String {
    std::env::args()
        .map(|arg| {
            if !arg.is_empty()
                && arg
                    .chars()
                    .all(|c| c.is_ascii_alphanumeric() || "-_./=:,+@%".contains(c))
            {
                arg
            } else {
                format!("'{}'", arg.replace('\'', "'\\''"))
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

pub(super) fn write_run_manifest(
    args: &Args,
    pa: &crate::parser::PresenceAbsence,
    stats: &SummaryStats,
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("run_manifest.json");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    let total_matrix_corrected = stats.matrix_corrected + stats.nsc_matrix_corrected;
    let rescued = stats.dropped
        + stats.missing_cds
        + stats.split
        + stats.fragmented
        + stats.alignment_split
        + stats.start_variant
        + stats.paralog
        + stats.pseudogene;

    writeln!(w, "{{")?;
    writeln!(w, "  \"tool\": \"recore\",")?;
    writeln!(
        w,
        "  \"tool_version\": {},",
        json_escape(env!("CARGO_PKG_VERSION"))
    )?;
    writeln!(w, "  \"schema_version\": {},", json_escape(OUTPUT_SCHEMA))?;
    writeln!(w, "  \"created_utc\": {},", json_escape(&utc_timestamp()))?;
    writeln!(w, "  \"command_line\": {},", json_escape(&command_line()))?;
    writeln!(w, "  \"inputs\": {{")?;
    writeln!(
        w,
        "    \"presence_absence\": {},",
        json_path(&args.presence_absence)
    )?;
    writeln!(
        w,
        "    \"pan_reference\": {},",
        json_path(&args.pan_reference)
    )?;
    writeln!(w, "    \"genomes\": {},", json_path(&args.genomes))?;
    if let Some(path) = args.annotations.as_deref() {
        writeln!(w, "    \"annotations\": {},", json_path(path))?;
    } else {
        writeln!(w, "    \"annotations\": null,")?;
    }
    if let Some(tree) = args.tree.as_deref() {
        writeln!(w, "    \"tree\": {},", json_path(tree))?;
    } else {
        writeln!(w, "    \"tree\": null,")?;
    }
    if let Some(mobile_ref) = args.mobile_element_reference.as_deref() {
        writeln!(
            w,
            "    \"mobile_element_reference\": {}",
            json_path(mobile_ref)
        )?;
    } else {
        writeln!(w, "    \"mobile_element_reference\": null")?;
    }
    writeln!(w, "  }},")?;

    writeln!(w, "  \"parameters\": {{")?;
    writeln!(w, "    \"min_identity\": {:.3},", args.min_identity)?;
    writeln!(w, "    \"min_coverage\": {:.3},", args.min_coverage)?;
    writeln!(w, "    \"min_prevalence\": {:.6},", args.min_prevalence)?;
    writeln!(w, "    \"preset\": {},", json_escape(&args.preset))?;
    writeln!(w, "    \"threads\": {},", args.threads)?;
    writeln!(
        w,
        "    \"matrix_policy\": {},",
        json_escape(&format!("{:?}", args.matrix_policy))
    )?;
    writeln!(w, "    \"multi_allele\": {},", args.multi_allele)?;
    writeln!(w, "    \"synteny\": {},", args.synteny)?;
    writeln!(
        w,
        "    \"write_corrected_gff\": {},",
        args.write_corrected_gff
    )?;
    writeln!(w, "    \"write_html_report\": {},", args.write_html_report)?;
    if let Some(reference_genome) = args.reference_genome.as_deref() {
        writeln!(
            w,
            "    \"reference_genome\": {},",
            json_escape(reference_genome)
        )?;
    } else {
        writeln!(w, "    \"reference_genome\": null,")?;
    }
    writeln!(w, "    \"nscluster\": {},", args.nscluster)?;
    writeln!(w, "    \"nsc_engine\": {},", json_escape(&args.nsc_engine))?;
    writeln!(w, "    \"nsc_kmer_size\": {},", args.nsc_kmer_size)?;
    writeln!(w, "    \"nsc_threshold\": {:.6},", args.nsc_threshold)?;
    writeln!(w, "    \"nsc_min_seq_id\": {:.6},", args.nsc_min_seq_id)?;
    writeln!(w, "    \"nsc_coverage\": {:.6},", args.nsc_coverage)?;
    writeln!(
        w,
        "    \"mobile_element_min_identity\": {:.3},",
        args.mobile_element_min_identity
    )?;
    writeln!(
        w,
        "    \"mobile_element_min_coverage\": {:.3},",
        args.mobile_element_min_coverage
    )?;
    writeln!(
        w,
        "    \"mobile_element_window\": {}",
        args.mobile_element_window
    )?;
    writeln!(w, "  }},")?;

    writeln!(w, "  \"counts\": {{")?;
    writeln!(w, "    \"clusters\": {},", pa.n_clusters())?;
    writeln!(w, "    \"genomes\": {},", pa.n_genomes())?;
    writeln!(w, "    \"total_gaps\": {},", stats.total_gaps)?;
    writeln!(w, "    \"gaps_passing_filter\": {},", stats.filtered_gaps)?;
    writeln!(w, "    \"gaps_searched\": {},", stats.searchable_gaps)?;
    writeln!(
        w,
        "    \"skipped_by_prevalence\": {},",
        stats.skipped_by_prevalence
    )?;
    writeln!(
        w,
        "    \"missing_fasta_gaps\": {},",
        stats.missing_fasta_gaps
    )?;
    writeln!(
        w,
        "    \"missing_reference_gaps\": {},",
        stats.skipped_no_reference
    )?;
    writeln!(w, "    \"remaining_gaps\": {},", stats.remaining_gaps)?;
    writeln!(w, "    \"rescued\": {},", rescued)?;
    writeln!(
        w,
        "    \"matrix_corrected_total\": {},",
        total_matrix_corrected
    )?;
    writeln!(
        w,
        "    \"matrix_corrected_rescue_cells\": {},",
        stats.matrix_corrected
    )?;
    writeln!(
        w,
        "    \"matrix_corrected_nsc_cells\": {},",
        stats.nsc_matrix_corrected
    )?;
    writeln!(w, "    \"synteny_rejected\": {},", stats.synteny_rejected)?;
    writeln!(w, "    \"paralog_rejected\": {},", stats.paralog_rejected)?;
    writeln!(w, "    \"nsc_clusters\": {},", stats.nsc_clusters)?;
    writeln!(w, "    \"nsc_merged_genes\": {}", stats.nsc_merged_genes)?;
    writeln!(w, "  }}")?;
    writeln!(w, "}}")?;
    w.flush()?;
    Ok(())
}

fn html_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(ch),
        }
    }
    out
}

fn format_count(value: usize) -> String {
    let text = value.to_string();
    let mut out = String::with_capacity(text.len() + text.len() / 3);
    for (idx, ch) in text.chars().rev().enumerate() {
        if idx > 0 && idx % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out.chars().rev().collect()
}

fn format_percent(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        return "NA".to_string();
    }
    format!("{:.1}%", numerator as f64 / denominator as f64 * 100.0)
}

fn count_tsv_column(path: &std::path::Path, column: &str) -> Result<BTreeMap<String, usize>> {
    let mut counts = BTreeMap::new();
    if !path.exists() {
        return Ok(counts);
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let Some(col_idx) = headers.iter().position(|header| header == column) else {
        return Ok(counts);
    };
    for record in reader.records() {
        let record = record?;
        let value = record.get(col_idx).unwrap_or("").trim();
        let key = if value.is_empty() { "NA" } else { value };
        *counts.entry(key.to_string()).or_default() += 1;
    }
    Ok(counts)
}

/// genome_qc.tsv as (status, flags) per matrix genome, for the genome card.
fn genome_qc_rows(
    path: &std::path::Path,
    pa: &crate::parser::PresenceAbsence,
) -> Result<Vec<(String, String)>> {
    let mut rows = vec![(String::new(), String::new()); pa.n_genomes()];
    if !path.exists() {
        return Ok(rows);
    }
    let index: HashMap<&str, usize> = pa
        .genome_names()
        .iter()
        .enumerate()
        .map(|(gi, name)| (name.as_str(), gi))
        .collect();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| headers.iter().position(|h| h == name);
    let (Some(genome_col), Some(status_col)) = (col("genome"), col("qc_status")) else {
        return Ok(rows);
    };
    let flags_col = col("qc_flags");
    for record in reader.records() {
        let record = record?;
        if let Some(&gi) = record.get(genome_col).and_then(|g| index.get(g)) {
            rows[gi] = (
                record.get(status_col).unwrap_or("").to_string(),
                flags_col
                    .and_then(|c| record.get(c))
                    .unwrap_or("")
                    .to_string(),
            );
        }
    }
    Ok(rows)
}

/// What the report says a reader can put in a methods section, built from the
/// run's own settings and counts so it cannot drift from what was done.
#[allow(clippy::too_many_arguments)]
fn methods_text(
    args: &Args,
    pa: &crate::parser::PresenceAbsence,
    stats: &SummaryStats,
    corrected: usize,
    kept_out: usize,
    rejected: usize,
    merged_by_nsc: usize,
) -> String {
    let mut text = format!(
        "False gene absences in the gene presence/absence matrix ({} clusters, {} genomes) were corrected with Recore v{}. ",
        format_count(pa.n_clusters()),
        format_count(pa.n_genomes()),
        env!("CARGO_PKG_VERSION")
    );
    text.push_str(&format!(
        "Each absent cell of a cluster present in at least {}% of genomes was searched by aligning the cluster's representative sequence{} to that genome's assembly with minimap2 (preset {}); a hit counted when it reached at least {}% identity over at least {}% of the query. ",
        trim_number(args.min_prevalence * 100.0),
        if args.multi_allele { " and every allele" } else { "" },
        args.preset,
        trim_number(args.min_identity),
        trim_number(args.min_coverage)
    ));
    if args.annotations.is_some() {
        text.push_str("Hits were classified against the genome annotation as missing annotations, split or fragmented genes, start variants or pseudogenes");
        text.push_str(if args.synteny {
            ", and kept only when the flanking gene context (synteny) agreed with the cluster's. "
        } else {
            ". "
        });
    }
    text.push_str(&format!(
        "Under the {} matrix policy, {} cells were set to present; {} further cells with sequence evidence were kept out of the binary matrix, and the synteny and paralog guards rejected {} candidates. ",
        format!("{:?}", stats.matrix_policy).to_lowercase(),
        format_count(corrected),
        format_count(kept_out),
        format_count(rejected)
    ));
    if merged_by_nsc > 0 {
        text.push_str(&format!(
            "Clusters identical at the nucleotide level were merged by similarity clustering ({} clusters merged). ",
            format_count(merged_by_nsc)
        ));
    }
    text.push_str("The evidence behind every corrected and kept-out cell is listed in orthology_evidence.tsv.");
    text
}

/// 90.0 as "90", 0.1 as "0.1": thresholds read as people write them.
pub(super) fn trim_number(value: f64) -> String {
    let text = format!("{value:.3}");
    text.trim_end_matches('0').trim_end_matches('.').to_string()
}

#[derive(Debug, Clone)]
struct CandidateReportRow {
    cluster: String,
    priority: String,
    score: f64,
    mechanism: String,
    support: String,
    phylo: String,
    wording: String,
}

#[derive(Debug, Clone, Default)]
struct ExplorerEvolutionRow {
    mechanism: String,
    detail: String,
    support: String,
    phylo: String,
    wording: String,
    reference_event_class: String,
    target_reference_cluster: String,
    target_reference_gene: String,
}

#[derive(Debug, Clone)]
struct ExplorerSyntenyBlock {
    rank: usize,
    n_genomes: usize,
    fraction: String,
    upstream: Vec<String>,
    downstream: Vec<String>,
    interpretation: String,
    examples: String,
    /// Matrix genome indices, sorted. Kept as indices from the moment the TSV is
    /// read: one String per carrier name was most of the memory a full-synteny
    /// run used, and merging blocks compared names pairwise.
    carrier_genomes: Vec<u32>,
    representative_genome: String,
    representative_contig: String,
    representative_start: String,
    representative_end: String,
    representative_length: String,
    genes: Vec<ExplorerSyntenyGene>,
    reference_genes: Vec<ExplorerSyntenyGene>,
}

#[derive(Debug, Clone)]
struct ExplorerSyntenyGene {
    cluster: String,
    gene_id: String,
    strand: String,
    contig: String,
    start: String,
    end: String,
    length: String,
}

fn explorer_evolution_rows(
    path: &std::path::Path,
) -> Result<HashMap<String, ExplorerEvolutionRow>> {
    let mut rows = HashMap::new();
    if !path.exists() {
        return Ok(rows);
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let idx = |name: &str| headers.iter().position(|header| header == name);
    let Some(cluster_idx) = idx("gene_cluster") else {
        return Ok(rows);
    };
    let mechanism_idx = idx("molecular_mechanism");
    let detail_idx = idx("mechanism_detail");
    let support_idx = idx("support_level");
    let phylo_idx = idx("phylogenetic_pattern");
    let wording_idx = idx("recommended_wording");
    let reference_event_class_idx = idx("reference_event_class");
    let target_reference_cluster_idx = idx("target_reference_cluster");
    let target_reference_gene_idx = idx("target_reference_gene");

    for record in reader.records() {
        let record = record?;
        let cluster = record.get(cluster_idx).unwrap_or("").trim();
        if cluster.is_empty() {
            continue;
        }
        rows.insert(
            cluster.to_string(),
            ExplorerEvolutionRow {
                mechanism: mechanism_idx
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string(),
                detail: detail_idx
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string(),
                support: support_idx
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string(),
                phylo: phylo_idx
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string(),
                wording: wording_idx
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string(),
                reference_event_class: reference_event_class_idx
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string(),
                target_reference_cluster: target_reference_cluster_idx
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string(),
                target_reference_gene: target_reference_gene_idx
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string(),
            },
        );
    }
    Ok(rows)
}

fn parse_synteny_side(value: &str) -> Vec<String> {
    value
        .split('|')
        .map(str::trim)
        .filter(|item| !item.is_empty() && *item != "NA")
        .map(ToString::to_string)
        .collect()
}

fn parse_synteny_gene_records(value: &str) -> Vec<ExplorerSyntenyGene> {
    value
        .split(';')
        .filter_map(|item| {
            let fields: Vec<_> = item.split(',').collect();
            if fields.len() != 7 {
                return None;
            }
            Some(ExplorerSyntenyGene {
                cluster: fields[0].trim().to_string(),
                gene_id: fields[1].trim().to_string(),
                strand: fields[2].trim().to_string(),
                contig: fields[3].trim().to_string(),
                start: fields[4].trim().to_string(),
                end: fields[5].trim().to_string(),
                length: fields[6].trim().to_string(),
            })
        })
        .collect()
}

fn parse_synteny_carriers(value: &str, genome_index: &HashMap<&str, u32>) -> Vec<u32> {
    let mut carriers: Vec<u32> = value
        .split(';')
        .map(str::trim)
        .filter_map(|name| genome_index.get(name).copied())
        .collect();
    carriers.sort_unstable();
    carriers.dedup();
    carriers
}

fn explorer_synteny_signature(
    genes: &[ExplorerSyntenyGene],
    upstream: &[String],
    downstream: &[String],
    cluster: &str,
) -> String {
    let signature: Vec<&str> = genes
        .iter()
        .filter_map(|gene| {
            let gene_id = gene.gene_id.trim();
            let identity = if !gene_id.is_empty() && gene_id != "H37Rv_context_slot" {
                gene_id
            } else {
                gene.cluster.trim()
            };
            match identity {
                "" | "NA" => None,
                value => Some(value),
            }
        })
        .collect();
    if !signature.is_empty() {
        return signature.join("|");
    }
    upstream
        .iter()
        .map(String::as_str)
        .chain(std::iter::once(cluster))
        .chain(downstream.iter().map(String::as_str))
        .filter(|item| !item.is_empty() && *item != "NA")
        .collect::<Vec<_>>()
        .join("|")
}

fn merge_explorer_synteny_block(
    existing: &mut ExplorerSyntenyBlock,
    incoming: ExplorerSyntenyBlock,
) {
    let incoming_n_genomes = incoming.n_genomes;
    existing.carrier_genomes.extend(incoming.carrier_genomes);
    existing.carrier_genomes.sort_unstable();
    existing.carrier_genomes.dedup();
    existing.n_genomes = if existing.carrier_genomes.is_empty() {
        existing.n_genomes + incoming_n_genomes
    } else {
        existing.carrier_genomes.len()
    };
    let merged_fraction = existing.fraction.parse::<f64>().ok().unwrap_or(0.0)
        + incoming.fraction.parse::<f64>().ok().unwrap_or(0.0);
    if merged_fraction > 0.0 {
        existing.fraction = format!("{:.3}", merged_fraction.min(1.0));
    }
    // The page shows one example locus; keep the first.
    if existing.examples.is_empty() {
        existing.examples = incoming.examples;
    }
}

fn explorer_synteny_blocks(
    path: &std::path::Path,
    limit_per_cluster: usize,
    genome_index: &HashMap<&str, u32>,
) -> Result<HashMap<String, Vec<ExplorerSyntenyBlock>>> {
    let mut rows: HashMap<String, Vec<ExplorerSyntenyBlock>> = HashMap::new();
    let mut signature_indices: HashMap<String, HashMap<String, usize>> = HashMap::new();
    if !path.exists() || limit_per_cluster == 0 {
        return Ok(rows);
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let idx = |name: &str| headers.iter().position(|header| header == name);
    let Some(cluster_idx) = idx("gene_cluster") else {
        return Ok(rows);
    };
    let rank_idx = idx("block_rank");
    let n_genomes_idx = idx("n_genomes");
    let fraction_idx = idx("fraction_of_present");
    let upstream_idx = idx("upstream_block_5");
    let downstream_idx = idx("downstream_block_5");
    let interpretation_idx = idx("block_interpretation");
    let examples_idx = idx("example_genomes");
    let carrier_genomes_idx = idx("carrier_genomes");
    let representative_genome_idx = idx("representative_genome");
    let representative_contig_idx = idx("representative_contig");
    let representative_start_idx = idx("representative_start");
    let representative_end_idx = idx("representative_end");
    let representative_length_idx = idx("representative_length_bp");
    let representative_genes_idx = idx("representative_block_genes");
    let reference_genes_idx = idx("reference_block_genes");

    for record in reader.records() {
        let record = record?;
        let cluster = record.get(cluster_idx).unwrap_or("").trim();
        if cluster.is_empty() {
            continue;
        }
        let upstream = upstream_idx
            .and_then(|i| record.get(i))
            .map(parse_synteny_side)
            .unwrap_or_default();
        let downstream = downstream_idx
            .and_then(|i| record.get(i))
            .map(parse_synteny_side)
            .unwrap_or_default();
        let genes = representative_genes_idx
            .and_then(|i| record.get(i))
            .map(parse_synteny_gene_records)
            .unwrap_or_default();
        let reference_genes = reference_genes_idx
            .and_then(|i| record.get(i))
            .map(parse_synteny_gene_records)
            .unwrap_or_default();
        let signature = explorer_synteny_signature(&genes, &upstream, &downstream, cluster);
        let block = ExplorerSyntenyBlock {
            rank: rank_idx
                .and_then(|i| record.get(i))
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(1),
            n_genomes: n_genomes_idx
                .and_then(|i| record.get(i))
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(0),
            fraction: fraction_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            upstream,
            downstream,
            interpretation: interpretation_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            examples: examples_idx
                .and_then(|i| record.get(i))
                .and_then(|value| value.split(';').map(str::trim).find(|v| !v.is_empty()))
                .unwrap_or("")
                .to_string(),
            carrier_genomes: carrier_genomes_idx
                .and_then(|i| record.get(i))
                .map(|value| parse_synteny_carriers(value, genome_index))
                .unwrap_or_default(),
            representative_genome: representative_genome_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            representative_contig: representative_contig_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            representative_start: representative_start_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            representative_end: representative_end_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            representative_length: representative_length_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            genes,
            reference_genes,
        };
        let entry = rows.entry(cluster.to_string()).or_default();
        let cluster_signatures = signature_indices.entry(cluster.to_string()).or_default();
        if let Some(existing_idx) = cluster_signatures.get(&signature).copied() {
            if let Some(existing) = entry.get_mut(existing_idx) {
                merge_explorer_synteny_block(existing, block);
            }
            continue;
        }
        if entry.len() >= limit_per_cluster {
            continue;
        }
        cluster_signatures.insert(signature, entry.len());
        entry.push(block);
    }
    for blocks in rows.values_mut() {
        blocks.sort_by_key(|block| block.rank);
        blocks.truncate(limit_per_cluster);
    }
    Ok(rows)
}

fn top_candidate_rows(path: &std::path::Path, limit: usize) -> Result<Vec<CandidateReportRow>> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let idx = |name: &str| headers.iter().position(|header| header == name);
    let Some(cluster_idx) = idx("gene_cluster") else {
        return Ok(Vec::new());
    };
    let priority_idx = idx("highlight_priority");
    let score_idx = idx("evolutionary_signal_score");
    let mechanism_idx = idx("molecular_mechanism");
    let support_idx = idx("support_level");
    let phylo_idx = idx("phylogenetic_pattern");
    let wording_idx = idx("recommended_wording");

    let mut rows = Vec::new();
    for record in reader.records() {
        let record = record?;
        let score = score_idx
            .and_then(|i| record.get(i))
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(0.0);
        rows.push(CandidateReportRow {
            cluster: record.get(cluster_idx).unwrap_or("").to_string(),
            priority: priority_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            score,
            mechanism: mechanism_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            support: support_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            phylo: phylo_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            wording: wording_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
        });
    }
    rows.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| a.cluster.cmp(&b.cluster))
    });
    rows.truncate(limit);
    Ok(rows)
}

fn count_items_html(
    title: &str,
    counts: &BTreeMap<String, usize>,
    empty: &str,
    bar_kind: &str,
) -> String {
    use std::fmt::Write as _;
    let mut html = String::new();
    write!(
        html,
        "<section class='panel'><div class='panel-head'><h2>{}</h2></div>",
        html_escape(title)
    )
    .expect("write to String is infallible");
    if counts.is_empty() {
        write!(
            html,
            "<p class='muted'>{}</p></section>",
            html_escape(empty)
        )
        .expect("write to String is infallible");
        return html;
    }
    let total: usize = counts.values().copied().sum::<usize>().max(1);
    html.push_str("<div class='bars'>");
    let mut rows: Vec<_> = counts.iter().collect();
    rows.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    for (label, value) in rows {
        let share = *value as f64 / total as f64 * 100.0;
        // Genome QC is coloured by what each status means, not one hue for all.
        let kind = match (bar_kind, label.as_str()) {
            ("qc", "pass") => "pos",
            ("qc", "review") => "warn",
            ("qc", "fail") => "danger",
            ("qc", _) => "",
            (kind, _) => kind,
        };
        write!(
            html,
            "<div class='bar-row'><div class='bar-label'><span data-tip='{}'>{}</span><strong class='num'>{} <small>{:.1}%</small></strong></div><div class='track'><div class='bar {}' style='width:{:.1}%'></div></div></div>",
            html_escape(label),
            html_escape(&plain_label(label)),
            format_count(*value),
            share,
            kind,
            share.clamp(1.5, 100.0)
        )
        .expect("write to String is infallible");
    }
    html.push_str("</div></section>");
    html
}

fn output_file_description(name: &str) -> &'static str {
    match name {
        "gene_presence_absence_recore.csv" => "Corrected binary matrix",
        "gene_presence_absence_recore_states.csv" => {
            "Multistate matrix: input, corrected, NSC and kept-out codes"
        }
        "orthology_evidence.tsv" => {
            "Every rescue and rejected candidate, with its decision and why"
        }
        "rescue_report.tsv" => "One row per rescue: type, identity, coverage and locus",
        "evolutionary_events.tsv" => "Cluster mechanisms, priorities and phylogenetic patterns",
        "genome_qc.tsv" => "Per-genome inputs, correction burden and QC status",
        "reference_genome_audit.tsv" => "Reference-absence audit (--reference-genome)",
        "cluster_synteny_blocks.tsv" => "Unique flanking-gene blocks per cluster",
        "cluster_locus_context.tsv" => "Locus context per cluster and genome",
        "cluster_locus_summary.tsv" => "Locus context summary per cluster",
        "cluster_relationships.tsv" => "Split, alias, paralog and NSC links between clusters",
        "cluster_components.tsv" => "Components of related clusters",
        "cluster_component_members.tsv" => "Members of each cluster component",
        "run_manifest.json" => "Inputs, parameters, counts and the command line",
        "summary.txt" => "Plain-text run summary",
        _ => "",
    }
}

/// One link per file actually in the output directory, with a description when
/// Recore knows the file and its row count (or size, for large files).
pub(super) fn output_file_entries(output_dir: &std::path::Path) -> Vec<String> {
    let Ok(read) = std::fs::read_dir(output_dir) else {
        return Vec::new();
    };
    let mut names: Vec<(String, bool, u64)> = read
        .filter_map(|entry| entry.ok())
        .filter_map(|entry| {
            let name = entry.file_name().to_str()?.to_string();
            let meta = entry.metadata().ok()?;
            (name != "recore_report.html" && !name.starts_with('.'))
                .then(|| (name, meta.is_dir(), meta.len()))
        })
        .collect();
    // Known outputs first, in the order a reader needs them, then the rest.
    let rank = |name: &str| {
        [
            "gene_presence_absence_recore.csv",
            "orthology_evidence.tsv",
            "gene_presence_absence_recore_states.csv",
            "rescue_report.tsv",
            "genome_qc.tsv",
            "evolutionary_events.tsv",
            "summary.txt",
            "run_manifest.json",
        ]
        .iter()
        .position(|known| *known == name)
        .unwrap_or(usize::MAX)
    };
    names.sort_by(|a, b| rank(&a.0).cmp(&rank(&b.0)).then_with(|| a.0.cmp(&b.0)));
    names
        .into_iter()
        .map(|(name, is_dir, len)| {
            let detail = if is_dir {
                "folder".to_string()
            } else if (name.ends_with(".tsv") || name.ends_with(".csv")) && len < 32 * 1024 * 1024 {
                let rows = std::fs::read(output_dir.join(&name))
                    .map(|bytes| bytecount_lines(&bytes).saturating_sub(1))
                    .unwrap_or(0);
                format!("{} {}", format_count(rows), if rows == 1 { "row" } else { "rows" })
            } else {
                format_bytes(len)
            };
            let description = output_file_description(&name);
            format!(
                "<a class='file' href='{0}' target='_blank' rel='noopener'><span class='file-name'>{1}</span><small>{2}{3}{4}</small></a>",
                html_escape(&name),
                html_escape(&name).replace('_', "_<wbr>"),
                html_escape(description),
                if description.is_empty() { "" } else { " &middot; " },
                html_escape(&detail)
            )
        })
        .collect()
}

fn bytecount_lines(bytes: &[u8]) -> usize {
    let newlines = bytes.iter().filter(|&&b| b == b'\n').count();
    newlines + usize::from(!bytes.is_empty() && !bytes.ends_with(b"\n"))
}

fn format_bytes(len: u64) -> String {
    const UNITS: [&str; 4] = ["B", "KB", "MB", "GB"];
    let mut value = len as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit < UNITS.len() - 1 {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{len} B")
    } else {
        format!("{value:.1} {}", UNITS[unit])
    }
}

const RECORE_REPO_URL: &str = "https://github.com/PathoGenOmics-Lab/recore";
// The documentation button in the report. The MkDocs site is built and deployed
// on every push to main, but Pages will not serve it while the repository is
// private on a free plan, so this points at the pages in the repository, which
// render on GitHub and resolve today. Once the site is live, this is the one
// line to change.
const RECORE_DOCS_URL: &str = "https://github.com/PathoGenOmics-Lab/recore/tree/main/docs";
const REPORT_LOGO: &str = include_str!("../../assets/recore_logo.svg");
const REPORT_CSS: &str = include_str!("report_assets/report.css");
const REPORT_UI_JS: &str = include_str!("report_assets/report-ui.js");
const EXPLORER_JS: &str = include_str!("report_assets/explorer.js");
/// Reads the embedded data before the viewer starts. A failure to parse it, or
/// to start the viewer, used to leave a blank workspace with no message; now it
/// says what happened and where the results are. The parsed text is removed
/// from the page, which frees a copy of the whole payload.
const EXPLORER_BOOT_JS: &str = r#"(function () {
  window.recoreExplorerFailed = function (error) {
    var main = document.getElementById("explorer");
    if (!main) return;
    var box = document.createElement("section");
    box.className = "panel full explorer explorer-fallback explorer-failed";
    box.setAttribute("role", "alert");
    var title = document.createElement("h2");
    title.textContent = "The interactive workspace could not start";
    var text = document.createElement("p");
    text.textContent = "Reason: " + (error && error.message ? error.message : String(error)) +
      ". Every number is still in Run details (press D) and in the TSV files of the output directory.";
    box.appendChild(title);
    box.appendChild(text);
    main.replaceWith(box);
    if (window.console) console.error(error);
  };
  var element = document.getElementById("explorerData");
  try {
    window.recoreExplorerData = JSON.parse(element.textContent);
    element.remove();
  } catch (error) {
    window.recoreExplorerFailed(error);
  }
})();"#;
// The viewer draws only its visible window and ships the matrix in compact
// per-cluster sets, so the ceiling is now the payload the browser has to parse
// rather than anything about the rendering.
const EXPLORER_MAX_TIPS: usize = 20_000;
const EXPLORER_MAX_CELLS: usize = 1_200_000_000;
const EXPLORER_MAX_SYNTENY_BLOCKS_PER_CLUSTER: usize = 10;
/// Above this the browser may fail to hold the embedded JSON at all, which left
/// a blank page; synteny detail, the largest part, is dropped first.
const EXPLORER_MAX_JSON_BYTES: usize = 256 * 1024 * 1024;

/// Runs in <head> so `data-theme` is set before the browser paints anything.
const THEME_BOOTSTRAP_JS: &str = "<script>(()=>{try{\
const stored=localStorage.getItem('recore.theme');\
const dark=window.matchMedia&&window.matchMedia('(prefers-color-scheme: dark)').matches;\
document.documentElement.setAttribute('data-theme',stored||(dark?'dark':'light'));\
}catch(e){document.documentElement.setAttribute('data-theme','light');}})();</script>";

/// The mark ships web-ready: `tools/make_logo_variants.py` crops it to the
/// artwork, strips the export's fixed size and `preserveAspectRatio="none"`, and
/// writes the dark-background twin whose inks the stylesheet applies. So there is
/// nothing to do here but hand it to the page.
fn report_logo_svg() -> String {
    REPORT_LOGO.trim().to_string()
}

fn mechanism_color(name: &str) -> String {
    format!("var(--mech-{})", mechanism_family(name))
}

fn priority_badge_class(priority: &str) -> &'static str {
    let p = priority.to_ascii_lowercase();
    if p.contains("filter") || p.contains("low") {
        "filter"
    } else if p.contains("review") || p.contains("candidate") || p.contains("medium") {
        "review"
    } else if p.contains("high") || p.contains("priority") || p.contains("pass") {
        "pass"
    } else {
        "muted"
    }
}

fn pretty_label(value: &str) -> String {
    value.replace('_', " ")
}

/// Pipeline vocabulary in words a reader knows; anything not listed reads with
/// spaces for underscores. The raw value stays in the element's tip.
fn plain_label(value: &str) -> String {
    match value {
        "correct_presence" => "corrected".to_string(),
        "audit_only_policy" => "kept out by the matrix policy".to_string(),
        "audit_only_quality" => "kept out by a quality floor".to_string(),
        "audit_only" => "rejected".to_string(),
        "high_confidence" => "high".to_string(),
        "moderate_confidence" => "moderate".to_string(),
        "low_confidence" => "low".to_string(),
        "filter_or_review_blinking" => "filter or review (unstable)".to_string(),
        "NA" => "not given".to_string(),
        other => pretty_label(other),
    }
}

/// Mechanisms are drawn in six families; the report's markup takes the
/// family's colour token so it follows the theme.
fn mechanism_family(name: &str) -> &'static str {
    match name {
        "ancestral_retained_no_disruption" => "retained",
        "clade_associated_gain" | "clean_insertion_or_acquisition_candidate" => "gain",
        "deletion_loss_candidate" | "reference_or_lineage_loss_candidate" => "loss",
        "IS_interruption_supported"
        | "pseudogenization"
        | "gene_fragmentation_or_structural_break"
        | "reference_gene_disruption_or_truncation" => "disruption",
        "duplication_or_paralog_expansion"
        | "translocation_or_context_shift_candidate"
        | "partial_reference_homology"
        | "reference_gene_allelic_derivative" => "paralog",
        _ => "unresolved",
    }
}

/// JSON made safe to sit inside a `<script>` element.
///
/// In JSON these characters can only occur inside strings, so writing them as
/// `\uXXXX` leaves every value exactly as it was. Escaping only `</` was not
/// enough: `<!--` followed by `<script` anywhere in a genome name or an
/// annotation switches the HTML parser into the script's escaped state, the
/// closing tag is never seen, and the whole report stayed blank without an
/// error.
pub(super) fn script_safe_json(json: &str) -> String {
    let mut out = String::with_capacity(json.len() + json.len() / 64);
    for ch in json.chars() {
        match ch {
            '<' => out.push_str("\\u003c"),
            '>' => out.push_str("\\u003e"),
            '&' => out.push_str("\\u0026"),
            '\u{2028}' => out.push_str("\\u2028"),
            '\u{2029}' => out.push_str("\\u2029"),
            c => out.push(c),
        }
    }
    out
}

fn json_script_escape(value: &str) -> String {
    json_escape(value).replace("</", "<\\/")
}

fn json_string_array(values: &[String]) -> String {
    let mut out = String::from("[");
    for (idx, value) in values.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str(&json_script_escape(value));
    }
    out.push(']');
    out
}

fn json_synteny_gene_array(values: &[ExplorerSyntenyGene]) -> String {
    use std::fmt::Write as _;
    let mut out = String::from("[");
    for (idx, gene) in values.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        write!(
            out,
            "[{},{},{},{},{},{},{}]",
            json_script_escape(&gene.cluster),
            json_script_escape(&gene.gene_id),
            json_script_escape(&gene.strand),
            json_script_escape(&gene.contig),
            json_script_escape(&gene.start),
            json_script_escape(&gene.end),
            json_script_escape(&gene.length)
        )
        .expect("write to String is infallible");
    }
    out.push(']');
    out
}

fn matched_phylo_tips(root: &PhyloNode) -> usize {
    let mut tips = 0;
    let mut stack = vec![root];
    while let Some(node) = stack.pop() {
        if node.children.is_empty() {
            tips += usize::from(node.genome_idx.is_some());
        }
        stack.extend(node.children.iter());
    }
    tips
}

/// Without --tree the rows still need an order: every genome hangs off one
/// root, in matrix order, so the workspace opens and the tree pane lists them.
fn star_tree(pa: &crate::parser::PresenceAbsence) -> PhyloTree {
    PhyloTree {
        root: PhyloNode {
            label: None,
            genome_idx: None,
            branch_length: None,
            children: (0..pa.n_genomes())
                .map(|gi| PhyloNode {
                    label: Some(pa.genome_name(GenomeIdx(gi)).to_string()),
                    genome_idx: Some(GenomeIdx(gi)),
                    branch_length: None,
                    children: Vec::new(),
                })
                .collect(),
        },
    }
}

/// The tree as the explorer reads it: only tips that are matrix genomes, and no
/// one-child nodes (a chain of them collapses into its child, branch lengths
/// summed), so a reroot in the browser cannot grow a nameless tip.
///
/// Iterative throughout: a ladder-shaped tree is as deep as it has tips, and the
/// recursive version overflowed the stack and copied every subtree's text once
/// per level above it.
pub(super) fn pruned_phylo_json(root: &PhyloNode) -> Option<String> {
    use std::fmt::Write as _;

    // Flatten in pre-order; children keep their order.
    let mut nodes: Vec<&PhyloNode> = Vec::new();
    let mut children: Vec<Vec<usize>> = Vec::new();
    let mut stack = vec![(root, usize::MAX)];
    while let Some((node, parent)) = stack.pop() {
        let idx = nodes.len();
        nodes.push(node);
        children.push(Vec::new());
        if parent != usize::MAX {
            children[parent].push(idx);
        }
        for child in node.children.iter().rev() {
            stack.push((child, idx));
        }
    }
    // Children come after their parent, so one backward pass sees every
    // subtree before the node that holds it.
    let mut kept: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    let mut keep = vec![false; nodes.len()];
    for idx in (0..nodes.len()).rev() {
        if children[idx].is_empty() {
            keep[idx] = nodes[idx].genome_idx.is_some();
        } else {
            kept[idx] = children[idx].iter().copied().filter(|&c| keep[c]).collect();
            keep[idx] = !kept[idx].is_empty();
        }
    }
    if !keep[0] {
        return None;
    }
    let collapse = |mut idx: usize| -> (usize, Option<f64>) {
        let mut length = nodes[idx].branch_length;
        while kept[idx].len() == 1 {
            let child = kept[idx][0];
            length = match (length, nodes[child].branch_length) {
                (None, None) => None,
                (a, b) => Some(a.unwrap_or(0.0) + b.unwrap_or(0.0)),
            };
            idx = child;
        }
        (idx, length)
    };

    enum Step {
        Enter(usize, Option<f64>, bool),
        Exit,
    }
    let mut out = String::new();
    let (top, top_length) = collapse(0);
    let mut steps = vec![Step::Enter(top, top_length, false)];
    while let Some(step) = steps.pop() {
        let (idx, length, comma) = match step {
            Step::Exit => {
                out.push_str("]}");
                continue;
            }
            Step::Enter(idx, length, comma) => (idx, length, comma),
        };
        if comma {
            out.push(',');
        }
        let node = nodes[idx];
        let label = json_script_escape(node.label.as_deref().unwrap_or(""));
        let length = length
            .map(|value| value.to_string())
            .unwrap_or_else(|| "null".to_string());
        if kept[idx].is_empty() {
            let gi = node.genome_idx.expect("only tips with a genome are kept");
            write!(
                out,
                "{{\"label\":{},\"genome\":{},\"length\":{},\"children\":[]}}",
                label,
                gi.val(),
                length
            )
            .expect("write to String is infallible");
            continue;
        }
        write!(
            out,
            "{{\"label\":{},\"genome\":null,\"length\":{},\"children\":[",
            label, length
        )
        .expect("write to String is infallible");
        steps.push(Step::Exit);
        for (k, &child) in kept[idx].iter().enumerate().rev() {
            let (child, child_length) = collapse(child);
            steps.push(Step::Enter(child, child_length, k > 0));
        }
    }
    Some(out)
}

fn metadata_column_index(pa: &crate::parser::PresenceAbsence, name: &str) -> Option<usize> {
    let mut meta_idx = 0usize;
    for (col_idx, header) in pa.header().iter().enumerate() {
        if pa.genome_col_indices().contains(&col_idx) {
            continue;
        }
        if header == name {
            return Some(meta_idx);
        }
        meta_idx += 1;
    }
    None
}

/// The state a cell shows in the explorer.
///
/// `K` is an absent cell for which the run found evidence that did not reach
/// the binary matrix: a pseudogene balanced keeps out, sequence-only evidence
/// under strict, a candidate rejected for its context. These used to be drawn
/// as plain absences, so the evidence Recore exists to surface was invisible on
/// the page while gene_presence_absence_recore_states.csv recorded it.
fn explorer_state_char(
    pa: &crate::parser::PresenceAbsence,
    ci: ClusterIdx,
    gi: GenomeIdx,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_presence: &ahash::HashSet<(usize, usize)>,
    kept_out: &ahash::HashSet<(usize, usize)>,
) -> char {
    if pa.is_present(ci, gi) {
        'P'
    } else if rescue_map.contains_key(&(ci, gi)) {
        'R'
    } else if nsc_presence.contains(&(ci.val(), gi.val())) {
        'N'
    } else if kept_out.contains(&(ci.val(), gi.val())) {
        'K'
    } else {
        'A'
    }
}

/// Evidence columns that take few distinct values. A row stores one index into
/// a table of their combinations instead of eleven string indices.
pub(super) const EXPLORER_EVIDENCE_PROFILE: [&str; 11] = [
    "decision",
    "evidence_state",
    "matrix_action",
    "publication_status",
    "recommended_wording",
    "strand",
    "synteny",
    "orf_integrity_signal",
    "rejection_reason",
    "orthology_confidence_tier",
    "risk_flags",
];
/// Free-text evidence columns, each stored as a string-table index.
pub(super) const EXPLORER_EVIDENCE_TEXT: [&str; 6] = [
    "rescue_type",
    "contig",
    "synteny_score",
    "annotation_match",
    "annotation_cluster",
    "fragment_note",
];

/// Every evidence row, accepted or rejected, for the clusters the explorer draws.
///
/// The rows are the ones orthology_evidence.tsv is written from, so a cell's
/// explanation on the page is word for word what the file says. Each drawn
/// cluster gets one base64 string, decoded only when the cluster is selected:
/// a row count, then per row (sorted by genome) the varints
/// `genome delta, profile, text columns.., start+1, end+1, identity*10+1,
/// coverage*10+1`, where 0 means the field was empty. About 25 bytes a row, so
/// a 10,000-genome run with 280,000 evidence rows adds a few MB, not 30.
pub(super) fn explorer_evidence_json(
    evidence: &[OrthologyEvidenceRow],
    rendered_index: &HashMap<usize, usize>,
    n_rendered: usize,
) -> String {
    let mut strings: Vec<String> = vec![String::new()];
    let mut string_index: HashMap<String, u64> = HashMap::default();
    string_index.insert(String::new(), 0);
    let mut intern = |value: &str| -> u64 {
        if let Some(&idx) = string_index.get(value) {
            return idx;
        }
        let idx = strings.len() as u64;
        strings.push(value.to_string());
        string_index.insert(value.to_string(), idx);
        idx
    };
    let mut profiles: Vec<Vec<u64>> = Vec::new();
    let mut profile_index: HashMap<Vec<u64>, u64> = HashMap::default();
    let mut by_cluster: Vec<Vec<(usize, Vec<u64>)>> = vec![Vec::new(); n_rendered];
    // Empty, unparsable or non-finite numbers are stored as 0, "not given".
    let plus_one = |text: &str, scale: f64| -> u64 {
        match text.parse::<f64>() {
            Ok(v) if v.is_finite() && v >= 0.0 => (v * scale).round() as u64 + 1,
            _ => 0,
        }
    };
    for row in evidence {
        let Some(&cluster) = rendered_index.get(&row.cluster_idx.val()) else {
            continue;
        };
        let profile: Vec<u64> = EXPLORER_EVIDENCE_PROFILE
            .iter()
            .map(|column| intern(row.field(column)))
            .collect();
        let profile_id = match profile_index.get(&profile) {
            Some(&id) => id,
            None => {
                let id = profiles.len() as u64;
                profiles.push(profile.clone());
                profile_index.insert(profile, id);
                id
            }
        };
        let mut values = vec![profile_id];
        values.extend(
            EXPLORER_EVIDENCE_TEXT
                .iter()
                .map(|column| intern(row.field(column))),
        );
        values.push(plus_one(row.field("start"), 1.0));
        values.push(plus_one(row.field("end"), 1.0));
        values.push(plus_one(row.field("identity"), 10.0));
        values.push(plus_one(row.field("coverage"), 10.0));
        by_cluster[cluster].push((row.genome_idx.val(), values));
    }

    let mut out = String::from("{\"profileColumns\":");
    out.push_str(&json_string_array(
        &EXPLORER_EVIDENCE_PROFILE.map(str::to_string),
    ));
    out.push_str(",\"textColumns\":");
    out.push_str(&json_string_array(
        &EXPLORER_EVIDENCE_TEXT.map(str::to_string),
    ));
    out.push_str(",\"strings\":");
    out.push_str(&json_string_array(&strings));
    out.push_str(",\"profiles\":[");
    for (idx, profile) in profiles.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push('[');
        for (j, value) in profile.iter().enumerate() {
            if j > 0 {
                out.push(',');
            }
            out.push_str(&value.to_string());
        }
        out.push(']');
    }
    out.push_str("],\"byCluster\":[");
    let mut bytes = Vec::new();
    for (cluster, rows) in by_cluster.iter_mut().enumerate() {
        if cluster > 0 {
            out.push(',');
        }
        if rows.is_empty() {
            out.push_str("\"\"");
            continue;
        }
        rows.sort_by_key(|(genome, _)| *genome);
        bytes.clear();
        cellcodec::push_varint(&mut bytes, rows.len() as u64);
        let mut previous = 0usize;
        for (genome, values) in rows.iter() {
            cellcodec::push_varint(&mut bytes, (genome - previous) as u64);
            previous = *genome;
            for &value in values {
                cellcodec::push_varint(&mut bytes, value);
            }
        }
        out.push('"');
        out.push_str(&cellcodec::base64(&bytes));
        out.push('"');
    }
    out.push_str("]}");
    out
}

/// Cells with evidence that the corrected matrix does not show, each with the
/// matrix_action of its first evidence row. The explorer draws these as "kept
/// out", and the run readout counts them by reason; both read this one map.
fn kept_out_cells<'a>(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_presence: &ahash::HashSet<(usize, usize)>,
    evidence: &'a [OrthologyEvidenceRow],
) -> ahash::HashMap<(usize, usize), &'a str> {
    let mut out = ahash::HashMap::default();
    for row in evidence {
        let (ci, gi) = (row.cluster_idx, row.genome_idx);
        if pa.is_present(ci, gi)
            || rescue_map.contains_key(&(ci, gi))
            || nsc_presence.contains(&(ci.val(), gi.val()))
        {
            continue;
        }
        out.entry((ci.val(), gi.val()))
            .or_insert_with(|| row.field("matrix_action"));
    }
    out
}

/// Kept-out cells by reason: (matrix policy, quality floors, rejected by a guard).
fn kept_out_counts(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
    evidence: &[OrthologyEvidenceRow],
) -> (usize, usize, usize) {
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let mut counts = (0, 0, 0);
    for action in kept_out_cells(pa, rescue_map, &nsc_presence, evidence).into_values() {
        match action {
            "audit_only_policy" => counts.0 += 1,
            "audit_only_quality" => counts.1 += 1,
            _ => counts.2 += 1,
        }
    }
    counts
}

#[allow(clippy::too_many_arguments)]
fn build_explorer_data_json(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
    tree: &PhyloTree,
    reference_gi: Option<GenomeIdx>,
    evolution_rows: &HashMap<String, ExplorerEvolutionRow>,
    synteny_blocks: &HashMap<String, Vec<ExplorerSyntenyBlock>>,
    evidence: &[OrthologyEvidenceRow],
    notices: &[String],
    genome_qc: &[(String, String)],
) -> Option<(String, String)> {
    use std::fmt::Write as _;

    let tree_json = pruned_phylo_json(&tree.root)?;
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let kept_out: ahash::HashSet<(usize, usize)> =
        kept_out_cells(pa, rescue_map, &nsc_presence, evidence)
            .into_keys()
            .collect();
    let clusters_to_skip: ahash::HashSet<usize> = nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();
    let annotation_idx = metadata_column_index(pa, "Annotation");

    let mut out = String::new();
    out.push_str("{\"genomes\":[");
    for (idx, genome) in pa.genome_names().iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str(&json_script_escape(genome));
    }
    out.push_str("],\"referenceGenome\":");
    if let Some(gi) = reference_gi {
        write!(
            out,
            "{{\"index\":{},\"name\":{}}}",
            gi.val(),
            json_script_escape(pa.genome_name(gi))
        )
        .expect("write to String is infallible");
    } else {
        out.push_str("null");
    }
    out.push_str(",\"clusters\":[");

    let mut states = Vec::new();
    let mut rendered_cluster_names = Vec::new();
    let mut first_cluster = true;
    // Per-genome state totals, used for the tip-burden colour modes. Counting
    // them here costs one pass we are already making; in the browser it would
    // mean decoding every cluster before the first paint.
    // [present, rescued, NSC, kept out, absent]
    let mut genome_totals = vec![[0usize; 5]; pa.n_genomes()];
    let mut rendered_index: HashMap<usize, usize> = HashMap::new();
    for ci_val in 0..pa.n_clusters() {
        if clusters_to_skip.contains(&ci_val) {
            continue;
        }
        let ci = ClusterIdx(ci_val);
        // Collect the non-absent genomes per state instead of a dense row: the
        // codec turns these into a few bytes each, where a dense row costs one
        // character per genome for every cluster in the matrix.
        let mut sets: [Vec<u32>; cellcodec::STATE_SETS] = Default::default();
        let mut present_after = 0usize;
        let mut reference_state = None;
        for (gi_val, totals) in genome_totals.iter_mut().enumerate() {
            let gi = GenomeIdx(gi_val);
            let state = explorer_state_char(pa, ci, gi, rescue_map, &nsc_presence, &kept_out);
            if Some(gi) == reference_gi {
                reference_state = Some(state);
            }
            let slot = match state {
                'P' => Some(0),
                'R' => Some(1),
                'N' => Some(2),
                'K' => Some(3),
                _ => None,
            };
            if let Some(slot) = slot {
                sets[slot].push(gi_val as u32);
                totals[slot] += 1;
                // Kept-out evidence is evidence, not presence.
                if slot < 3 {
                    present_after += 1;
                }
            } else {
                totals[4] += 1;
            }
        }
        if !first_cluster {
            out.push(',');
        }
        first_cluster = false;
        let annotation = annotation_idx
            .and_then(|idx| pa.metadata(ci).get(idx))
            .map(String::as_str)
            .unwrap_or("");
        let prevalence = if pa.n_genomes() == 0 {
            0.0
        } else {
            present_after as f64 / pa.n_genomes() as f64
        };
        let cluster_name = pa.cluster_name(ci);
        let evolution = evolution_rows
            .get(cluster_name)
            .cloned()
            .unwrap_or_default();
        let reference_cell = reference_gi.map(|gi| pa.cell(ci, gi)).unwrap_or("");
        if let Some(reference_state) = reference_state {
            write!(
                out,
                "{{\"name\":{},\"annotation\":{},\"present\":{},\"prevalence\":{:.5},\"referenceState\":{},\"referenceAbsent\":{},\"referenceCell\":{},\"mechanism\":{},\"mechanismDetail\":{},\"support\":{},\"phylo\":{},\"wording\":{},\"referenceEventClass\":{},\"targetReferenceCluster\":{},\"targetReferenceGene\":{}}}",
                json_script_escape(cluster_name),
                json_script_escape(annotation),
                present_after,
                prevalence,
                json_script_escape(&reference_state.to_string()),
                if reference_state == 'A' { "true" } else { "false" },
                json_script_escape(reference_cell),
                json_script_escape(&evolution.mechanism),
                json_script_escape(&evolution.detail),
                json_script_escape(&evolution.support),
                json_script_escape(&evolution.phylo),
                json_script_escape(&evolution.wording),
                json_script_escape(&evolution.reference_event_class),
                json_script_escape(&evolution.target_reference_cluster),
                json_script_escape(&evolution.target_reference_gene)
            )
            .expect("write to String is infallible");
        } else {
            write!(
                out,
                "{{\"name\":{},\"annotation\":{},\"present\":{},\"prevalence\":{:.5},\"referenceState\":null,\"referenceAbsent\":false,\"referenceCell\":{},\"mechanism\":{},\"mechanismDetail\":{},\"support\":{},\"phylo\":{},\"wording\":{},\"referenceEventClass\":{},\"targetReferenceCluster\":{},\"targetReferenceGene\":{}}}",
                json_script_escape(cluster_name),
                json_script_escape(annotation),
                present_after,
                prevalence,
                json_script_escape(reference_cell),
                json_script_escape(&evolution.mechanism),
                json_script_escape(&evolution.detail),
                json_script_escape(&evolution.support),
                json_script_escape(&evolution.phylo),
                json_script_escape(&evolution.wording),
                json_script_escape(&evolution.reference_event_class),
                json_script_escape(&evolution.target_reference_cluster),
                json_script_escape(&evolution.target_reference_gene)
            )
            .expect("write to String is infallible");
        }
        states.push(cellcodec::base64(&cellcodec::encode_cluster(
            &sets,
            pa.n_genomes(),
        )));
        rendered_index.insert(ci_val, rendered_cluster_names.len());
        rendered_cluster_names.push(cluster_name.to_string());
    }

    write!(out, "],\"nscEnabled\":{}", nsc_result.is_some()).expect("write to String");
    out.push_str(",\"cellEncoding\":\"recore.cells.v2\",\"cells\":[");
    for (idx, row) in states.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str(&json_script_escape(row));
    }
    out.push_str("],\"genomeStats\":[");
    for (idx, totals) in genome_totals.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        write!(
            out,
            "[{},{},{},{},{}]",
            totals[0], totals[1], totals[2], totals[3], totals[4]
        )
        .expect("write to String is infallible");
    }
    out.push_str("],\"evidence\":");
    out.push_str(&explorer_evidence_json(
        evidence,
        &rendered_index,
        states.len(),
    ));
    out.push_str(",\"notices\":");
    out.push_str(&json_string_array(notices));
    out.push_str(",\"genomeQc\":[");
    for (idx, (status, flags)) in genome_qc.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        write!(
            out,
            "[{},{}]",
            json_script_escape(status),
            json_script_escape(flags)
        )
        .expect("write to String is infallible");
    }
    out.push(']');
    out.push_str(",\"tree\":");
    out.push_str(&tree_json);
    // Synteny travels in its own script element, parsed the first time the
    // panel opens: at 10,000 x 30,000 it was most of the payload, and every
    // load paid for parsing it. The page only needs each cluster's block count
    // up front.
    out.push_str(",\"syntenyCounts\":[");
    for (idx, cluster_name) in rendered_cluster_names.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        let n = synteny_blocks.get(cluster_name).map_or(0, Vec::len);
        write!(out, "{n}").expect("write to String is infallible");
    }
    out.push_str("]}");

    let mut synteny = String::from("{");
    let mut first_synteny_cluster = true;
    for cluster_name in rendered_cluster_names {
        let Some(blocks) = synteny_blocks.get(&cluster_name) else {
            continue;
        };
        if blocks.is_empty() {
            continue;
        }
        if !first_synteny_cluster {
            synteny.push(',');
        }
        first_synteny_cluster = false;
        synteny.push_str(&json_script_escape(&cluster_name));
        synteny.push_str(":[");
        for (idx, block) in blocks.iter().enumerate() {
            if idx > 0 {
                synteny.push(',');
            }
            let carriers = cellcodec::base64(&cellcodec::encode_one_set(
                &block.carrier_genomes,
                pa.n_genomes(),
            ));
            write!(
                synteny,
                "{{\"rank\":{},\"nGenomes\":{},\"fraction\":{},\"interpretation\":{},\"examples\":{},\"carriers\":\"{}\",\"representativeGenome\":{},\"representativeContig\":{},\"representativeStart\":{},\"representativeEnd\":{},\"representativeLength\":{},\"genes\":{}",
                block.rank,
                block.n_genomes,
                json_script_escape(&block.fraction),
                json_script_escape(&block.interpretation),
                json_script_escape(&block.examples),
                carriers,
                json_script_escape(&block.representative_genome),
                json_script_escape(&block.representative_contig),
                json_script_escape(&block.representative_start),
                json_script_escape(&block.representative_end),
                json_script_escape(&block.representative_length),
                json_synteny_gene_array(&block.genes),
            )
            .expect("write to String is infallible");
            // The flanking names are only drawn when there are no gene records,
            // and the reference context only for the first block.
            if block.genes.is_empty() {
                write!(
                    synteny,
                    ",\"upstream\":{},\"downstream\":{}",
                    json_string_array(&block.upstream),
                    json_string_array(&block.downstream)
                )
                .expect("write to String is infallible");
            }
            if idx == 0 {
                write!(
                    synteny,
                    ",\"referenceGenes\":{}",
                    json_synteny_gene_array(&block.reference_genes)
                )
                .expect("write to String is infallible");
            }
            synteny.push('}');
        }
        synteny.push(']');
    }
    synteny.push('}');
    Some((out, synteny))
}

/// The workspace markup, and whether the viewer was actually embedded. Without a
/// tree (or above the size ceilings) only the fallback notice is written, and no
/// script comes with it, so the chrome that drives the viewer must not be drawn
/// either.
struct ExplorerWorkspace {
    html: String,
    embedded: bool,
}

fn explorer_unavailable_html(reason: &str) -> String {
    format!(
        "<section class='panel full explorer explorer-fallback'><h2>Phylogeny and matrix explorer</h2><p class='muted'>{}</p></section>",
        html_escape(reason)
    )
}

/// Inline 16px icon set. A dense tool reads as controls, not as a wall of words,
/// and every icon button carries a data-tip so nothing is a guess.
fn icon(name: &str) -> String {
    let path = match name {
        "book" => "<path d='M4 19.5A2.5 2.5 0 0 1 6.5 17H20'/><path d='M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z'/>",
        "github" => "<path d='M9 19c-5 1.5-5-2.5-7-3m14 6v-3.9a3.4 3.4 0 0 0-.9-2.6c3-.3 6.2-1.5 6.2-6.7A5.2 5.2 0 0 0 20 5.8a4.9 4.9 0 0 0-.1-3.6s-1.1-.3-3.6 1.4a12.3 12.3 0 0 0-6.6 0C7.2 1.9 6.1 2.2 6.1 2.2A4.9 4.9 0 0 0 6 5.8a5.2 5.2 0 0 0-1.4 3.6c0 5.2 3.2 6.4 6.2 6.7a3.4 3.4 0 0 0-.9 2.6V22'/>",
        "info" => "<circle cx='12' cy='12' r='10'/><path d='M12 16v-4M12 8h.01'/>",
        "moon" => "<path d='M21 12.79A9 9 0 1 1 11.21 3a7 7 0 0 0 9.79 9.79z'/>",
        "download" => "<path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3'/>",
        "expand" => "<path d='M8 3H5a2 2 0 0 0-2 2v3M21 8V5a2 2 0 0 0-2-2h-3M16 21h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'/>",
        "filter" => "<path d='M22 3H2l8 9.5V19l4 2v-8.5L22 3z'/>",
        "target" => "<circle cx='12' cy='12' r='10'/><circle cx='12' cy='12' r='6'/><circle cx='12' cy='12' r='2'/>",
        "layers" => "<path d='M12 2 2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5'/>",
        "tree" => "<path d='M4 4v16M4 8h6M4 16h6M10 5v6M10 13v6M10 5h5M10 11h5M10 13h5M10 19h5'/>",
        "search" => "<circle cx='11' cy='11' r='7'/><path d='m20 20-3.5-3.5'/>",
        "sliders" => "<path d='M4 6h10M18 6h2M4 12h4M12 12h8M4 18h10M18 18h2'/><circle cx='16' cy='6' r='2'/><circle cx='10' cy='12' r='2'/><circle cx='16' cy='18' r='2'/>",
        "left" => "<path d='m14 6-6 6 6 6'/>",
        "right" => "<path d='m10 6 6 6-6 6'/>",
        "keyboard" => "<rect x='2' y='6' width='20' height='12' rx='2'/><path d='M6 10h.01M10 10h.01M14 10h.01M18 10h.01M8 14h8'/>",
        "close" => "<path d='M6 6l12 12M18 6 6 18'/>",
        _ => "<circle cx='12' cy='12' r='9'/>",
    };
    // A String rather than Box::leak: this is called once per chrome control per
    // report, and leaking to satisfy a &'static str signature buys nothing.
    format!("<svg class='glyph' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round' aria-hidden='true'>{path}</svg>")
}

/// The workspace is the report: phylogeny, matrix and synteny fill the viewport
/// and everything else lives in the run-details drawer.
#[allow(clippy::too_many_arguments)]
fn explorer_workspace_html(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
    phylo_tree: Option<&PhyloTree>,
    reference_gi: Option<GenomeIdx>,
    evolution_rows: &HashMap<String, ExplorerEvolutionRow>,
    synteny_blocks: &HashMap<String, Vec<ExplorerSyntenyBlock>>,
    evidence: &[OrthologyEvidenceRow],
    genome_qc: &[(String, String)],
) -> ExplorerWorkspace {
    use std::fmt::Write as _;

    let unavailable = |reason: &str| ExplorerWorkspace {
        html: explorer_unavailable_html(reason),
        embedded: false,
    };

    let mut notices: Vec<String> = Vec::new();
    let star;
    let tree = match phylo_tree {
        Some(tree) if matched_phylo_tips(&tree.root) > 0 => tree,
        supplied => {
            star = star_tree(pa);
            notices.push(if supplied.is_some() {
                "The supplied tree has no tips matching matrix genomes, so rows follow the matrix order and the tree pane only lists the genomes.".to_string()
            } else {
                "No --tree was given, so rows follow the matrix order and the tree pane only lists the genomes. Ancestral readings need a tree.".to_string()
            });
            &star
        }
    };
    let matched_tips = matched_phylo_tips(&tree.root);
    if matched_tips > EXPLORER_MAX_TIPS {
        return unavailable(&format!(
            "The tree has {} matched tips. The workspace is embedded up to {} tips to keep the report responsive.",
            format_count(matched_tips),
            format_count(EXPLORER_MAX_TIPS)
        ));
    }

    let rendered_clusters = nsc_result
        .as_ref()
        .map(|r| pa.n_clusters().saturating_sub(r.merge_map.len()))
        .unwrap_or_else(|| pa.n_clusters());
    let cell_count = matched_tips.saturating_mul(rendered_clusters);
    if cell_count > EXPLORER_MAX_CELLS {
        return unavailable(&format!(
            "The matrix would embed {} cells ({} tips x {} clusters), above the {} cell limit. Use gene_presence_absence_recore_states.csv for the full matrix.",
            format_count(cell_count),
            format_count(matched_tips),
            format_count(rendered_clusters),
            format_count(EXPLORER_MAX_CELLS)
        ));
    }

    let build = |synteny: &HashMap<String, Vec<ExplorerSyntenyBlock>>, notices: &[String]| {
        build_explorer_data_json(
            pa,
            rescue_map,
            nsc_result,
            tree,
            reference_gi,
            evolution_rows,
            synteny,
            evidence,
            notices,
            genome_qc,
        )
    };
    let Some((mut data_json, mut synteny_json)) = build(synteny_blocks, &notices) else {
        return unavailable(
            "The tree could not be converted into a matched-tip browser workspace.",
        );
    };
    let total = data_json.len() + synteny_json.len();
    if total > EXPLORER_MAX_JSON_BYTES && !synteny_blocks.is_empty() {
        notices.push(format!(
            "Synteny detail was left out of this report: with it the embedded data would be {} MB, more than a browser reliably holds. cluster_synteny_blocks.tsv has every block.",
            format_count(total / (1024 * 1024))
        ));
        match build(&HashMap::new(), &notices) {
            Some(json) => (data_json, synteny_json) = json,
            None => {
                return unavailable(
                    "The tree could not be converted into a matched-tip browser workspace.",
                )
            }
        }
    }
    if data_json.len() + synteny_json.len() > EXPLORER_MAX_JSON_BYTES {
        return unavailable(&format!(
            "The workspace data would be {} MB, more than a browser reliably holds. The run details drawer and the TSV files have every result.",
            format_count((data_json.len() + synteny_json.len()) / (1024 * 1024))
        ));
    }

    let mut html = String::new();

    html.push_str(
        "<div id='workspaceNotices' class='workspace-notices' role='status' hidden></div>",
    );
    html.push_str("<div id='lineageLegend' class='lineage-legend'></div>");

    // ---- Workspace grid.
    html.push_str("<main id='explorer' class='workspace synteny-hidden'>");
    html.push_str("<div id='viewerSplit' class='panes-wrap'>");
    html.push_str("<div id='viewerBody' class='panes'>");

    // Tree pane. Everything that shapes the tree lives in the tree's own menu:
    // a control belongs beside the thing it changes, not in a shared strip.
    html.push_str("<section class='pane tree-pane' id='treePane' aria-label='Phylogeny'>");
    write!(
        html,
        "<header class='pane-head'><span class='pane-name' role='heading' aria-level='2'>Phylogeny</span><small id='topTreeMeta'></small><span class='head-zoom'><input id='rowHeight' type='range' min='-3.9120' max='3.6889' step='any' value='2.0794' aria-label='Row zoom' data-tip='Row height (Alt and the wheel)'></span><button id='treeMenuBtn' type='button' class='icon-btn ghost-icon' data-tip='Tree display: shape, root, labels and colour' aria-label='Tree display' aria-haspopup='dialog' aria-expanded='false' aria-controls='treeMenu'>{}</button></header>",
        icon("sliders")
    )
    .expect("write to String is infallible");
    html.push_str(
        "<div id='treeMenu' class='pane-menu' role='dialog' aria-label='Tree display' hidden>",
    );
    html.push_str("<div class='menu-row'><span class='menu-label'>Shape</span><div class='segmented-control' role='group' aria-label='Tree mode'><button id='phylogramMode' type='button' class='active' data-tip='Branch lengths drawn to scale'>Phylogram</button><button id='cladogramMode' type='button' data-tip='Topology only, branch lengths ignored'>Cladogram</button></div></div>");
    html.push_str("<div class='menu-row'><label class='menu-label' for='rootMode'>Root</label><select id='rootMode' data-tip='Where the displayed tree is rooted. Every ancestral reading depends on this.'><option value='l8'>L8 outgroup</option><option value='input'>As supplied</option><option value='midpoint'>Midpoint</option></select></div>");
    html.push_str("<p class='menu-hint'>Alt-click any tip to root the tree on that genome.</p>");
    html.push_str("<div class='menu-sep'></div>");
    html.push_str("<div class='menu-row'><label class='menu-label' for='tipColorMode'>Colour tips by</label><select id='tipColorMode'><option value='cluster'>Selected cluster</option><option value='lineage'>Lineage</option><option value='process'>Mechanism</option><option value='rescue'>Rescue burden</option><option value='absence'>Absence burden</option><option value='correction'>Any correction</option><option value='qc'>Genome QC status</option></select></div>");
    html.push_str("<div class='menu-row'><label class='menu-label' for='tipLabels'>Tip labels</label><input id='tipLabels' type='checkbox' class='switch' checked data-tip='Draw genome names beside the tips'></div>");
    html.push_str("</div>");
    // Tree, lineage strip and matrix all draw only their visible window onto a
    // viewport-sized canvas; the spacers carry the virtual size so the browser
    // still owns the scrollbars.
    html.push_str("<div class='tree-head-spacer' id='treeHeadSpacer' aria-hidden='true'></div>");
    html.push_str("<div class='pane-viewport tree-viewport'><canvas id='treeCanvas' class='tree-canvas' role='img' aria-label='Matched phylogeny'></canvas><canvas id='lineageStrip' class='lineage-strip' aria-hidden='true'></canvas><canvas id='treeOverlay' class='pane-overlay' aria-hidden='true'></canvas><div class='pane-scroll' id='treeScroll'><div class='pane-spacer' id='treeSpacer'></div></div></div>");
    html.push_str("<div id='treeScaleBar' class='scale-bar' aria-hidden='true'></div>");
    html.push_str("</section>");

    html.push_str("<div class='split-handle tree-handle' id='treeHandle' role='separator' aria-orientation='vertical' tabindex='0' aria-label='Resize phylogeny pane' title='Drag to resize'></div>");

    // Matrix pane.
    html.push_str("<section class='pane matrix-region' aria-label='Corrected matrix'>");
    write!(
        html,
        "<header class='pane-head'><span class='pane-name' role='heading' aria-level='2'>Corrected matrix</span><small id='topMatrixMeta'></small><span id='cellLegend' class='legend' aria-label='Cell legend'><span class='key'><i class='sw p'></i>input</span><span class='key'><i class='sw r'></i>corrected</span><span class='key'><i class='sw a'></i>absent</span></span><span class='head-zoom'><input id='cellSize' type='range' min='-4.6052' max='4.4998' step='any' value='2.0794' aria-label='Column zoom' data-tip='Column width (Ctrl or Cmd and the wheel)'></span><button id='matrixMenuBtn' type='button' class='icon-btn ghost-icon' data-tip='Matrix display: cell colour, column order and filters' aria-label='Matrix display' aria-haspopup='dialog' aria-expanded='false' aria-controls='matrixMenu'>{}</button></header>",
        icon("sliders")
    )
    .expect("write to String is infallible");
    html.push_str("<div id='matrixMenu' class='pane-menu wide' role='dialog' aria-label='Matrix display' hidden>");
    html.push_str("<div class='menu-row'><label class='menu-label' for='matrixColorMode'>Colour cells by</label><select id='matrixColorMode'><option value='state'>Evidence state</option><option value='process'>Mechanism</option></select></div>");
    html.push_str("<div class='menu-row'><label class='menu-label' for='columnOrder'>Column order</label><select id='columnOrder'><option value='matrix'>Input matrix</option><option value='prevalence'>Prevalence</option><option value='mechanism'>Mechanism</option><option value='corrected'>Corrected cells</option><option value='name'>Cluster name</option></select></div>");
    html.push_str("<div class='menu-row'><label class='menu-label' for='referenceOnly'>Reference-absent only</label><input id='referenceOnly' type='checkbox' class='switch' data-tip='Show only clusters absent from the reference genome'></div>");
    html.push_str("<div class='menu-sep'></div>");
    html.push_str("<div class='menu-title'>Filter columns by mechanism</div>");
    html.push_str("<div id='processLegend' class='process-legend'></div>");
    html.push_str("</div>");
    html.push_str("<div class='matrix-header'><canvas id='columnHeaderCanvas' class='column-head'></canvas><canvas id='headerOverlay' class='pane-overlay' aria-hidden='true'></canvas><div class='pane-scroll header-scroll' id='headerScroller'><div class='pane-spacer' id='headerSpacer'></div></div></div>");
    html.push_str("<div class='pane-viewport matrix-viewport'><canvas id='matrixCanvas' role='img' aria-label='Corrected presence-absence matrix'></canvas><canvas id='matrixOverlay' class='pane-overlay' aria-hidden='true'></canvas><div class='pane-scroll' id='matrixPane' tabindex='0' role='application' aria-roledescription='matrix' aria-label='Corrected matrix. Arrow keys move a cell cursor, Shift moves a page, Home and End go to the first and last row, Enter opens the cell.'><div class='pane-spacer' id='matrixSpacer'></div></div><div id='matrixTooltip' class='matrix-tooltip'></div></div>");
    // Whole-matrix overview: at 30,000 columns the scrollbar tells you nothing
    // about where the interesting columns are. Prevalence per column does.
    html.push_str("<div class='overview' id='overviewStrip' data-tip='Column prevalence across the whole matrix. Click or drag to jump.'><canvas id='overviewCanvas'></canvas><div class='overview-window' id='overviewWindow'></div></div>");
    html.push_str("</section>");

    html.push_str("<div class='split-handle inspector-handle' id='inspectorHandle' role='separator' aria-orientation='vertical' tabindex='0' aria-label='Resize inspector' title='Drag to resize'></div>");

    // Inspector. Tabbed, because a tip click should be able to answer a question
    // about a genome without evicting what you know about the cluster.
    html.push_str("<aside class='pane inspector' id='inspectorPane' aria-label='Inspector'>");
    html.push_str("<header class='pane-head inspector-head'><div class='tab-bar' role='tablist' aria-label='Inspector'><button type='button' id='clusterTab' class='tab active' role='tab' aria-selected='true' aria-controls='clusterDetails' data-panel='cluster'>Cluster</button><button type='button' id='genomeTab' class='tab' role='tab' aria-selected='false' aria-controls='genomeDetails' data-panel='genome'>Genome</button></div><button id='closeInspector' type='button' class='icon-btn ghost-icon sheet-close' aria-label='Close the inspector' data-tip='Close (Esc)'>&times;</button></header>");
    html.push_str("<div id='clusterDetails' class='cluster-details' role='tabpanel' aria-labelledby='clusterTab'></div>");
    html.push_str("<div id='genomeDetails' class='cluster-details' role='tabpanel' aria-labelledby='genomeTab' hidden></div>");
    html.push_str("</aside>");

    html.push_str("</div></div>");

    html.push_str("<div class='split-handle synteny-handle' id='syntenyHandle' role='separator' aria-orientation='horizontal' tabindex='0' aria-label='Resize synteny panel' title='Drag to resize'></div>");

    // Synteny, full width under everything else: gene arrows need horizontal room.
    html.push_str("<section class='synteny-pane' id='syntenyPanel' aria-label='Unique syntenies'>");
    write!(
        html,
        "<header class='pane-head'><span class='pane-name' role='heading' aria-level='2'>Unique syntenies</span><small id='topSyntenyMeta'></small><small id='syntenySummary' class='synteny-summary'></small><span class='zoom-inline'><label for='syntenyZoom'>Scale</label><input id='syntenyZoom' type='range' min='0.45' max='2.4' step='0.05' value='1' aria-label='Synteny zoom'></span><button id='closeSynteny' type='button' class='icon-btn ghost-icon' data-tip='Hide the synteny panel (S)' aria-label='Hide synteny'>{}</button></header>",
        icon("close")
    )
    .expect("write to String is infallible");
    html.push_str("<div class='synteny-scroll'><svg id='syntenySvg' class='synteny-svg' role='img' aria-label='Unique synteny blocks'></svg></div>");
    html.push_str("</section>");
    html.push_str("</main>");

    html.push_str("<footer class='statusbar'><span id='matrixStatus'></span><span class='status-hint'><kbd>&#8984;K</kbd> find<span class='dot'></span><kbd>?</kbd> shortcuts</span></footer>");

    html.push_str("<script id='syntenyData' type='application/json'>");
    html.push_str(&script_safe_json(&synteny_json));
    html.push_str("</script>");
    html.push_str("<script id='explorerData' type='application/json'>");
    html.push_str(&script_safe_json(&data_json));
    html.push_str("</script><script>");
    html.push_str(EXPLORER_BOOT_JS);
    html.push_str("</script><script>try{");
    html.push_str(EXPLORER_JS);
    html.push_str("}catch(error){window.recoreExplorerFailed(error);}</script>");
    ExplorerWorkspace {
        html,
        embedded: true,
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn write_html_report(
    args: &Args,
    pa: &crate::parser::PresenceAbsence,
    stats: &SummaryStats,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
    phylo_tree: Option<&PhyloTree>,
    evidence: &[OrthologyEvidenceRow],
    output_dir: &std::path::Path,
) -> Result<()> {
    use std::fmt::Write as _;

    let matrix_corrected_total = stats.matrix_corrected + stats.nsc_matrix_corrected;
    let rejected_total = stats.synteny_rejected + stats.paralog_rejected;
    let (kept_policy, kept_quality, _kept_rejected) =
        kept_out_counts(pa, rescue_map, nsc_result, evidence);
    let kept_total = kept_policy + kept_quality;
    let rescued = stats.dropped
        + stats.missing_cds
        + stats.split
        + stats.fragmented
        + stats.alignment_split
        + stats.start_variant
        + stats.paralog
        + stats.pseudogene;
    let core_after = count_core_after(pa, rescue_map, nsc_result);
    let soft_before = pa.soft_core_size();
    let soft_after = count_soft_core_after(pa, rescue_map, nsc_result);
    let core_delta = core_after as i64 - stats.core_before as i64;
    let soft_delta = soft_after as i64 - soft_before as i64;

    let orthology_actions =
        count_tsv_column(&output_dir.join("orthology_evidence.tsv"), "matrix_action")?;
    let orthology_tiers = count_tsv_column(
        &output_dir.join("orthology_evidence.tsv"),
        "orthology_confidence_tier",
    )?;
    let priorities = count_tsv_column(
        &output_dir.join("evolutionary_events.tsv"),
        "highlight_priority",
    )?;
    let mechanisms = count_tsv_column(
        &output_dir.join("evolutionary_events.tsv"),
        "molecular_mechanism",
    )?;
    let qc = count_tsv_column(&output_dir.join("genome_qc.tsv"), "qc_status")?;
    let top_candidates = top_candidate_rows(&output_dir.join("evolutionary_events.tsv"), 25)?;
    let explorer_evolution = explorer_evolution_rows(&output_dir.join("evolutionary_events.tsv"))?;
    let genome_index: HashMap<&str, u32> = pa
        .genome_names()
        .iter()
        .enumerate()
        .map(|(gi, name)| (name.as_str(), gi as u32))
        .collect();
    let explorer_synteny = explorer_synteny_blocks(
        &output_dir.join("cluster_synteny_blocks.tsv"),
        EXPLORER_MAX_SYNTENY_BLOCKS_PER_CLUSTER,
        &genome_index,
    )?;

    let priority_set: std::collections::BTreeSet<String> = top_candidates
        .iter()
        .map(|row| row.priority.clone())
        .filter(|p| !p.is_empty())
        .collect();

    let mut html = String::new();
    html.push_str("<!doctype html><html lang='en'><head><meta charset='utf-8'>");
    html.push_str("<meta name='viewport' content='width=device-width,initial-scale=1'>");
    html.push_str("<title>Recore run report</title>");
    write!(
        html,
        "<meta name='generator' content='Recore v{}'>",
        html_escape(env!("CARGO_PKG_VERSION"))
    )
    .expect("write to String is infallible");
    html.push_str("<meta name='description' content='Recore audit of a corrected pangenome presence/absence matrix.'>");
    // Resolve the theme before the first paint, otherwise a dark-mode reader gets a
    // full-page white flash while the rest of the body streams in.
    html.push_str(THEME_BOOTSTRAP_JS);
    // The report has to open from a results directory with no network, so every
    // asset is inlined and the font stack falls back to what the system has.
    html.push_str("<style>");
    html.push_str(REPORT_CSS);
    html.push_str("</style></head><body class='app'>");
    // A heading for assistive technology to start from, and one polite region
    // for what the keyboard cursor lands on and what an action did.
    html.push_str("<h1 class='sr-only'>Recore run report</h1><div id='liveRegion' class='sr-only' role='status' aria-live='polite'></div>");

    // Build the workspace first: without a tree the viewer is not embedded and no
    // script comes with it, so the chrome that drives the viewer must not be
    // drawn either. It used to be emitted unconditionally, which left a topbar of
    // controls that looked live and did nothing.
    let reference_gi = args
        .reference_genome
        .as_deref()
        .and_then(|name| resolve_reference_genome(pa, name));
    let genome_qc = genome_qc_rows(&output_dir.join("genome_qc.tsv"), pa)?;
    let workspace = explorer_workspace_html(
        pa,
        rescue_map,
        nsc_result,
        phylo_tree,
        reference_gi,
        &explorer_evolution,
        &explorer_synteny,
        evidence,
        &genome_qc,
    );
    let has_viewer = workspace.embedded;

    // ---- One chrome row. Everything that shapes a single pane moved into that
    // pane's own menu, which left identity, the find field, the run readout and
    // the global actions -- and gave the data back the 40px the second row took.
    html.push_str("<header class='topbar'>");
    html.push_str("<div class='topbar-brand'><span class='topbar-logo'>");
    html.push_str(&report_logo_svg());
    html.push_str("</span><span class='brand-rule' aria-hidden='true'></span><span class='topbar-title'>Corrected pangenome matrix</span></div>");

    // The find field is the primary control, so it gets primary placement and a
    // shortcut. A 30,000-option <select> was neither reachable nor fast.
    if has_viewer {
        html.push_str("<div class='topbar-find'>");
        write!(
        html,
        "<button id='clusterOmnibox' type='button' class='omnibox' aria-haspopup='dialog' aria-expanded='false' aria-controls='palette' data-tip='Find any cluster or genome by name or annotation'>{}<span class='omnibox-body'><span class='omnibox-name' id='omniboxName'>Search</span><span class='omnibox-note' id='omniboxNote'>clusters and genomes</span></span><kbd class='omnibox-key' id='omniboxKey'>&#8984;K</kbd></button>",
        icon("search")
    )
    .expect("write to String is infallible");
        write!(
        html,
        "<div class='step-group' role='group' aria-label='Step through clusters'><button id='prevCluster' type='button' class='icon-btn flush' data-tip='Previous cluster (\u{2190})' aria-label='Previous cluster'>{}</button><span id='clusterPosition' class='step-count num'></span><button id='nextCluster' type='button' class='icon-btn flush' data-tip='Next cluster (\u{2192})' aria-label='Next cluster'>{}</button></div>",
        icon("left"),
        icon("right")
    )
    .expect("write to String is infallible");
        html.push_str("</div>");
    }

    // A readout, not a row of pills: the number carries the weight and the unit
    // sits under it in small caps, so six figures read at a glance.
    html.push_str("<div class='run-readout'>");
    // The matrix draws clusters after NSC merging; say how the two counts relate
    // instead of showing 30,000 here and 29,320 over the matrix.
    let merged_by_nsc = nsc_result.as_ref().map_or(0, |r| r.merge_map.len());
    let clusters_tip = if merged_by_nsc > 0 {
        format!(
            "{} clusters after NSC merging: {} in the input, {} merged into near-identical clusters",
            format_count(pa.n_clusters() - merged_by_nsc),
            format_count(pa.n_clusters()),
            format_count(merged_by_nsc)
        )
    } else {
        "Gene clusters in the input matrix".to_string()
    };
    for (value, label, kind, tip) in [
        (
            format_count(pa.n_clusters() - merged_by_nsc),
            "clusters",
            "",
            clusters_tip.as_str(),
        ),
        (
            format_count(pa.n_genomes()),
            "genomes",
            "",
            "Genomes in the input matrix",
        ),
        (
            format_count(matrix_corrected_total),
            "corrected",
            "rescue",
            "Absent cells Recore turned into a supported presence",
        ),
        (
            // Rescue cells only: each one is a searchable gap, so the share is
            // bounded by 100%. NSC cells were never searched, and counting them
            // here printed shares like 220.0%.
            format_percent(stats.matrix_corrected, stats.searchable_gaps),
            "of gaps",
            "",
            "Share of the searchable absence cells that were corrected",
        ),
        (
            format_count(kept_total),
            "kept out",
            if kept_total > 0 { "kept" } else { "" },
            "Cells where Recore found the gene but the matrix policy or a quality floor kept it out of the corrected matrix",
        ),
        (
            format_count(rejected_total),
            "rejected",
            if rejected_total > 0 { "warn" } else { "" },
            "Candidate corrections rejected by the synteny and paralog guards",
        ),
        (
            format!("{:+}", core_delta),
            "core",
            if core_delta >= 0 { "pos" } else { "warn" },
            "Change in core genome size after correction",
        ),
    ] {
        write!(
            html,
            "<div class='readout {}' data-tip='{}'><b class='num'>{}</b><span>{}</span></div>",
            kind,
            html_escape(tip),
            html_escape(&value),
            html_escape(label)
        )
        .expect("write to String is infallible");
    }
    write!(
        html,
        "<div class='readout' data-tip='Matrix correction policy this run used'><b class='tag'>{:?}</b><span>policy</span></div>",
        stats.matrix_policy
    )
    .expect("write to String is infallible");
    html.push_str("</div>");

    html.push_str("<div class='topbar-actions'>");
    if has_viewer {
        write!(
            html,
            "<button id='toggleSynteny' type='button' class='icon-btn' data-tip='Show the flanking-gene context for the selected cluster (S)' aria-label='Show synteny' aria-pressed='false'>{}</button>",
            icon("layers")
        )
        .expect("write to String is infallible");
        write!(
            html,
            "<button id='shortcutsButton' type='button' class='icon-btn' data-tip='Keyboard shortcuts (?)' aria-label='Keyboard shortcuts'>{}</button>",
            icon("keyboard")
        )
        .expect("write to String is infallible");
        html.push_str("<span class='action-rule' aria-hidden='true'></span>");
    }
    write!(
        html,
        "<a class='version-chip' data-tip='Release notes for this version' href='{}/releases/tag/v{}' target='_blank' rel='noopener'>v{}</a>",
        RECORE_REPO_URL,
        html_escape(env!("CARGO_PKG_VERSION")),
        html_escape(env!("CARGO_PKG_VERSION"))
    )
    .expect("write to String is infallible");
    write!(
        html,
        "<a class='icon-btn' data-tip='Documentation' aria-label='Documentation' href='{}' target='_blank' rel='noopener'>{}</a>",
        RECORE_DOCS_URL,
        icon("book")
    )
    .expect("write to String is infallible");
    write!(
        html,
        "<a class='icon-btn' data-tip='Source on GitHub' aria-label='Source on GitHub' href='{}' target='_blank' rel='noopener'>{}</a>",
        RECORE_REPO_URL,
        icon("github")
    )
    .expect("write to String is infallible");
    write!(
        html,
        "<button id='openDetails' type='button' class='icon-btn' data-tip='Run details: summary, breakdowns, candidates and output files (D)' aria-label='Run details' aria-expanded='false' aria-controls='runDetails'>{}</button>",
        icon("info")
    )
    .expect("write to String is infallible");
    if has_viewer {
        write!(
            html,
            "<button id='exportView' type='button' class='icon-btn' data-tip='Download this view as PNG, and the synteny panel as SVG (E)' aria-label='Export view'>{}</button>",
            icon("download")
        )
        .expect("write to String is infallible");
    }
    write!(
        html,
        "<button id='themeToggle' type='button' class='icon-btn' data-tip='Switch between light and dark' aria-pressed='false' aria-label='Toggle colour theme'>{}<span class='label sr-only'>Dark</span></button>",
        icon("moon")
    )
    .expect("write to String is infallible");
    if has_viewer {
        write!(
            html,
            "<button id='explorerFullscreen' type='button' class='icon-btn' data-tip='Fullscreen (F)' aria-label='Fullscreen' aria-pressed='false'>{}</button>",
            icon("expand")
        )
        .expect("write to String is infallible");
    }
    html.push_str("</div>");
    html.push_str("</header>");

    // ---- Find-anything palette and shortcut sheet. Both are driven from
    // explorer.js, so neither is written when the viewer is not embedded.
    if has_viewer {
        html.push_str("<div id='palette' class='palette' hidden><div class='palette-scrim' id='paletteScrim'></div>");
        html.push_str("<div class='palette-box' role='dialog' aria-modal='true' aria-label='Find a cluster or genome'>");
        write!(
        html,
        "<div class='palette-input'>{}<input id='clusterSearch' type='search' autocomplete='off' spellcheck='false' placeholder='Search clusters and genomes' aria-label='Search clusters and genomes' aria-controls='paletteList'><kbd>esc</kbd></div>",
        icon("search")
    )
    .expect("write to String is infallible");
        html.push_str("<div class='palette-tabs' role='tablist' aria-label='Search scope'><button type='button' role='tab' data-scope='all' class='active' aria-selected='true'>All</button><button type='button' role='tab' data-scope='cluster' aria-selected='false'>Clusters</button><button type='button' role='tab' data-scope='genome' aria-selected='false'>Genomes</button></div>");
        html.push_str(
            "<ul id='paletteList' class='palette-list' role='listbox' aria-label='Results'></ul>",
        );
        html.push_str("<div class='palette-foot'><span><kbd>&#8593;</kbd><kbd>&#8595;</kbd> move</span><span><kbd>&#8629;</kbd> open</span><span id='paletteCount' class='num'></span></div>");
        html.push_str("</div></div>");

        // ---- Shortcut sheet.
        write!(
        html,
        "<div id='shortcutSheet' class='sheet' hidden><div class='sheet-scrim' id='sheetScrim'></div><div class='sheet-box' role='dialog' aria-modal='true' aria-labelledby='sheetTitle'><div class='sheet-head'><h2 id='sheetTitle'>Keyboard shortcuts</h2><button id='closeShortcuts' type='button' class='icon-btn' data-tip='Close (Esc)' aria-label='Close shortcuts'>{}</button></div><div class='sheet-grid'>",
        icon("close")
    )
    .expect("write to String is infallible");
        for (keys, action) in [
            ("\u{2318}K / Ctrl K", "Find a cluster or genome"),
            ("/", "Find a cluster or genome"),
            ("\u{2190} \u{2192}", "Previous / next cluster"),
            ("\u{21e7} \u{2190} \u{2192}", "Jump ten clusters"),
            ("S", "Show or hide the synteny panel"),
            ("D", "Open the run details"),
            ("E", "Export the view"),
            ("F", "Fullscreen"),
            ("T", "Switch theme"),
            ("C", "Phylogram or cladogram"),
            ("R", "Reroot at the midpoint"),
            ("Z", "Fit the matrix to the window"),
            ("\u{2318} wheel", "Zoom columns, under the pointer"),
            ("\u{2325} wheel", "Zoom rows, under the pointer"),
            ("\u{21e7} wheel", "Scroll sideways"),
            ("Alt-click a tip", "Root the tree on that genome"),
            ("Esc", "Close whatever is open"),
        ] {
            write!(
                html,
                "<div class='sheet-row'><kbd>{}</kbd><span>{}</span></div>",
                html_escape(keys),
                html_escape(action)
            )
            .expect("write to String is infallible");
        }
        html.push_str("</div></div></div>");
    }

    html.push_str(&workspace.html);

    // ---- Run details drawer: everything that used to push the viewer off screen.
    write!(
        html,
        "<aside id='runDetails' class='drawer' aria-hidden='true' inert><div class='drawer-head'><h2>Run details</h2><button id='closeDetails' type='button' class='icon-btn' data-tip='Close (Esc)' aria-label='Close run details'>{}</button></div><div class='drawer-body'>",
        icon("close")
    )
    .expect("write to String is infallible");

    write!(
        html,
        "<section class='drawer-section'><h3>Summary</h3><p><strong class='num'>{}</strong> binary matrix cells were corrected under the <strong>{:?}</strong> policy from <strong class='num'>{}</strong> searchable gaps. Recore found the gene in <strong class='num'>{}</strong> more cells that it kept out of the matrix, <strong class='num'>{}</strong> by the {:?} policy and <strong class='num'>{}</strong> by the quality floors, and its synteny and paralog guards rejected <strong class='num'>{}</strong> candidates. All <strong class='num'>{}</strong> accepted evidence rows are in orthology_evidence.tsv. Core changed from <strong class='num'>{}</strong> to <strong class='num'>{}</strong> clusters (<strong class='num'>{:+}</strong>); soft-core changed from <strong class='num'>{}</strong> to <strong class='num'>{}</strong> (<strong class='num'>{:+}</strong>).</p><p class='path' data-tip='{}'>Output directory: {}</p></section>",
        format_count(matrix_corrected_total),
        stats.matrix_policy,
        format_count(stats.searchable_gaps),
        format_count(kept_total),
        format_count(kept_policy),
        stats.matrix_policy,
        format_count(kept_quality),
        format_count(rejected_total),
        format_count(rescued),
        format_count(stats.core_before),
        format_count(core_after),
        core_delta,
        format_count(soft_before),
        format_count(soft_after),
        soft_delta,
        html_escape(&output_dir.display().to_string()),
        // The name, not the whole path: an absolute path printed into a PDF
        // carries the user's home or a temp directory into supplementary data.
        html_escape(
            &output_dir
                .file_name()
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_else(|| output_dir.display().to_string())
        ),
    )
    .expect("write to String is infallible");

    // ---- Using these results: what a reader needs to write the run up.
    let orf_review = evidence
        .iter()
        .filter(|row| {
            row.field("matrix_action") == "correct_presence"
                && row.field("orf_integrity_signal") == "out_of_frame_span_review"
        })
        .count();
    let mut caveats: Vec<String> = Vec::new();
    if kept_total > 0 {
        caveats.push(format!(
            "{} cells have sequence evidence but are not in the corrected matrix ({} by the {:?} policy, {} by the quality floors). They are drawn as kept out and carry publication_status held_out_policy or held_out_quality.",
            format_count(kept_total),
            format_count(kept_policy),
            stats.matrix_policy,
            format_count(kept_quality)
        ));
    }
    if orf_review > 0 {
        caveats.push(format!(
            "{} corrections span a length that is not a whole number of codons (publication_status corrected_review); check them before calling the gene intact.",
            format_count(orf_review)
        ));
    }
    if rejected_total > 0 {
        caveats.push(format!(
            "The synteny and paralog guards rejected {} candidates; they stay absent and are listed with decision rejected.",
            format_count(rejected_total)
        ));
    }
    if !args.synteny {
        caveats.push(
            "Synteny was off, so corrections rest on sequence similarity and annotation alone."
                .to_string(),
        );
    }
    if phylo_tree.is_none() {
        caveats
            .push("No tree was supplied, so no ancestral or lineage reading was made.".to_string());
    }
    let methods = methods_text(
        args,
        pa,
        stats,
        matrix_corrected_total,
        kept_total,
        rejected_total,
        merged_by_nsc,
    );
    html.push_str("<section class='drawer-section using'><h3>Using these results</h3>");
    write!(
        html,
        "<div class='methods'><p id='methodsText'>{}</p><button id='copyMethods' type='button' class='ghost-btn' data-tip='Copy this paragraph for a methods section'>Copy methods text</button></div>",
        html_escape(&methods)
    )
    .expect("write to String is infallible");
    if !caveats.is_empty() {
        html.push_str("<ul class='caveats'>");
        for caveat in &caveats {
            write!(html, "<li>{}</li>", html_escape(caveat))
                .expect("write to String is infallible");
        }
        html.push_str("</ul>");
    }
    html.push_str("<div class='corrections-head'><h4>Changed cells</h4><span class='corrections-actions'><button id='copyCorrections' type='button' class='ghost-btn' data-tip='Every evidence row, as tab-separated text'>Copy as TSV</button><button id='downloadCorrections' type='button' class='ghost-btn' data-tip='Every evidence row, as a TSV file'>Download TSV</button></span></div>");
    html.push_str("<div id='correctionsTable' class='corrections-table'><p class='muted'>The list of changed cells needs the interactive workspace; orthology_evidence.tsv has every row.</p></div>");
    let parameters: Vec<(&str, String)> = vec![
        (
            "Recore",
            format!("v{} ({})", env!("CARGO_PKG_VERSION"), OUTPUT_SCHEMA),
        ),
        ("Report written", utc_timestamp()),
        ("Command line", command_line()),
        (
            "Minimum identity",
            format!("{}%", trim_number(args.min_identity)),
        ),
        (
            "Minimum query coverage",
            format!("{}%", trim_number(args.min_coverage)),
        ),
        (
            "Minimum prevalence",
            format!("{}% of genomes", trim_number(args.min_prevalence * 100.0)),
        ),
        ("minimap2 preset", args.preset.clone()),
        ("Matrix policy", format!("{:?}", args.matrix_policy)),
        (
            "All alleles as queries",
            if args.multi_allele { "yes" } else { "no" }.to_string(),
        ),
        (
            "Annotations",
            if args.annotations.is_some() {
                "supplied"
            } else {
                "not supplied"
            }
            .to_string(),
        ),
        (
            "Synteny",
            if args.synteny { "on" } else { "off" }.to_string(),
        ),
        ("NSC", if args.nscluster { "on" } else { "off" }.to_string()),
        (
            "Reference genome",
            args.reference_genome
                .clone()
                .unwrap_or_else(|| "not set".to_string()),
        ),
        (
            "Tree",
            if phylo_tree.is_some() {
                "supplied"
            } else {
                "not supplied"
            }
            .to_string(),
        ),
    ];
    html.push_str("<details class='parameters'><summary>Parameters and command line</summary><dl class='provenance'>");
    for (label, value) in &parameters {
        write!(
            html,
            "<dt>{}</dt><dd>{}</dd>",
            html_escape(label),
            html_escape(value)
        )
        .expect("write to String is infallible");
    }
    html.push_str("</dl></details></section>");

    html.push_str("<section class='drawer-section'><h3>Run statistics</h3><div class='stat-grid'>");
    for (label, value, sub, accent) in [
        (
            "Searchable gaps",
            format_count(stats.searchable_gaps),
            "after filters and FASTA matching",
            "violet",
        ),
        (
            "Matrix corrected",
            format_count(matrix_corrected_total),
            "rescues + NSC merged cells",
            "pos",
        ),
        (
            "Accepted rescues",
            format_count(rescued),
            "all accepted evidence rows",
            "pos",
        ),
        (
            "Rejected evidence",
            format_count(rejected_total),
            "synteny/paralog guardrails",
            "warn",
        ),
        (
            "Core genome",
            format!(
                "{} \u{2192} {}",
                format_count(stats.core_before),
                format_count(core_after)
            ),
            "strict core clusters",
            "accent",
        ),
        (
            "Soft-core",
            format!(
                "{} \u{2192} {}",
                format_count(soft_before),
                format_count(soft_after)
            ),
            "\u{2265}95% prevalence",
            "accent",
        ),
    ] {
        write!(
            html,
            "<div class='stat {}'><div class='label'>{}</div><div class='value num'>{}</div><div class='sub'>{}</div></div>",
            accent,
            html_escape(label),
            html_escape(&value),
            html_escape(sub)
        )
        .expect("write to String is infallible");
    }
    html.push_str("</div></section>");

    html.push_str("<section class='drawer-section'><div class='section-grid'>");
    html.push_str(&count_items_html(
        "Orthology actions",
        &orthology_actions,
        "No orthology evidence rows found.",
        "",
    ));
    html.push_str(&count_items_html(
        "Confidence tiers",
        &orthology_tiers,
        "No confidence tiers found.",
        "violet",
    ));
    html.push_str(&count_items_html(
        "Evolutionary priorities",
        &priorities,
        "No evolutionary priority rows found.",
        "warn",
    ));
    html.push_str(&count_items_html(
        "Genome QC",
        &qc,
        "No genome QC rows found.",
        "qc",
    ));
    html.push_str("</div></section>");

    html.push_str("<section class='drawer-section'><h3>Dominant mechanisms</h3>");
    if mechanisms.is_empty() {
        html.push_str("<p class='muted'>No mechanism counts found.</p>");
    } else {
        let mut mechanism_rows: Vec<_> = mechanisms.iter().collect();
        mechanism_rows.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
        let max = mechanism_rows
            .iter()
            .map(|(_, count)| **count)
            .max()
            .unwrap_or(1)
            .max(1);
        // Bars are scaled to the most frequent mechanism so the ranking is readable
        // even in a long tail, and each row prints its share of all clusters so the
        // reader still sees the real magnitude behind a full-width bar.
        let total: usize = mechanisms.values().copied().sum::<usize>().max(1);
        html.push_str("<div class='mech-bars'>");
        for (label, value) in mechanism_rows.into_iter().take(12) {
            let pct = (*value as f64 / max as f64 * 100.0).clamp(1.5, 100.0);
            let share = *value as f64 / total as f64 * 100.0;
            let color = mechanism_color(label);
            write!(
                html,
                "<div class='mech-row'><div class='mech-name'><span class='swatch' style='background:{}'></span><span title='{}'>{}</span></div><div class='mech-bar-wrap'><div class='track'><div class='bar' style='width:{:.1}%;background:{}'></div></div></div><div class='mech-count num'>{} <small>{:.1}%</small></div></div>",
                color,
                html_escape(label),
                html_escape(&pretty_label(label)),
                pct,
                color,
                format_count(*value),
                share,
            )
            .expect("write to String is infallible");
        }
        html.push_str("</div>");
    }
    html.push_str("</section>");

    html.push_str("<section class='drawer-section'><h3>Top candidate clusters</h3>");
    if top_candidates.is_empty() {
        html.push_str("<p class='muted'>No evolutionary event table found.</p>");
    } else {
        html.push_str("<div class='toolbar-inline'>");
        html.push_str(
            "<input id='candidatesSearch' type='search' placeholder='Filter by cluster, wording, mechanism'>",
        );
        html.push_str(
            "<select id='candidatesPriority' aria-label='Filter by priority'><option value=''>All priorities</option>",
        );
        for prio in &priority_set {
            write!(
                html,
                "<option value='{}'>{}</option>",
                html_escape(prio),
                html_escape(&plain_label(prio))
            )
            .expect("write to String is infallible");
        }
        html.push_str("</select>");
        html.push_str("<span id='candidatesCount' class='count num'></span>");
        html.push_str("</div>");
        html.push_str("<div class='table-wrap'>");
        html.push_str("<table id='candidatesTable' data-sortable>");
        html.push_str("<thead><tr>");
        for (label, sort_type, extra_class) in [
            ("Cluster", "string", ""),
            ("Priority", "string", ""),
            ("Score", "number", " num"),
            ("Mechanism and wording", "string", ""),
        ] {
            write!(
                html,
                "<th scope='col' class='sortable{}' data-sort-type='{}' tabindex='0' role='columnheader' aria-sort='none' title='Sort by {}'>{}</th>",
                extra_class, sort_type, label, label
            )
            .expect("write to String is infallible");
        }
        html.push_str("</tr></thead><tbody>");
        for row in top_candidates {
            let badge_class = priority_badge_class(&row.priority);
            let mech_color = mechanism_color(&row.mechanism);
            write!(
                html,
                "<tr data-priority='{}' data-cluster='{}'><td><strong>{}</strong></td><td><span class='badge {}' data-tip='{}'>{}</span></td><td class='num' data-sort-value='{:.6}'>{:.3}</td><td data-sort-value='{}'><span class='mech-badge' style='--mech-color:{}' data-tip='{}'><i></i>{}</span><span class='wording'>{}</span></td></tr>",
                html_escape(&row.priority),
                html_escape(&row.cluster),
                html_escape(&row.cluster),
                badge_class,
                html_escape(&row.priority),
                html_escape(&plain_label(&row.priority)),
                row.score,
                row.score,
                html_escape(&row.mechanism),
                mech_color,
                html_escape(&format!(
                    "{}; phylogenetic pattern: {}; support: {}",
                    row.mechanism,
                    pretty_label(&row.phylo),
                    row.support
                )),
                html_escape(&pretty_label(&row.mechanism)),
                html_escape(&row.wording)
            )
            .expect("write to String is infallible");
        }
        html.push_str("</tbody></table></div>");
        html.push_str("<p class='muted candidates-foot'><button id='candidatesMore' type='button' class='ghost-btn' hidden>Show all</button> Click a row, or Tab to it and press Enter, to open that cluster in the matrix.</p>");
    }
    html.push_str("</section>");

    // Listed from the directory itself: a hand-written list linked an empty
    // file and left eight outputs out.
    html.push_str("<section class='drawer-section'><h3>Files in the output directory</h3><p class='muted'>The links open the files beside this report, so they work while the report sits in its output directory.</p><div class='files'>");
    for entry in output_file_entries(output_dir) {
        html.push_str(&entry);
    }
    html.push_str("</div></section>");

    html.push_str("</div></aside>");
    html.push_str("<div id='drawerScrim' class='drawer-scrim' hidden></div>");

    html.push_str("<script>");
    html.push_str(REPORT_UI_JS);
    html.push_str("</script>");
    html.push_str("</body></html>\n");

    let path = output_dir.join("recore_report.html");
    fs::write(&path, html).with_context(|| format!("Cannot write {}", path.display()))?;
    info!("HTML report written to {}", path.display());
    Ok(())
}
