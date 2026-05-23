use super::evolution::{PhyloNode, PhyloTree};
use super::*;

pub(super) fn write_corrected_csv(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    output_dir: &std::path::Path,
    nsc_result: &Option<nscluster::NsClusterResult>,
) -> Result<()> {
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
            } else if nsc_presence.contains(&(ci_val, gi_val)) {
                // Gene absent in this cluster but present in a merged NSC member
                row[col_idx].clear();
                row[col_idx].push_str("NSC_MERGED");
                isolate_count += 1;
                sequence_count += 1;
            } else {
                row[col_idx].clear();
                let cell = pa.cell(ci, gi);
                row[col_idx].push_str(cell);
                if !cell.trim().is_empty() {
                    isolate_count += 1;
                    sequence_count += cell_sequence_count(cell);
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
    Ok(())
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
    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

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
    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

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
    writeln!(
        w,
        "Remaining gaps:        {}",
        stats.searchable_gaps.saturating_sub(total_matrix_corrected)
    )?;
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
    writeln!(w, "  \"schema_version\": \"recore.outputs.v6\",")?;
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

fn count_items_html(title: &str, counts: &BTreeMap<String, usize>, empty: &str) -> String {
    use std::fmt::Write as _;
    let mut html = String::new();
    write!(
        html,
        "<section class='panel'><h2>{}</h2>",
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
    let max = counts.values().copied().max().unwrap_or(1).max(1);
    html.push_str("<div class='bars'>");
    let mut rows: Vec<_> = counts.iter().collect();
    rows.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    for (label, value) in rows {
        let pct = (*value as f64 / max as f64 * 100.0).clamp(2.0, 100.0);
        write!(
            html,
            "<div class='bar-row'><div class='bar-label'><span>{}</span><strong>{}</strong></div><div class='track'><div class='bar' style='width:{:.1}%'></div></div></div>",
            html_escape(label),
            format_count(*value),
            pct
        )
        .expect("write to String is infallible");
    }
    html.push_str("</div></section>");
    html
}

fn file_link_html(path: &str, description: &str) -> String {
    format!(
        "<a class='file' href='{0}'><span>{0}</span><small>{1}</small></a>",
        html_escape(path),
        html_escape(description)
    )
}

const REPORT_LOGO: &str = include_str!("../../assets/recore_logo.svg");
const EXPLORER_MAX_TIPS: usize = 250;
const EXPLORER_MAX_CELLS: usize = 500_000;

fn json_script_escape(value: &str) -> String {
    json_escape(value).replace("</", "<\\/")
}

fn matched_phylo_tips(node: &PhyloNode) -> usize {
    if node.children.is_empty() {
        usize::from(node.genome_idx.is_some())
    } else {
        node.children.iter().map(matched_phylo_tips).sum()
    }
}

fn pruned_phylo_json(node: &PhyloNode) -> Option<String> {
    use std::fmt::Write as _;

    if node.children.is_empty() {
        let gi = node.genome_idx?;
        let mut out = String::new();
        write!(
            out,
            "{{\"label\":{},\"genome\":{},\"children\":[]}}",
            json_script_escape(node.label.as_deref().unwrap_or("")),
            gi.val()
        )
        .expect("write to String is infallible");
        return Some(out);
    }

    let children: Vec<String> = node.children.iter().filter_map(pruned_phylo_json).collect();
    if children.is_empty() {
        return None;
    }

    let mut out = String::new();
    write!(
        out,
        "{{\"label\":{},\"genome\":null,\"children\":[{}]}}",
        json_script_escape(node.label.as_deref().unwrap_or("")),
        children.join(",")
    )
    .expect("write to String is infallible");
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

fn explorer_state_char(
    pa: &crate::parser::PresenceAbsence,
    ci: ClusterIdx,
    gi: GenomeIdx,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_presence: &ahash::HashSet<(usize, usize)>,
) -> char {
    if pa.is_present(ci, gi) {
        'P'
    } else if rescue_map.contains_key(&(ci, gi)) {
        'R'
    } else if nsc_presence.contains(&(ci.val(), gi.val())) {
        'N'
    } else {
        'A'
    }
}

fn build_explorer_data_json(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
    tree: &PhyloTree,
) -> Option<String> {
    use std::fmt::Write as _;

    let tree_json = pruned_phylo_json(&tree.root)?;
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
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
    out.push_str("],\"clusters\":[");

    let mut states = Vec::new();
    let mut first_cluster = true;
    for ci_val in 0..pa.n_clusters() {
        if clusters_to_skip.contains(&ci_val) {
            continue;
        }
        let ci = ClusterIdx(ci_val);
        let mut row_states = String::with_capacity(pa.n_genomes());
        let mut present_after = 0usize;
        for gi_val in 0..pa.n_genomes() {
            let state = explorer_state_char(pa, ci, GenomeIdx(gi_val), rescue_map, &nsc_presence);
            if state != 'A' {
                present_after += 1;
            }
            row_states.push(state);
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
        write!(
            out,
            "{{\"name\":{},\"annotation\":{},\"present\":{},\"prevalence\":{:.5}}}",
            json_script_escape(pa.cluster_name(ci)),
            json_script_escape(annotation),
            present_after,
            prevalence
        )
        .expect("write to String is infallible");
        states.push(row_states);
    }

    out.push_str("],\"states\":[");
    for (idx, state_row) in states.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str(&json_script_escape(state_row));
    }
    out.push_str("],\"tree\":");
    out.push_str(&tree_json);
    out.push('}');
    Some(out)
}

fn explorer_unavailable_html(reason: &str) -> String {
    format!(
        "<section class='panel full explorer' style='margin-top:14px'><h2>Phylogeny and matrix explorer</h2><p class='muted'>{}</p></section>",
        html_escape(reason)
    )
}

fn phylo_matrix_explorer_html(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
    phylo_tree: Option<&PhyloTree>,
) -> String {
    let Some(tree) = phylo_tree else {
        return explorer_unavailable_html(
            "Add --tree to include an interactive phylogeny plus matrix explorer in this report.",
        );
    };
    let matched_tips = matched_phylo_tips(&tree.root);
    if matched_tips == 0 {
        return explorer_unavailable_html(
            "The supplied tree has no tips matching matrix genomes, so the explorer was not embedded.",
        );
    }
    if matched_tips > EXPLORER_MAX_TIPS {
        return explorer_unavailable_html(&format!(
            "The tree has {} matched tips. The interactive explorer is embedded only up to {} tips to keep the report responsive.",
            matched_tips, EXPLORER_MAX_TIPS
        ));
    }

    let rendered_clusters = nsc_result
        .as_ref()
        .map(|r| pa.n_clusters().saturating_sub(r.merge_map.len()))
        .unwrap_or_else(|| pa.n_clusters());
    let cell_count = matched_tips.saturating_mul(rendered_clusters);
    if cell_count > EXPLORER_MAX_CELLS {
        return explorer_unavailable_html(&format!(
            "The tree and matrix would embed {} cells ({} tips x {} clusters). The explorer limit is {} cells; use gene_presence_absence_recore_states.csv for the full matrix.",
            format_count(cell_count),
            matched_tips,
            rendered_clusters,
            format_count(EXPLORER_MAX_CELLS)
        ));
    }

    let Some(data_json) = build_explorer_data_json(pa, rescue_map, nsc_result, tree) else {
        return explorer_unavailable_html(
            "The tree could not be converted into a matched-tip browser explorer.",
        );
    };

    let mut html = String::new();
    html.push_str("<section class='panel full explorer' style='margin-top:14px'>");
    html.push_str("<div class='explorer-head'><div><h2>Phylogeny and matrix explorer</h2><p class='muted'>Click matrix columns or choose a cluster to recolor tree tips dynamically. Matrix rows follow the pruned tree tip order.</p></div><div class='legend'><span><i class='sw p'></i>input</span><span><i class='sw r'></i>rescued</span><span><i class='sw n'></i>NSC</span><span><i class='sw a'></i>absent</span></div></div>");
    html.push_str("<div class='explorer-controls'><input id='clusterSearch' type='search' placeholder='Search cluster or annotation'><select id='clusterSelect'></select><span id='clusterInfo' class='muted'></span></div>");
    html.push_str("<div class='explorer-grid'><div class='tree-pane'><svg id='treeSvg' role='img' aria-label='Matched phylogeny'></svg></div><div class='matrix-pane'><canvas id='matrixCanvas'></canvas></div></div><div id='matrixStatus' class='matrix-status muted'>Hover over the matrix to inspect a cell.</div>");
    html.push_str("<script>window.recoreExplorerData=");
    html.push_str(&data_json);
    html.push_str(
        r##";(() => {
const data = window.recoreExplorerData;
if (!data || !data.clusters.length) return;
const colors = {P:"#155eef", R:"#12b76a", N:"#7a5af8", A:"#ffffff"};
const strokes = {P:"#155eef", R:"#067647", N:"#6941c6", A:"#d0d5dd"};
const labels = {P:"present in input", R:"rescued by Recore", N:"supported by NSC", A:"absent"};
const select = document.getElementById("clusterSelect");
const search = document.getElementById("clusterSearch");
const info = document.getElementById("clusterInfo");
const status = document.getElementById("matrixStatus");
const svg = document.getElementById("treeSvg");
const canvas = document.getElementById("matrixCanvas");
const ctx = canvas.getContext("2d");
let selected = 0;
let leafOrder = [];
let maxDepth = 0;
let treeHeight = 0;
let treeWidth = 360;
let cellH = 10;
let cellW = 5;

function stateAt(clusterIdx, genomeIdx) {
  return (data.states[clusterIdx] && data.states[clusterIdx][genomeIdx]) || "A";
}
function layout(node, depth) {
  node.depth = depth;
  maxDepth = Math.max(maxDepth, depth);
  if (!node.children || node.children.length === 0) {
    node.row = leafOrder.length;
    leafOrder.push(node.genome);
    return node.row;
  }
  let sum = 0;
  let count = 0;
  node.children.forEach(child => { sum += layout(child, depth + 1); count += 1; });
  node.row = count ? sum / count : leafOrder.length;
  return node.row;
}
function stateClass(state) {
  if (state === "R") return "r";
  if (state === "N") return "n";
  if (state === "A") return "a";
  return "p";
}
function populateOptions() {
  const q = search.value.trim().toLowerCase();
  select.innerHTML = "";
  data.clusters.forEach((cluster, idx) => {
    const haystack = (cluster.name + " " + (cluster.annotation || "")).toLowerCase();
    if (q && !haystack.includes(q)) return;
    const opt = document.createElement("option");
    opt.value = String(idx);
    opt.textContent = cluster.annotation ? `${cluster.name} - ${cluster.annotation}` : cluster.name;
    select.appendChild(opt);
  });
  if (![...select.options].some(opt => Number(opt.value) === selected)) {
    selected = select.options.length ? Number(select.options[0].value) : 0;
  }
  select.value = String(selected);
}
function drawNode(node, parts, xScale) {
  const x = 14 + node.depth * xScale;
  const y = 18 + node.row * cellH;
  if (node.children && node.children.length) {
    node.children.forEach(child => {
      const cx = 14 + child.depth * xScale;
      const cy = 18 + child.row * cellH;
      parts.push(`<line x1="${x}" y1="${cy}" x2="${cx}" y2="${cy}" class="tree-edge"/>`);
      parts.push(`<line x1="${x}" y1="${y}" x2="${x}" y2="${cy}" class="tree-edge"/>`);
      drawNode(child, parts, xScale);
    });
  } else {
    const state = stateAt(selected, node.genome);
    const label = data.genomes[node.genome] || node.label || "";
    parts.push(`<circle cx="${x}" cy="${y}" r="4" fill="${colors[state]}" stroke="${strokes[state]}" stroke-width="1.2"><title>${label}: ${labels[state]}</title></circle>`);
    parts.push(`<text x="${x + 9}" y="${y + 3}" class="tip-label">${escapeHtml(label)}</text>`);
  }
}
function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, ch => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[ch]));
}
function renderTree() {
  const xScale = Math.max(24, (treeWidth - 155) / Math.max(1, maxDepth));
  const parts = [];
  drawNode(data.tree, parts, xScale);
  svg.setAttribute("width", String(treeWidth));
  svg.setAttribute("height", String(treeHeight));
  svg.setAttribute("viewBox", `0 0 ${treeWidth} ${treeHeight}`);
  svg.innerHTML = parts.join("");
}
function renderMatrix(scrollToSelected) {
  const rows = leafOrder.length;
  const cols = data.clusters.length;
  cellH = rows <= 80 ? 13 : rows <= 160 ? 9 : 6;
  cellW = cols <= 120 ? 8 : cols <= 400 ? 5 : 3;
  treeHeight = Math.max(70, rows * cellH + 36);
  canvas.width = Math.max(1, cols * cellW);
  canvas.height = Math.max(1, rows * cellH);
  canvas.style.width = `${canvas.width}px`;
  canvas.style.height = `${canvas.height}px`;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  for (let r = 0; r < rows; r += 1) {
    const genomeIdx = leafOrder[r];
    for (let c = 0; c < cols; c += 1) {
      ctx.fillStyle = colors[stateAt(c, genomeIdx)] || colors.A;
      ctx.fillRect(c * cellW, r * cellH, Math.max(1, cellW - 1), Math.max(1, cellH - 1));
    }
  }
  ctx.strokeStyle = "#111827";
  ctx.lineWidth = Math.max(1, Math.min(2, cellW));
  ctx.strokeRect(selected * cellW + 0.5, 0.5, Math.max(1, cellW - 1), Math.max(1, canvas.height - 1));
  if (scrollToSelected) {
    const pane = canvas.parentElement;
    pane.scrollLeft = Math.max(0, selected * cellW - 80);
  }
  const cluster = data.clusters[selected];
  info.textContent = cluster ? `${cluster.name}: ${cluster.present}/${data.genomes.length} genomes (${Math.round(cluster.prevalence * 1000) / 10}%)` : "";
}
function render(scrollToSelected = false) {
  renderMatrix(scrollToSelected);
  renderTree();
}
layout(data.tree, 0);
populateOptions();
render(true);
search.addEventListener("input", () => { populateOptions(); render(true); });
select.addEventListener("change", () => { selected = Number(select.value || 0); render(true); });
canvas.addEventListener("mousemove", event => {
  const rect = canvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  const r = Math.floor((event.clientY - rect.top) / cellH);
  if (c < 0 || r < 0 || c >= data.clusters.length || r >= leafOrder.length) return;
  const genomeIdx = leafOrder[r];
  const state = stateAt(c, genomeIdx);
  status.textContent = `${data.clusters[c].name} / ${data.genomes[genomeIdx]}: ${labels[state]}`;
});
canvas.addEventListener("click", event => {
  const rect = canvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  if (c < 0 || c >= data.clusters.length) return;
  selected = c;
  search.value = "";
  populateOptions();
  render(true);
});
})();</script>"##,
    );
    html.push_str("</section>");
    html
}

pub(super) fn write_html_report(
    args: &Args,
    pa: &crate::parser::PresenceAbsence,
    stats: &SummaryStats,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
    phylo_tree: Option<&PhyloTree>,
    output_dir: &std::path::Path,
) -> Result<()> {
    use std::fmt::Write as _;

    let matrix_corrected_total = stats.matrix_corrected + stats.nsc_matrix_corrected;
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
    let top_candidates = top_candidate_rows(&output_dir.join("evolutionary_events.tsv"), 12)?;

    let mut html = String::new();
    html.push_str("<!doctype html><html lang='en'><head><meta charset='utf-8'>");
    html.push_str("<meta name='viewport' content='width=device-width,initial-scale=1'>");
    html.push_str("<title>Recore run report</title>");
    html.push_str(
        "<style>
        :root{color-scheme:light;--ink:#17202a;--muted:#667085;--line:#d0d5dd;--bg:#f6f8fb;--panel:#ffffff;--blue:#155eef;--green:#067647;--amber:#b54708;--red:#b42318}
        *{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 -apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif}
        main{max-width:1180px;margin:0 auto;padding:28px}
        header{background:linear-gradient(135deg,#102a43,#155eef);color:white;border-radius:18px;padding:28px;margin-bottom:18px;box-shadow:0 16px 38px rgba(16,42,67,.18)}
        .brand{display:flex;align-items:center;gap:18px}.logo{width:108px;height:62px;border-radius:14px;overflow:hidden;background:rgba(255,255,255,.94);padding:7px;flex:0 0 auto}.logo svg{width:100%;height:100%;display:block}.hero-text{min-width:0}.path{word-break:break-all}
        .kicker{text-transform:uppercase;letter-spacing:.14em;font-size:11px;opacity:.82;font-weight:700}h1{margin:6px 0 8px;font-size:30px;line-height:1.1}h2{font-size:16px;margin:0 0 14px}p{margin:0}.muted{color:var(--muted)}header .muted{color:rgba(255,255,255,.75)}
        .meta{display:flex;gap:10px;flex-wrap:wrap;margin-top:18px}.pill{border:1px solid rgba(255,255,255,.28);border-radius:999px;padding:6px 10px;background:rgba(255,255,255,.12);font-size:12px}
        .grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px}.panel,.card{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:16px;box-shadow:0 8px 20px rgba(16,24,40,.05)}
        .card .label{color:var(--muted);font-size:12px}.card .value{font-size:26px;font-weight:760;margin-top:4px}.card .sub{color:var(--muted);font-size:12px;margin-top:2px}
        .wide{grid-column:span 2}.full{grid-column:1/-1}.section-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin-top:14px}
        .bars{display:grid;gap:10px}.bar-label{display:flex;justify-content:space-between;gap:12px;font-size:12px;margin-bottom:4px}.track{height:9px;background:#eef2f6;border-radius:999px;overflow:hidden}.bar{height:100%;background:linear-gradient(90deg,#155eef,#12b76a);border-radius:999px}
        table{width:100%;border-collapse:collapse;font-size:12px}th{text-align:left;color:var(--muted);font-weight:700;border-bottom:1px solid var(--line);padding:8px 6px}td{border-bottom:1px solid #eaecf0;padding:8px 6px;vertical-align:top}tr:last-child td{border-bottom:0}
        .badge{display:inline-flex;border-radius:999px;padding:2px 8px;font-size:11px;font-weight:700;background:#eef4ff;color:#1849a9}.badge.review{background:#fff6ed;color:var(--amber)}.badge.filter{background:#fef3f2;color:var(--red)}.badge.pass{background:#ecfdf3;color:var(--green)}
        .files{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}.file{text-decoration:none;color:var(--ink);border:1px solid var(--line);border-radius:12px;padding:12px;background:#fff;display:block}.file span{font-weight:750;display:block}.file small{color:var(--muted)}
        .explorer-head{display:flex;justify-content:space-between;gap:16px;align-items:flex-start}.legend{display:flex;gap:10px;flex-wrap:wrap;font-size:12px;color:var(--muted)}.sw{display:inline-block;width:10px;height:10px;border-radius:2px;border:1px solid var(--line);margin-right:5px;vertical-align:-1px}.sw.p{background:#155eef}.sw.r{background:#12b76a}.sw.n{background:#7a5af8}.sw.a{background:#fff}
        .explorer-controls{display:grid;grid-template-columns:minmax(160px,1fr) minmax(220px,2fr) minmax(170px,1fr);gap:10px;align-items:center;margin-top:14px}.explorer input,.explorer select{border:1px solid var(--line);border-radius:10px;padding:8px 10px;background:#fff;color:var(--ink);font:inherit;min-width:0}
        .explorer-grid{display:grid;grid-template-columns:380px minmax(0,1fr);gap:12px;margin-top:14px}.tree-pane,.matrix-pane{border:1px solid var(--line);border-radius:12px;background:#fff;overflow:auto;min-height:180px;max-height:520px}.tree-edge{stroke:#98a2b3;stroke-width:1.2;fill:none}.tip-label{font-size:10px;fill:#344054}.matrix-pane{position:relative}.matrix-pane canvas{display:block;image-rendering:pixelated}.matrix-status{margin-top:8px;font-size:12px}
        footer{color:var(--muted);font-size:12px;margin:18px 4px}
        @media(max-width:900px){.grid,.section-grid,.files,.explorer-grid,.explorer-controls{grid-template-columns:1fr}.wide{grid-column:auto}main{padding:16px}h1{font-size:24px}.brand{align-items:flex-start}.logo{width:86px;height:50px}}
        </style>",
    );
    html.push_str("</head><body><main>");

    html.push_str("<header><div class='brand'><div class='logo'>");
    html.push_str(REPORT_LOGO);
    write!(
        html,
        "</div><div class='hero-text'><div class='kicker'>Recore run report</div><h1>Correction audit</h1><p class='muted'>Publication-oriented audit of the corrected pangenome matrix.</p><p class='muted path'>Output: {}</p></div></div><div class='meta'><span class='pill'>schema recore.outputs.v6</span><span class='pill'>policy {:?}</span><span class='pill'>synteny {}</span><span class='pill'>NSC {}</span><span class='pill'>reference {}</span></div></header>",
        html_escape(&output_dir.display().to_string()),
        stats.matrix_policy,
        if args.synteny { "on" } else { "off" },
        if args.nscluster { "on" } else { "off" },
        html_escape(args.reference_genome.as_deref().unwrap_or("not assessed")),
    )
    .expect("write to String is infallible");

    html.push_str("<section class='grid'>");
    for (label, value, sub) in [
        ("Clusters", format_count(pa.n_clusters()), "input matrix"),
        ("Genomes", format_count(pa.n_genomes()), "input matrix"),
        (
            "Searchable gaps",
            format_count(stats.searchable_gaps),
            "after filters and FASTA matching",
        ),
        (
            "Matrix corrected",
            format_count(matrix_corrected_total),
            "rescues + NSC merged cells",
        ),
        (
            "Accepted rescues",
            format_count(rescued),
            "all accepted evidence rows",
        ),
        (
            "Rejected evidence",
            format_count(stats.synteny_rejected + stats.paralog_rejected),
            "synteny/paralog guardrails",
        ),
        (
            "Core genome",
            format!(
                "{} -> {}",
                format_count(stats.core_before),
                format_count(core_after)
            ),
            "strict core clusters",
        ),
        (
            "Soft-core",
            format!(
                "{} -> {}",
                format_count(soft_before),
                format_count(soft_after)
            ),
            ">=95% prevalence",
        ),
    ] {
        write!(
            html,
            "<div class='card'><div class='label'>{}</div><div class='value'>{}</div><div class='sub'>{}</div></div>",
            html_escape(label),
            html_escape(&value),
            html_escape(sub)
        )
        .expect("write to String is infallible");
    }
    html.push_str("</section>");

    html.push_str("<section class='section-grid'>");
    html.push_str(&count_items_html(
        "Orthology actions",
        &orthology_actions,
        "No orthology evidence rows found.",
    ));
    html.push_str(&count_items_html(
        "Confidence tiers",
        &orthology_tiers,
        "No confidence tiers found.",
    ));
    html.push_str(&count_items_html(
        "Evolutionary priorities",
        &priorities,
        "No evolutionary priority rows found.",
    ));
    html.push_str(&count_items_html(
        "Genome QC",
        &qc,
        "No genome QC rows found.",
    ));
    html.push_str("</section>");

    html.push_str(
        "<section class='panel full' style='margin-top:14px'><h2>Dominant mechanisms</h2>",
    );
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
        html.push_str("<div class='bars'>");
        for (label, value) in mechanism_rows.into_iter().take(12) {
            let pct = (*value as f64 / max as f64 * 100.0).clamp(2.0, 100.0);
            write!(
                html,
                "<div class='bar-row'><div class='bar-label'><span>{}</span><strong>{}</strong></div><div class='track'><div class='bar' style='width:{:.1}%'></div></div></div>",
                html_escape(label),
                format_count(*value),
                pct
            )
            .expect("write to String is infallible");
        }
        html.push_str("</div>");
    }
    html.push_str("</section>");

    html.push_str(
        "<section class='panel full' style='margin-top:14px'><h2>Top candidate clusters</h2>",
    );
    if top_candidates.is_empty() {
        html.push_str("<p class='muted'>No evolutionary event table found.</p>");
    } else {
        html.push_str("<table><thead><tr><th>Cluster</th><th>Priority</th><th>Score</th><th>Mechanism</th><th>Phylogeny</th><th>Support</th><th>Recommended wording</th></tr></thead><tbody>");
        for row in top_candidates {
            let badge_class = if row.priority.contains("filter") {
                "filter"
            } else if row.priority.contains("review") || row.priority.contains("candidate") {
                "review"
            } else {
                "pass"
            };
            write!(
                html,
                "<tr><td><strong>{}</strong></td><td><span class='badge {}'>{}</span></td><td>{:.3}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                html_escape(&row.cluster),
                badge_class,
                html_escape(&row.priority),
                row.score,
                html_escape(&row.mechanism),
                html_escape(&row.phylo),
                html_escape(&row.support),
                html_escape(&row.wording)
            )
            .expect("write to String is infallible");
        }
        html.push_str("</tbody></table>");
    }
    html.push_str("</section>");

    html.push_str(&phylo_matrix_explorer_html(
        pa, rescue_map, nsc_result, phylo_tree,
    ));

    html.push_str("<section class='panel full' style='margin-top:14px'><h2>Files to inspect</h2><div class='files'>");
    for (path, description) in [
        (
            "gene_presence_absence_recore.csv",
            "Corrected binary matrix",
        ),
        (
            "gene_presence_absence_recore_states.csv",
            "Multistate evidence matrix",
        ),
        (
            "orthology_evidence.tsv",
            "Per-rescue decisions, adaptive scores, and risk flags",
        ),
        (
            "evolutionary_events.tsv",
            "Cluster-level mechanisms, priorities, and phylogenetic flags",
        ),
        (
            "reference_genome_audit.tsv",
            "Reference-absence audit when --reference-genome is set",
        ),
        (
            "genome_qc.tsv",
            "Per-genome FASTA/GFF/rescue/rejection burden",
        ),
        (
            "cluster_components.tsv",
            "Related split, alias, paralog, and NSC components",
        ),
        ("summary.txt", "Plain-text run summary"),
    ] {
        html.push_str(&file_link_html(path, description));
    }
    html.push_str("</div></section>");

    write!(
        html,
        "<footer>Generated by Recore v{} with <code>--write-html-report</code>. This report is an audit summary; use the TSV files for exact downstream analyses.</footer>",
        html_escape(env!("CARGO_PKG_VERSION"))
    )
    .expect("write to String is infallible");
    html.push_str("</main></body></html>\n");

    let path = output_dir.join("recore_report.html");
    fs::write(&path, html).with_context(|| format!("Cannot write {}", path.display()))?;
    info!("HTML report written to {}", path.display());
    Ok(())
}
