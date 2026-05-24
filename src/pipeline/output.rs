use super::evolution::{PhyloNode, PhyloTree};
use super::reference_audit::resolve_reference_genome;
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
    carrier_genomes: Vec<String>,
    representative_genome: String,
    representative_gene: String,
    representative_contig: String,
    representative_start: String,
    representative_end: String,
    representative_strand: String,
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

fn parse_synteny_carriers(value: &str) -> Vec<String> {
    value
        .split(';')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(ToString::to_string)
        .collect()
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
    for genome in incoming.carrier_genomes {
        if !existing.carrier_genomes.contains(&genome) {
            existing.carrier_genomes.push(genome);
        }
    }
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
    if !incoming.examples.is_empty() {
        if existing.examples.is_empty() {
            existing.examples = incoming.examples;
        } else {
            existing.examples.push(';');
            existing.examples.push_str(&incoming.examples);
        }
    }
}

fn explorer_synteny_blocks(
    path: &std::path::Path,
    limit_per_cluster: usize,
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
    let representative_gene_idx = idx("representative_gene");
    let representative_contig_idx = idx("representative_contig");
    let representative_start_idx = idx("representative_start");
    let representative_end_idx = idx("representative_end");
    let representative_strand_idx = idx("representative_strand");
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
                .unwrap_or("")
                .to_string(),
            carrier_genomes: carrier_genomes_idx
                .and_then(|i| record.get(i))
                .map(parse_synteny_carriers)
                .unwrap_or_default(),
            representative_genome: representative_genome_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            representative_gene: representative_gene_idx
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
            representative_strand: representative_strand_idx
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
        let pct = (*value as f64 / total as f64 * 100.0).clamp(1.5, 100.0);
        write!(
            html,
            "<div class='bar-row'><div class='bar-label'><span>{}</span><strong class='num'>{}</strong></div><div class='track'><div class='bar {}' style='width:{:.1}%'></div></div></div>",
            html_escape(label),
            format_count(*value),
            bar_kind,
            pct
        )
        .expect("write to String is infallible");
    }
    html.push_str("</div></section>");
    html
}

fn file_link_html(path: &str, description: &str) -> String {
    format!(
        "<a class='file' href='{0}'><span class='file-name'>{0}</span><small>{1}</small></a>",
        html_escape(path),
        html_escape(description)
    )
}

const REPORT_LOGO: &str = include_str!("../../assets/recore_logo.svg");
const REPORT_CSS: &str = include_str!("report_assets/report.css");
const REPORT_UI_JS: &str = include_str!("report_assets/report-ui.js");
const EXPLORER_JS: &str = include_str!("report_assets/explorer.js");
const EXPLORER_MAX_TIPS: usize = 5_000;
const EXPLORER_MAX_CELLS: usize = 25_000_000;
const EXPLORER_MAX_SYNTENY_BLOCKS_PER_CLUSTER: usize = 10;

fn mechanism_color(name: &str) -> &'static str {
    match name {
        "clade_associated_gain" | "clean_insertion_or_acquisition_candidate" => "#079455",
        "ancestral_retained_no_disruption" => "#155eef",
        "duplication_or_paralog_expansion" => "#7a5af8",
        "deletion_loss_candidate" => "#b42318",
        "reference_or_lineage_loss_candidate" => "#d92d20",
        "IS_interruption_supported" => "#f79009",
        "pseudogenization" => "#b54708",
        "gene_fragmentation_or_structural_break" => "#dc6803",
        "partial_reference_homology" => "#0e9384",
        "reference_gene_allelic_derivative" => "#2970ff",
        "reference_gene_disruption_or_truncation" => "#c11574",
        "translocation_or_context_shift_candidate" => "#9e77ed",
        "annotation_or_clustering_artifact" => "#667085",
        "unstable_or_homoplastic_accessory_signal" => "#dd2590",
        "unresolved_accessory_variation" => "#475467",
        _ => "#475467",
    }
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

fn genome_bitset_hex<I>(indices: I, n_genomes: usize) -> String
where
    I: IntoIterator<Item = usize>,
{
    let mut nibbles = vec![0u8; (n_genomes + 3) / 4];
    for idx in indices {
        if idx < n_genomes {
            nibbles[idx / 4] |= 1 << (idx % 4);
        }
    }
    let mut out = String::with_capacity(nibbles.len());
    for nibble in nibbles {
        let ch = std::char::from_digit(nibble as u32, 16).unwrap_or('0');
        out.push(ch);
    }
    out
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
        let length = node
            .branch_length
            .map(|value| value.to_string())
            .unwrap_or_else(|| "null".to_string());
        write!(
            out,
            "{{\"label\":{},\"genome\":{},\"length\":{},\"children\":[]}}",
            json_script_escape(node.label.as_deref().unwrap_or("")),
            gi.val(),
            length
        )
        .expect("write to String is infallible");
        return Some(out);
    }

    let children: Vec<String> = node.children.iter().filter_map(pruned_phylo_json).collect();
    if children.is_empty() {
        return None;
    }

    let mut out = String::new();
    let length = node
        .branch_length
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string());
    write!(
        out,
        "{{\"label\":{},\"genome\":null,\"length\":{},\"children\":[{}]}}",
        json_script_escape(node.label.as_deref().unwrap_or("")),
        length,
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
    reference_gi: Option<GenomeIdx>,
    evolution_rows: &HashMap<String, ExplorerEvolutionRow>,
    synteny_blocks: &HashMap<String, Vec<ExplorerSyntenyBlock>>,
) -> Option<String> {
    use std::fmt::Write as _;

    let tree_json = pruned_phylo_json(&tree.root)?;
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let clusters_to_skip: ahash::HashSet<usize> = nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();
    let annotation_idx = metadata_column_index(pa, "Annotation");
    let genome_name_to_idx: HashMap<&str, usize> = pa
        .genome_names()
        .iter()
        .enumerate()
        .map(|(idx, name)| (name.as_str(), idx))
        .collect();

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
    for ci_val in 0..pa.n_clusters() {
        if clusters_to_skip.contains(&ci_val) {
            continue;
        }
        let ci = ClusterIdx(ci_val);
        let mut row_states = String::with_capacity(pa.n_genomes());
        let mut present_after = 0usize;
        let mut reference_state = None;
        for gi_val in 0..pa.n_genomes() {
            let gi = GenomeIdx(gi_val);
            let state = explorer_state_char(pa, ci, gi, rescue_map, &nsc_presence);
            if Some(gi) == reference_gi {
                reference_state = Some(state);
            }
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
        states.push(row_states);
        rendered_cluster_names.push(cluster_name.to_string());
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
    out.push_str(",\"synteny\":{");
    let mut first_synteny_cluster = true;
    for cluster_name in rendered_cluster_names {
        let Some(blocks) = synteny_blocks.get(&cluster_name) else {
            continue;
        };
        if blocks.is_empty() {
            continue;
        }
        if !first_synteny_cluster {
            out.push(',');
        }
        first_synteny_cluster = false;
        out.push_str(&json_script_escape(&cluster_name));
        out.push_str(":[");
        for (idx, block) in blocks.iter().enumerate() {
            if idx > 0 {
                out.push(',');
            }
            let carrier_bitset = genome_bitset_hex(
                block
                    .carrier_genomes
                    .iter()
                    .filter_map(|genome| genome_name_to_idx.get(genome.as_str()).copied()),
                pa.n_genomes(),
            );
            write!(
                out,
                "{{\"rank\":{},\"nGenomes\":{},\"fraction\":{},\"upstream\":{},\"downstream\":{},\"interpretation\":{},\"examples\":{},\"genomeBits\":{},\"representativeGenome\":{},\"representativeGene\":{},\"representativeContig\":{},\"representativeStart\":{},\"representativeEnd\":{},\"representativeStrand\":{},\"representativeLength\":{},\"genes\":{},\"referenceGenes\":{}}}",
                block.rank,
                block.n_genomes,
                json_script_escape(&block.fraction),
                json_string_array(&block.upstream),
                json_string_array(&block.downstream),
                json_script_escape(&block.interpretation),
                json_script_escape(&block.examples),
                json_script_escape(&carrier_bitset),
                json_script_escape(&block.representative_genome),
                json_script_escape(&block.representative_gene),
                json_script_escape(&block.representative_contig),
                json_script_escape(&block.representative_start),
                json_script_escape(&block.representative_end),
                json_script_escape(&block.representative_strand),
                json_script_escape(&block.representative_length),
                json_synteny_gene_array(&block.genes),
                json_synteny_gene_array(&block.reference_genes)
            )
            .expect("write to String is infallible");
        }
        out.push(']');
    }
    out.push('}');
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
    reference_gi: Option<GenomeIdx>,
    evolution_rows: &HashMap<String, ExplorerEvolutionRow>,
    synteny_blocks: &HashMap<String, Vec<ExplorerSyntenyBlock>>,
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

    let Some(data_json) = build_explorer_data_json(
        pa,
        rescue_map,
        nsc_result,
        tree,
        reference_gi,
        evolution_rows,
        synteny_blocks,
    ) else {
        return explorer_unavailable_html(
            "The tree could not be converted into a matched-tip browser explorer.",
        );
    };

    let mut html = String::new();
    html.push_str(
        "<section id='explorer' class='panel explorer'>",
    );
    html.push_str("<div class='explorer-head'><div><div class='kicker'>Interactive review</div><h2>Phylogeny and corrected matrix explorer</h2><p class='muted'>Click matrix columns or pick a cluster to recolor tree tips. Matrix rows follow the matched tree-tip order.</p></div><div class='legend' aria-label='Cell legend'><span><i class='sw p'></i>input</span><span><i class='sw r'></i>rescued</span><span><i class='sw n'></i>NSC</span><span><i class='sw a'></i>absent</span></div></div>");
    html.push_str("<div class='explorer-controls cluster-controls'><input id='clusterSearch' type='search' placeholder='Search cluster or annotation'><select id='clusterSelect' aria-label='Cluster selector'></select><span id='clusterInfo' class='muted num'></span></div>");
    html.push_str("<div id='clusterDetails' class='cluster-details'></div><div id='processLegend' class='process-legend'></div><div class='viewer-controls'><div class='ctrl-group'><span class='ctrl-label'>Color</span><select id='tipColorMode' aria-label='Tip color mode'><option value='cluster'>Tips: by selected cluster</option><option value='lineage'>Tips: by lineage</option><option value='process'>Tips: by selected process</option><option value='rescue'>Tips: by rescue burden</option><option value='absence'>Tips: by absence burden</option><option value='correction'>Tips: by any correction</option></select><select id='matrixColorMode' aria-label='Matrix color mode'><option value='state'>Matrix: by evidence state</option><option value='process'>Matrix: by evolutionary process</option></select></div><div class='ctrl-group'><span class='ctrl-label'>Tree</span><div class='segmented-control' role='group' aria-label='Tree mode'><button id='phylogramMode' type='button' class='active'>Phylogram</button><button id='cladogramMode' type='button'>Cladogram</button></div><label class='toggle-control'><input id='tipLabels' type='checkbox' checked>Tip labels</label></div><div class='ctrl-group'><span class='ctrl-label'>Zoom</span><label class='scale-control'>Column <input id='cellSize' type='range' min='0.25' max='18' step='0.25' value='8'></label><label class='scale-control'>Row <input id='rowHeight' type='range' min='1' max='18' value='8'></label><label class='scale-control synteny-control'>Synteny <input id='syntenyZoom' type='range' min='0.45' max='2.4' step='0.05' value='1'></label></div><div class='ctrl-group ctrl-group-right'><button id='referenceOnly' type='button' class='action-button'>Reference-absent only</button><button id='toggleSynteny' type='button' class='action-button'>Show syntenies</button><button id='explorerFullscreen' type='button' class='action-button'>Fullscreen</button></div></div><div id='lineageLegend' class='lineage-legend'></div><div id='viewerSplit' class='viewer-layout synteny-hidden'><div class='viewer-top-bar'><div class='top-label top-tree' id='topTree'><span>Phylogeny</span><small id='topTreeMeta'>matched tips</small></div><div class='top-label top-matrix'><span>Corrected matrix by gene cluster</span><small id='topMatrixMeta'></small></div><div class='top-label top-synteny synteny-control'><span>Unique syntenies</span><small id='topSyntenyMeta'>select a cluster</small></div></div><div class='viewer-body' id='viewerBody'><div class='tree-pane' id='treePane'><div class='tree-scroll' id='treeScroll'><svg id='treeSvg' role='img' aria-label='Matched phylogeny'></svg><canvas id='lineageStrip' class='lineage-strip' aria-hidden='true'></canvas></div></div><div class='split-handle tree-handle' id='treeHandle' role='separator' aria-orientation='vertical' tabindex='0' aria-label='Resize tree pane' title='Drag to resize'></div><div class='matrix-region'><div class='matrix-header' id='headerScroller'><canvas id='columnHeaderCanvas' class='column-head'></canvas></div><div class='matrix-pane' id='matrixPane'><canvas id='matrixCanvas'></canvas><div id='matrixTooltip' class='matrix-tooltip'></div></div></div><div class='split-handle synteny-handle synteny-control' id='syntenyHandle' role='separator' aria-orientation='vertical' tabindex='0' aria-label='Resize synteny pane' title='Drag to resize'></div><aside class='synteny-pane' id='syntenyPanel'><div class='synteny-head'><small id='syntenySummary'>Select a gene cluster to inspect flanking loci.</small></div><div class='synteny-scroll'><svg id='syntenySvg' class='synteny-svg' role='img' aria-label='Unique synteny blocks'></svg></div></aside></div></div><div id='matrixStatus' class='matrix-status muted'>Wheel = zoom column; Shift + wheel = zoom row. Drag the gray handles between panes to resize them.</div>");
    html.push_str("<script id='explorerData' type='application/json'>");
    html.push_str(&data_json);
    html.push_str(
        "</script><script>window.recoreExplorerData=JSON.parse(document.getElementById('explorerData').textContent);</script><script>",
    );
    html.push_str(EXPLORER_JS);
    html.push_str("</script></section>");
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
    let rejected_total = stats.synteny_rejected + stats.paralog_rejected;
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
    let top_candidates = top_candidate_rows(&output_dir.join("evolutionary_events.tsv"), 12)?;
    let explorer_evolution = explorer_evolution_rows(&output_dir.join("evolutionary_events.tsv"))?;
    let explorer_synteny = explorer_synteny_blocks(
        &output_dir.join("cluster_synteny_blocks.tsv"),
        EXPLORER_MAX_SYNTENY_BLOCKS_PER_CLUSTER,
    )?;

    let priority_set: std::collections::BTreeSet<String> = top_candidates
        .iter()
        .map(|row| row.priority.clone())
        .filter(|p| !p.is_empty())
        .collect();

    let mut html = String::new();
    html.push_str("<!doctype html><html lang='en' data-theme='light'><head><meta charset='utf-8'>");
    html.push_str("<meta name='viewport' content='width=device-width,initial-scale=1'>");
    html.push_str("<title>Recore run report</title>");
    html.push_str("<link rel='preconnect' href='https://rsms.me/'>");
    html.push_str("<link rel='stylesheet' href='https://rsms.me/inter/inter.css'>");
    html.push_str("<style>");
    html.push_str(REPORT_CSS);
    html.push_str("</style></head><body>");

    html.push_str("<div class='banner'><div class='banner-inner'><div class='banner-brand'>");
    html.push_str("<div class='banner-logo' aria-hidden='true'>");
    html.push_str(REPORT_LOGO);
    html.push_str("</div><div>");
    html.push_str("<div class='kicker'>Recore audit report</div>");
    html.push_str("<h1>Corrected pangenome matrix review</h1>");
    html.push_str("<p class='banner-tagline'>A reproducible audit of rescued false absences, guarded matrix corrections, and the downstream interpretation signals that explain them.</p>");
    write!(
        html,
        "<div class='pills'><span class='pill'><span class='dot {}'></span>synteny {}</span><span class='pill'><span class='dot {}'></span>NSC {}</span><span class='pill'><span class='dot {}'></span>tree {}</span><span class='pill'>policy <code>{:?}</code></span><span class='pill'>reference <code>{}</code></span></div>",
        if args.synteny { "" } else { "off" },
        if args.synteny { "on" } else { "off" },
        if args.nscluster { "" } else { "off" },
        if args.nscluster { "on" } else { "off" },
        if phylo_tree.is_some() { "" } else { "off" },
        if phylo_tree.is_some() { "interactive" } else { "static" },
        stats.matrix_policy,
        html_escape(args.reference_genome.as_deref().unwrap_or("not set")),
    )
    .expect("write to String is infallible");
    write!(
        html,
        "<p class='path'>Output dir: {}</p>",
        html_escape(&output_dir.display().to_string()),
    )
    .expect("write to String is infallible");
    html.push_str("</div></div><div class='banner-meta'>");
    write!(
        html,
        "<span class='version'>Recore v{}</span><span class='scale num'>{} clusters &middot; {} genomes</span><span class='scale'>schema recore.outputs.v6</span>",
        html_escape(env!("CARGO_PKG_VERSION")),
        format_count(pa.n_clusters()),
        format_count(pa.n_genomes()),
    )
    .expect("write to String is infallible");
    html.push_str("<div class='banner-actions'>");
    html.push_str("<button id='themeToggle' type='button' class='icon-btn' aria-pressed='false' title='Toggle color theme'>");
    html.push_str("<svg class='glyph' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'><path d='M21 12.79A9 9 0 1 1 11.21 3a7 7 0 0 0 9.79 9.79z'/></svg><span class='label'>Dark mode</span>");
    html.push_str("</button>");
    html.push_str("</div></div></div></div>");

    html.push_str("<nav class='nav-bar' aria-label='Section navigation'><div class='nav-inner'>");
    let mut nav_items = vec![
        ("summary", "Summary"),
        ("stats", "Stats"),
        ("breakdowns", "Breakdowns"),
        ("mechanisms", "Mechanisms"),
        ("candidates", "Candidates"),
    ];
    if phylo_tree.is_some() {
        nav_items.push(("explorer", "Explorer"));
    }
    nav_items.push(("files", "Files"));
    for (target, label) in &nav_items {
        write!(
            html,
            "<a href='#{0}' data-target='{0}'>{1}</a>",
            target, label,
        )
        .expect("write to String is infallible");
    }
    html.push_str("</div></nav>");

    html.push_str("<main>");

    html.push_str("<section id='summary' class='panel'>");
    html.push_str("<div class='executive'>");
    html.push_str("<div>");
    html.push_str("<div class='kicker'>Executive summary</div>");
    write!(
        html,
        "<h2 style='margin:6px 0 10px;font-size:18px'>Pangenome correction at a glance</h2><p><strong class='num'>{}</strong> binary matrix cells were corrected under the <strong>{:?}</strong> policy from <strong class='num'>{}</strong> searchable gaps. Recore retained <strong class='num'>{}</strong> accepted evidence rows for audit and kept <strong class='num'>{}</strong> candidates out of the binary matrix through synteny/paralog guardrails. Core changed from <strong class='num'>{}</strong> to <strong class='num'>{}</strong> clusters (<strong class='num'>{:+}</strong>); soft-core changed from <strong class='num'>{}</strong> to <strong class='num'>{}</strong> (<strong class='num'>{:+}</strong>).</p>",
        format_count(matrix_corrected_total),
        stats.matrix_policy,
        format_count(stats.searchable_gaps),
        format_count(rescued),
        format_count(rejected_total),
        format_count(stats.core_before),
        format_count(core_after),
        core_delta,
        format_count(soft_before),
        format_count(soft_after),
        soft_delta,
    )
    .expect("write to String is infallible");
    html.push_str("</div>");
    html.push_str("<div class='headline-grid'>");
    write!(
        html,
        "<div class='headline'><span>Correction rate</span><strong>{}</strong></div>",
        html_escape(&format_percent(matrix_corrected_total, stats.searchable_gaps))
    )
    .expect("write to String is infallible");
    write!(
        html,
        "<div class='headline pos'><span>Evidence retained</span><strong>{}</strong></div>",
        format_count(rescued)
    )
    .expect("write to String is infallible");
    write!(
        html,
        "<div class='headline warn'><span>Guardrail rejections</span><strong>{}</strong></div>",
        format_count(rejected_total)
    )
    .expect("write to String is infallible");
    write!(
        html,
        "<div class='headline violet'><span>Review mode</span><strong>{}</strong></div>",
        if phylo_tree.is_some() { "tree-linked" } else { "table-only" }
    )
    .expect("write to String is infallible");
    html.push_str("</div></div></section>");

    html.push_str("<section id='stats' class='panel'>");
    html.push_str("<div class='panel-head'><div><h2>Run statistics</h2><small>Headline counts behind the executive summary.</small></div></div>");
    html.push_str("<div class='stat-grid'>");
    for (label, value, sub, accent) in [
        ("Clusters", format_count(pa.n_clusters()), "input matrix", "accent"),
        ("Genomes", format_count(pa.n_genomes()), "input matrix", "accent"),
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
            format_count(stats.synteny_rejected + stats.paralog_rejected),
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

    html.push_str("<section id='breakdowns'>");
    html.push_str("<div class='section-grid'>");
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
        "pos",
    ));
    html.push_str("</div></section>");

    html.push_str("<section id='mechanisms' class='panel'>");
    html.push_str("<div class='panel-head'><div><h2>Dominant mechanisms</h2><small>Top 12 molecular mechanisms across all clusters. Click matrix cells in the explorer to filter.</small></div></div>");
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
        html.push_str("<div class='mech-bars'>");
        for (label, value) in mechanism_rows.into_iter().take(12) {
            let pct = (*value as f64 / max as f64 * 100.0).clamp(1.5, 100.0);
            let color = mechanism_color(label);
            write!(
                html,
                "<div class='mech-row'><div class='mech-name'><span class='swatch' style='background:{}'></span><span title='{}'>{}</span></div><div class='mech-bar-wrap'><div class='track'><div class='bar' style='width:{:.1}%;background:linear-gradient(90deg,{},{})'></div></div></div><div class='mech-count num'>{}</div></div>",
                color,
                html_escape(label),
                html_escape(&pretty_label(label)),
                pct,
                color,
                color,
                format_count(*value),
            )
            .expect("write to String is infallible");
        }
        html.push_str("</div>");
    }
    html.push_str("</section>");

    html.push_str("<section id='candidates' class='panel'>");
    html.push_str("<div class='panel-head'><div><h2>Top candidate clusters</h2><small>Ranked by evolutionary_signal_score. Search annotation or cluster name and sort columns.</small></div></div>");
    if top_candidates.is_empty() {
        html.push_str("<p class='muted'>No evolutionary event table found.</p>");
    } else {
        html.push_str("<div class='toolbar'>");
        html.push_str(
            "<input id='candidatesSearch' type='search' placeholder='Filter by cluster, wording, mechanism'>",
        );
        html.push_str(
            "<select id='candidatesPriority' aria-label='Filter by priority'><option value=''>All priorities</option>",
        );
        for prio in &priority_set {
            write!(
                html,
                "<option value='{0}'>{0}</option>",
                html_escape(prio)
            )
            .expect("write to String is infallible");
        }
        html.push_str("</select>");
        html.push_str("<span id='candidatesCount' class='count num'></span>");
        html.push_str("</div>");
        html.push_str("<div class='table-wrap'>");
        html.push_str("<table id='candidatesTable' data-sortable>");
        html.push_str("<thead><tr>");
        html.push_str("<th class='sortable' data-sort-type='string'>Cluster</th>");
        html.push_str("<th class='sortable' data-sort-type='string'>Priority</th>");
        html.push_str("<th class='sortable num' data-sort-type='number'>Score</th>");
        html.push_str("<th class='sortable' data-sort-type='string'>Mechanism</th>");
        html.push_str("<th class='sortable' data-sort-type='string'>Phylogeny</th>");
        html.push_str("<th class='sortable' data-sort-type='string'>Support</th>");
        html.push_str("<th>Recommended wording</th>");
        html.push_str("</tr></thead><tbody>");
        for row in top_candidates {
            let badge_class = priority_badge_class(&row.priority);
            let mech_color = mechanism_color(&row.mechanism);
            write!(
                html,
                "<tr data-priority='{}'><td><strong>{}</strong></td><td><span class='badge {}'>{}</span></td><td class='num' data-sort-value='{:.6}'>{:.3}</td><td><span class='mech-badge' style='--mech-color:{}'><i></i>{}</span></td><td>{}</td><td>{}</td><td>{}</td></tr>",
                html_escape(&row.priority),
                html_escape(&row.cluster),
                badge_class,
                html_escape(&row.priority),
                row.score,
                row.score,
                mech_color,
                html_escape(&pretty_label(&row.mechanism)),
                html_escape(&row.phylo),
                html_escape(&row.support),
                html_escape(&row.wording)
            )
            .expect("write to String is infallible");
        }
        html.push_str("</tbody></table></div>");
    }
    html.push_str("</section>");

    let reference_gi = args
        .reference_genome
        .as_deref()
        .and_then(|name| resolve_reference_genome(pa, name));
    html.push_str(&phylo_matrix_explorer_html(
        pa,
        rescue_map,
        nsc_result,
        phylo_tree,
        reference_gi,
        &explorer_evolution,
        &explorer_synteny,
    ));

    html.push_str("<section id='files' class='panel'>");
    html.push_str("<div class='panel-head'><div><h2>Files to inspect</h2><small>Authoritative outputs alongside this report.</small></div></div>");
    html.push_str("<div class='files'>");
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
        "<footer><span>Generated by <strong>Recore v{}</strong> with <code>--write-html-report</code>. This report is an audit summary; use the TSV files for exact downstream analyses.</span><span class='num'>{} clusters &middot; {} genomes</span></footer>",
        html_escape(env!("CARGO_PKG_VERSION")),
        format_count(pa.n_clusters()),
        format_count(pa.n_genomes()),
    )
    .expect("write to String is infallible");
    html.push_str("</main>");
    html.push_str("<script>");
    html.push_str(REPORT_UI_JS);
    html.push_str("</script>");
    html.push_str("</body></html>\n");

    let path = output_dir.join("recore_report.html");
    fs::write(&path, html).with_context(|| format!("Cannot write {}", path.display()))?;
    info!("HTML report written to {}", path.display());
    Ok(())
}
