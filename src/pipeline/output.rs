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

fn explorer_synteny_blocks(
    path: &std::path::Path,
    limit_per_cluster: usize,
) -> Result<HashMap<String, Vec<ExplorerSyntenyBlock>>> {
    let mut rows: HashMap<String, Vec<ExplorerSyntenyBlock>> = HashMap::new();
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

    for record in reader.records() {
        let record = record?;
        let cluster = record.get(cluster_idx).unwrap_or("").trim();
        if cluster.is_empty() {
            continue;
        }
        let entry = rows.entry(cluster.to_string()).or_default();
        if entry.len() >= limit_per_cluster {
            continue;
        }
        entry.push(ExplorerSyntenyBlock {
            rank: rank_idx
                .and_then(|i| record.get(i))
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(entry.len() + 1),
            n_genomes: n_genomes_idx
                .and_then(|i| record.get(i))
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(0),
            fraction: fraction_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            upstream: upstream_idx
                .and_then(|i| record.get(i))
                .map(parse_synteny_side)
                .unwrap_or_default(),
            downstream: downstream_idx
                .and_then(|i| record.get(i))
                .map(parse_synteny_side)
                .unwrap_or_default(),
            interpretation: interpretation_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
            examples: examples_idx
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .to_string(),
        });
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
const EXPLORER_MAX_TIPS: usize = 5_000;
const EXPLORER_MAX_CELLS: usize = 25_000_000;
const EXPLORER_MAX_SYNTENY_BLOCKS_PER_CLUSTER: usize = 10;

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
            write!(
                out,
                "{{\"rank\":{},\"nGenomes\":{},\"fraction\":{},\"upstream\":{},\"downstream\":{},\"interpretation\":{},\"examples\":{}}}",
                block.rank,
                block.n_genomes,
                json_script_escape(&block.fraction),
                json_string_array(&block.upstream),
                json_string_array(&block.downstream),
                json_script_escape(&block.interpretation),
                json_script_escape(&block.examples)
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
        "<section id='matrixExplorer' class='panel full explorer' style='margin-top:14px'>",
    );
    html.push_str("<div class='explorer-head'><div><div class='kicker'>Interactive review</div><h2>Phylogeny and corrected matrix explorer</h2><p class='muted'>Click matrix columns or choose a cluster to recolor tree tips dynamically. Matrix rows follow the matched tree tip order.</p></div><div class='legend'><span><i class='sw p'></i>input</span><span><i class='sw r'></i>rescued</span><span><i class='sw n'></i>NSC</span><span><i class='sw a'></i>absent</span></div></div>");
    html.push_str("<div class='explorer-controls cluster-controls'><input id='clusterSearch' type='search' placeholder='Search cluster or annotation'><select id='clusterSelect'></select><span id='clusterInfo' class='muted'></span></div>");
    html.push_str("<div id='clusterDetails' class='cluster-details'></div><div id='processLegend' class='process-legend'></div><div class='viewer-controls'><select id='tipColorMode' aria-label='Tip color mode'><option value='cluster'>Color tips by selected cluster</option><option value='lineage'>Color tips by lineage</option><option value='process'>Color tips by selected process</option><option value='rescue'>Color tips by rescue burden</option><option value='absence'>Color tips by absence burden</option><option value='correction'>Color tips by any correction</option></select><select id='matrixColorMode' aria-label='Matrix color mode'><option value='state'>Matrix by evidence state</option><option value='process'>Matrix by evolutionary process</option></select><div class='segmented-control' role='group' aria-label='Tree mode'><button id='phylogramMode' type='button' class='active'>Phylogram</button><button id='cladogramMode' type='button'>Cladogram</button></div><label class='scale-control'>Column zoom <input id='cellSize' type='range' min='0.25' max='18' step='0.25' value='8'></label><label class='scale-control'>Y collapse/zoom <input id='rowHeight' type='range' min='1' max='18' value='8'></label><label class='toggle-control'><input id='tipLabels' type='checkbox' checked>Tip labels</label><button id='referenceOnly' type='button' class='action-button'>Reference-absent only</button><button id='toggleSynteny' type='button' class='action-button'>Hide syntenies</button><button id='explorerFullscreen' type='button' class='action-button'>Fullscreen explorer</button></div><div id='lineageLegend' class='lineage-legend'></div><div id='viewerSplit' class='viewer-split'><div class='viewer'><div class='viewer-top'><div class='tree-title'>Matched phylogeny</div><div class='matrix-scroller' id='headerScroller'><div class='matrix-title'>Corrected matrix by gene cluster</div><canvas id='columnHeaderCanvas' class='column-head'></canvas></div></div><div class='viewer-body' id='viewerBody'><div class='tree-pane' id='treePane'><svg id='treeSvg' role='img' aria-label='Matched phylogeny'></svg></div><div class='matrix-pane' id='matrixPane'><canvas id='matrixCanvas'></canvas><div id='matrixTooltip' class='matrix-tooltip'></div></div></div></div><aside id='syntenyPanel' class='synteny-panel'><div class='synteny-head'><span>Unique syntenies</span><small id='syntenySummary'>Select a gene cluster to inspect flanking loci.</small></div><svg id='syntenySvg' class='synteny-svg' role='img' aria-label='Unique synteny blocks'></svg></aside></div><div id='matrixStatus' class='matrix-status muted'>Wheel over tree zooms Y; wheel over the column header zooms columns. In the matrix body, Ctrl/Option + wheel zooms Y and Shift + wheel zooms columns.</div>");
    html.push_str("<script>window.recoreExplorerData=");
    html.push_str(&data_json);
    html.push_str(
        r##";(() => {
const data = window.recoreExplorerData;
if (!data || !data.clusters.length) return;
const clusterByName = new Map(data.clusters.map(cluster => [cluster.name, cluster]));
const colors = {P:"#155eef", R:"#12b76a", N:"#7a5af8", A:"#ffffff"};
const strokes = {P:"#155eef", R:"#067647", N:"#6941c6", A:"#d0d5dd"};
const labels = {P:"present in input", R:"rescued by Recore", N:"supported by NSC", A:"absent"};
const processColors = {
  clade_associated_gain:"#079455",
  ancestral_retained_no_disruption:"#155eef",
  duplication_or_paralog_expansion:"#7a5af8",
  deletion_loss_candidate:"#b42318",
  reference_or_lineage_loss_candidate:"#d92d20",
  clean_insertion_or_acquisition_candidate:"#079455",
  IS_interruption_supported:"#f79009",
  pseudogenization:"#b54708",
  gene_fragmentation_or_structural_break:"#dc6803",
  partial_reference_homology:"#0e9384",
  reference_gene_allelic_derivative:"#2970ff",
  reference_gene_disruption_or_truncation:"#c11574",
  translocation_or_context_shift_candidate:"#9e77ed",
  annotation_or_clustering_artifact:"#667085",
  unstable_or_homoplastic_accessory_signal:"#dd2590",
  unresolved_accessory_variation:"#475467",
};
const processDescriptions = {
  clade_associated_gain:"Legacy label for a clade-associated accessory gain; inspect current mechanism fields for final wording.",
  ancestral_retained_no_disruption:"Conserved ancestral presence across the assessed tree; not an obvious gain/loss event.",
  duplication_or_paralog_expansion:"Signal compatible with copy-number expansion or paralogous locus complexity.",
  deletion_loss_candidate:"Cluster is absent in a pattern compatible with deletion or lineage loss.",
  reference_or_lineage_loss_candidate:"Reference or a lineage appears to have lost a locus present elsewhere.",
  clean_insertion_or_acquisition_candidate:"Presence pattern is consistent with a gain/acquisition/insertion candidate.",
  IS_interruption_supported:"Loss/disruption is supported by nearby mobile element or insertion-sequence evidence.",
  pseudogenization:"Evidence suggests gene disruption, frameshift, internal stop, or degraded ORF.",
  gene_fragmentation_or_structural_break:"Partial homology or split alignment suggests fragmentation or a local structural break.",
  partial_reference_homology:"Cluster has partial similarity to a reference locus rather than being fully novel.",
  reference_gene_allelic_derivative:"Cluster likely represents an allelic or event-derived version of a reference gene.",
  reference_gene_disruption_or_truncation:"Reference-like sequence is present but length/frame evidence suggests disruption or truncation.",
  translocation_or_context_shift_candidate:"Sequence appears real but its genomic context differs from the reference locus.",
  annotation_or_clustering_artifact:"Pattern is more likely driven by annotation or clustering behavior than biology.",
  unstable_or_homoplastic_accessory_signal:"Patchy pattern with repeated appearances/disappearances; treat cautiously.",
  unresolved_accessory_variation:"Accessory signal remains unresolved from the available matrix/tree/evidence.",
};
const lineageColors = {
  L4:"#C94434",
  L2:"#3B5BA9",
  L3:"#8B6DAF",
  L1:"#D44D6E",
  L7:"#D4C847",
  L5:"#A07040",
  L6:"#4A9B68",
  L9:"#5CB85C",
  L10:"#7DA094",
  A1:"#C4A546",
  A2:"#6DB89B",
  A3:"#7AAED4",
  A4:"#D88EAF",
  L8:"#E08B38",
  "M. canettii":"#8C8C8C",
  unknown:"#667085",
};
const lineageDescriptions = {
  L4:"Euro-American",
  L2:"East Asian (Beijing)",
  L3:"East African-Indian (CAS)",
  L1:"Indo-Oceanic",
  L7:"Ethiopian",
  L5:"West African 1 (M. africanum)",
  L6:"West African 2 (M. africanum)",
  L9:"East African",
  L10:"Horn of Africa",
  A1:"Animal-adapted",
  A2:"Animal-adapted",
  A3:"Animal-adapted",
  A4:"Animal-adapted",
  L8:"African Great Lakes",
  "M. canettii":"Outgroup",
  unknown:"Lineage not parsed from genome name",
};
const lineageOrder = ["L4","L2","L3","L1","L7","L5","L6","L9","L10","A1","A2","A3","A4","L8","M. canettii","unknown"];
const select = document.getElementById("clusterSelect");
const search = document.getElementById("clusterSearch");
const sizeControl = document.getElementById("cellSize");
const rowHeightControl = document.getElementById("rowHeight");
const colorMode = document.getElementById("tipColorMode");
const matrixColorMode = document.getElementById("matrixColorMode");
const labelToggle = document.getElementById("tipLabels");
const referenceButton = document.getElementById("referenceOnly");
const syntenyToggle = document.getElementById("toggleSynteny");
const fullscreenButton = document.getElementById("explorerFullscreen");
const phylogramButton = document.getElementById("phylogramMode");
const cladogramButton = document.getElementById("cladogramMode");
const explorerSection = document.getElementById("matrixExplorer");
const info = document.getElementById("clusterInfo");
const details = document.getElementById("clusterDetails");
const processLegend = document.getElementById("processLegend");
const lineageLegend = document.getElementById("lineageLegend");
const status = document.getElementById("matrixStatus");
const svg = document.getElementById("treeSvg");
const headerCanvas = document.getElementById("columnHeaderCanvas");
const headerCtx = headerCanvas.getContext("2d");
const canvas = document.getElementById("matrixCanvas");
const ctx = canvas.getContext("2d");
const baseCanvas = document.createElement("canvas");
const baseCtx = baseCanvas.getContext("2d");
const matrixPane = document.getElementById("matrixPane");
const treePane = document.getElementById("treePane");
const headerScroller = document.getElementById("headerScroller");
const viewerBody = document.getElementById("viewerBody");
const viewerSplit = document.getElementById("viewerSplit");
const tooltip = document.getElementById("matrixTooltip");
const syntenyPanel = document.getElementById("syntenyPanel");
const syntenySummary = document.getElementById("syntenySummary");
const syntenySvg = document.getElementById("syntenySvg");
const defaultStatus = "Wheel over tree zooms Y; wheel over the column header zooms columns. In the matrix body, Ctrl/Option + wheel zooms Y and Shift + wheel zooms columns.";
let selected = 0;
let hoverRow = -1;
let hoverCol = -1;
let leafOrder = [];
let visibleClusters = [];
let genomeStats = [];
let maxGenomeStats = {R:1, N:1, A:1, correction:1};
let maxDepth = 0;
let maxDistance = 0;
let treeHeight = 0;
const treeWidth = 410;
const headerH = 24;
const maxCanvasDim = 32000;
const maxCanvasPixels = 90000000;
let cellH = 12;
let cellW = 8;
let baseCellH = 0;
let baseCellW = 0;
let baseValid = false;
let treeDirty = true;
let renderQueued = false;
let detailsKey = "";
let syntenyKey = "";
let treeMode = "phylogram";
let referenceOnly = false;
let syntenyHidden = false;
let mechanismFilters = new Set();
let treeRoot = data.tree;

function stateAt(clusterIdx, genomeIdx) {
  return (data.states[clusterIdx] && data.states[clusterIdx][genomeIdx]) || "A";
}
function countsForCluster(clusterIdx) {
  const counts = {P:0, R:0, N:0, A:0};
  const row = data.states[clusterIdx] || "";
  for (let i = 0; i < data.genomes.length; i += 1) {
    const state = row[i] || "A";
    counts[state] = (counts[state] || 0) + 1;
  }
  return counts;
}
function passesActiveFilters(idx, includeMechanism = true) {
  const cluster = data.clusters[idx];
  if (referenceOnly && !cluster.referenceAbsent) return false;
  if (includeMechanism && mechanismFilters.size && !mechanismFilters.has(processKind(cluster))) return false;
  return true;
}
function mechanismFilterKey() {
  return [...mechanismFilters].sort().join("|");
}
function mechanismFilterSummary() {
  const filters = [...mechanismFilters].sort();
  if (!filters.length) return "";
  const labels = filters.map(processLabel);
  if (labels.length <= 2) return labels.join(" + ");
  return `${labels.length} mechanisms`;
}
function refreshVisibleClusters() {
  visibleClusters = data.clusters
    .map((_, idx) => idx)
    .filter(idx => passesActiveFilters(idx, true));
  if (!visibleClusters.length) {
    if (referenceOnly) {
      referenceOnly = false;
      visibleClusters = data.clusters
        .map((_, idx) => idx)
        .filter(idx => passesActiveFilters(idx, true));
    }
    if (!visibleClusters.length && mechanismFilters.size) {
      mechanismFilters.clear();
      visibleClusters = data.clusters.map((_, idx) => idx);
    }
  }
  if (!visibleClusters.includes(selected)) {
    selected = visibleClusters.length ? visibleClusters[0] : 0;
  }
  if (referenceButton) {
    referenceButton.classList.toggle("active", referenceOnly);
    referenceButton.textContent = data.referenceGenome
      ? `Reference-absent only (${data.referenceGenome.name})`
      : "Reference-absent only";
    referenceButton.disabled = !data.referenceGenome;
  }
}
function selectedVisibleColumn() {
  return visibleClusters.indexOf(selected);
}
function originalClusterAtVisibleColumn(col) {
  return visibleClusters[col] ?? -1;
}
function computeGenomeStats() {
  genomeStats = data.genomes.map(() => ({P:0, R:0, N:0, A:0}));
  data.states.forEach(row => {
    for (let gi = 0; gi < data.genomes.length; gi += 1) {
      const state = row[gi] || "A";
      genomeStats[gi][state] = (genomeStats[gi][state] || 0) + 1;
    }
  });
  maxGenomeStats = genomeStats.reduce((acc, s) => ({
    R: Math.max(acc.R, s.R),
    N: Math.max(acc.N, s.N),
    A: Math.max(acc.A, s.A),
    correction: Math.max(acc.correction, s.R + s.N),
  }), {R:1, N:1, A:1, correction:1});
}
function gradientColor(value, maxValue, low, high) {
  if (!maxValue) return low;
  const t = Math.max(0, Math.min(1, value / maxValue));
  const parse = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
  const a = parse(low);
  const b = parse(high);
  const rgb = a.map((v, i) => Math.round(v + (b[i] - v) * t));
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}
function processColor(mechanism) {
  if (processColors[mechanism]) return processColors[mechanism];
  let hash = 0;
  const text = mechanism || "unresolved_accessory_variation";
  for (let i = 0; i < text.length; i += 1) hash = (hash * 31 + text.charCodeAt(i)) >>> 0;
  const hue = hash % 360;
  return `hsl(${hue} 58% 43%)`;
}
function processKind(cluster) {
  return (cluster && cluster.mechanism) || "unresolved_accessory_variation";
}
function processLabel(mechanism) {
  return (mechanism || "unresolved_accessory_variation").replace(/_/g, " ");
}
function processDescription(mechanism) {
  return processDescriptions[mechanism] || "Evolutionary process not explicitly classified; inspect TSV evidence before interpretation.";
}
function lineageOfGenomeName(name) {
  const text = String(name || "");
  if (text.startsWith("M. canettii")) return "M. canettii";
  const match = text.match(/^(L10|L[1-9]|A[1-4])(?=[_-])/);
  return match ? match[1] : "unknown";
}
function lineageColor(lineage) {
  return lineageColors[lineage] || lineageColors.unknown;
}
function lineageLabel(lineage) {
  const description = lineageDescriptions[lineage] || "";
  return description ? `${lineage} - ${description}` : lineage;
}
function renderLineageLegend() {
  if (!lineageLegend) return;
  if (colorMode.value !== "lineage") {
    lineageLegend.innerHTML = "";
    lineageLegend.style.display = "none";
    return;
  }
  const counts = new Map();
  leafOrder.forEach(genomeIdx => {
    const lineage = lineageOfGenomeName(data.genomes[genomeIdx]);
    counts.set(lineage, (counts.get(lineage) || 0) + 1);
  });
  const items = lineageOrder
    .filter(lineage => counts.has(lineage))
    .map(lineage => `<span class="lineage-chip"><i style="background:${lineageColor(lineage)}"></i><strong>${escapeHtml(lineage)}</strong><small>${counts.get(lineage)} tips</small></span>`)
    .join("");
  lineageLegend.style.display = "flex";
  lineageLegend.innerHTML = `<span class="lineage-title">MTBC lineage colors</span>${items}`;
}
function visibleProcessCounts() {
  const counts = new Map();
  data.clusters
    .map((_, idx) => idx)
    .filter(idx => passesActiveFilters(idx, false))
    .forEach(idx => {
    const mechanism = processKind(data.clusters[idx]);
    const row = counts.get(mechanism) || {clusters:0, present:0};
    row.clusters += 1;
    row.present += Number(data.clusters[idx]?.present || 0);
    counts.set(mechanism, row);
  });
  return [...counts.entries()].sort((a, b) =>
    b[1].clusters - a[1].clusters || processLabel(a[0]).localeCompare(processLabel(b[0]))
  );
}
function renderProcessLegend() {
  if (!processLegend) return;
  const cluster = data.clusters[selected] || {};
  const selectedMechanism = processKind(cluster);
  const filterLabel = mechanismFilters.size
    ? `<button type="button" class="process-clear" data-clear-mechanism="1">Showing ${escapeHtml(mechanismFilterSummary())} (${visibleClusters.length} clusters). Clear filters</button>`
    : `<span>Click mechanism cards to add/remove matrix-column filters.</span>`;
  const items = visibleProcessCounts().map(([mechanism, row]) => {
    const active = mechanism === selectedMechanism ? " active" : "";
    const filtered = mechanismFilters.has(mechanism) ? " filtered" : "";
    return `<button type="button" class="process-item${active}${filtered}" data-mechanism="${escapeHtml(mechanism)}"><span class="process-swatch" style="background:${processColor(mechanism)}"></span><div><strong>${escapeHtml(processLabel(mechanism))}</strong><small>${row.clusters} clusters; ${row.present} present cells. ${escapeHtml(processDescription(mechanism))}</small></div></button>`;
  }).join("");
  processLegend.innerHTML = `<div class="process-selected"><span class="process-chip" style="--process-color:${processColor(selectedMechanism)}"><i></i>${escapeHtml(processLabel(selectedMechanism))}</span><span>${escapeHtml(processDescription(selectedMechanism))}</span>${filterLabel}</div><div class="process-list">${items}</div>`;
}
function tipStyle(genomeIdx) {
  const mode = colorMode.value;
  if (mode === "cluster") {
    const state = stateAt(selected, genomeIdx);
    return {fill: colors[state], stroke: strokes[state], title: labels[state]};
  }
  if (mode === "process") {
    const cluster = data.clusters[selected] || {};
    const state = stateAt(selected, genomeIdx);
    const mechanism = processKind(cluster);
    const color = processColor(mechanism);
    return state === "A"
      ? {fill: "#ffffff", stroke: "#d0d5dd", title: `${processLabel(mechanism)} absent`}
      : {fill: color, stroke: color, title: `${processLabel(mechanism)}; ${labels[state]}`};
  }
  if (mode === "lineage") {
    const genome = data.genomes[genomeIdx] || "";
    const lineage = lineageOfGenomeName(genome);
    const color = lineageColor(lineage);
    return {fill: color, stroke: color, title: lineageLabel(lineage)};
  }
  const stats = genomeStats[genomeIdx] || {R:0, N:0, A:0};
  if (mode === "rescue") {
    return {fill: gradientColor(stats.R, maxGenomeStats.R, "#e0f2fe", "#0369a1"), stroke: "#0369a1", title: `rescued cells ${stats.R}`};
  }
  if (mode === "absence") {
    return {fill: gradientColor(stats.A, maxGenomeStats.A, "#f8fafc", "#b42318"), stroke: "#b42318", title: `absent cells ${stats.A}`};
  }
  const correction = stats.R + stats.N;
  return {fill: gradientColor(correction, maxGenomeStats.correction, "#ecfdf3", "#067647"), stroke: "#067647", title: `corrected/NSC cells ${correction}`};
}
function cloneNode(node) {
  return {
    label: node.label || "",
    genome: node.genome,
    length: Number(node.length || 0),
    children: (node.children || []).map(cloneNode),
  };
}
function tipName(node) {
  return node.genome === null || node.genome === undefined
    ? (node.label || "")
    : (data.genomes[node.genome] || node.label || "");
}
function isL8Tip(node) {
  return (!node.children || !node.children.length) && /^L8[_-]/.test(tipName(node));
}
function assignTreeIds(node, parent, nodes, edges) {
  const id = nodes.length;
  node._id = id;
  node._parent = parent;
  nodes.push(node);
  (node.children || []).forEach(child => {
    const childId = assignTreeIds(child, node, nodes, edges);
    const length = Number(child.length || 0);
    edges.push([id, childId, length]);
  });
  return id;
}
function collectL8Leaves(node, out) {
  if (isL8Tip(node)) out.push(node);
  (node.children || []).forEach(child => collectL8Leaves(child, out));
}
function commonAncestor(nodes) {
  if (!nodes.length) return null;
  const paths = nodes.map(node => {
    const path = [];
    let cursor = node;
    while (cursor) {
      path.push(cursor);
      cursor = cursor._parent;
    }
    return path.reverse();
  });
  let mrca = paths[0][0];
  for (let i = 0; ; i += 1) {
    const candidate = paths[0][i];
    if (!candidate || paths.some(path => path[i] !== candidate)) break;
    mrca = candidate;
  }
  return mrca;
}
function rootByL8(node) {
  const root = cloneNode(node);
  const nodes = [];
  const edges = [];
  assignTreeIds(root, null, nodes, edges);
  const leaves = [];
  collectL8Leaves(root, leaves);
  if (!leaves.length) return root;
  const rootNode = leaves.length === 1 && leaves[0]._parent ? leaves[0]._parent : commonAncestor(leaves);
  const adj = nodes.map(() => []);
  edges.forEach(([a, b, length]) => {
    adj[a].push({id:b, length});
    adj[b].push({id:a, length});
  });
  function build(id, parentId, length) {
    const source = nodes[id];
    const children = adj[id]
      .filter(edge => edge.id !== parentId)
      .map(edge => build(edge.id, id, edge.length));
    children.sort((a, b) => Number(hasL8Tip(b)) - Number(hasL8Tip(a)) || countTips(b) - countTips(a));
    return {label: source.label || "", genome: source.genome, length, children};
  }
  return build(rootNode ? rootNode._id : root._id, null, 0);
}
function hasL8Tip(node) {
  return isL8Tip(node) || (node.children || []).some(hasL8Tip);
}
function countTips(node) {
  if (!node.children || !node.children.length) return 1;
  return node.children.reduce((sum, child) => sum + countTips(child), 0);
}
function findLeafByGenome(node, genomeIdx) {
  if ((!node.children || !node.children.length) && node.genome === genomeIdx) return node;
  for (const child of node.children || []) {
    const found = findLeafByGenome(child, genomeIdx);
    if (found) return found;
  }
  return null;
}
function layout(node, depth, distance = 0) {
  node.depth = depth;
  node.distance = distance;
  maxDepth = Math.max(maxDepth, depth);
  maxDistance = Math.max(maxDistance, distance);
  if (!node.children || node.children.length === 0) {
    node.row = leafOrder.length;
    leafOrder.push(node.genome);
    return node.row;
  }
  let sum = 0;
  let count = 0;
  node.children.forEach(child => { sum += layout(child, depth + 1, distance + Number(child.length || 0)); count += 1; });
  node.row = count ? sum / count : leafOrder.length;
  return node.row;
}
function populateOptions() {
  const q = search.value.trim().toLowerCase();
  select.innerHTML = "";
  visibleClusters.forEach(idx => {
    const cluster = data.clusters[idx];
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
function currentCellSize() {
  const value = Number(sizeControl.value || 8);
  return Math.max(0.25, Math.min(18, value));
}
function currentRowHeight() {
  const value = Number(rowHeightControl.value || 8);
  return Math.max(1, Math.min(18, value));
}
function updateCellGeometry() {
  const rows = leafOrder.length;
  const cols = visibleClusters.length;
  const colSize = currentCellSize();
  const maxWByDim = Math.max(1, Math.floor(maxCanvasDim / Math.max(1, cols)));
  const maxHByDim = Math.max(1, Math.floor(maxCanvasDim / Math.max(1, rows)));
  cellH = Math.min(currentRowHeight(), maxHByDim);
  cellW = Math.min(Math.max(1, Math.round(colSize * 0.82)), maxWByDim);
  const pixelBudget = rows * cols * cellH * cellW;
  if (pixelBudget > maxCanvasPixels) {
    const scale = Math.sqrt(maxCanvasPixels / pixelBudget);
    cellH = Math.max(1, Math.floor(cellH * scale));
    cellW = Math.max(1, Math.floor(cellW * scale));
  }
  treeHeight = Math.max(70, Math.ceil(rows * cellH));
  return {rows, cols};
}
function clusterLabel(cluster) {
  return cluster.annotation ? `${cluster.name} - ${cluster.annotation}` : cluster.name;
}
function referenceInterpretation(cluster) {
  const referenceCell = (cluster.referenceCell || "").trim();
  const targetGene = (cluster.targetReferenceGene || "").trim();
  const targetCluster = (cluster.targetReferenceCluster || "").trim();
  const eventClass = (cluster.referenceEventClass || "").trim();
  const parts = [];
  if (referenceCell) {
    parts.push(`<span>reference gene(s): ${escapeHtml(referenceCell)}</span>`);
  }
  if (targetGene || targetCluster || eventClass) {
    const sourceBits = [];
    if (targetGene) sourceBits.push(`source gene ${targetGene}`);
    if (targetCluster) sourceBits.push(`source cluster ${targetCluster}`);
    if (eventClass) sourceBits.push(eventClass);
    parts.push(`<span>reference-derived/event source: ${escapeHtml(sourceBits.join("; "))}</span>`);
  }
  if (!parts.length) {
    const mechanism = processKind(cluster);
    const text = mechanism === "clean_insertion_or_acquisition_candidate"
      ? "no reference equivalent detected; candidate de novo/reference-absent locus"
      : "no reference equivalent detected by Recore";
    parts.push(`<span>${escapeHtml(text)}</span>`);
  }
  const state = cluster.referenceState ? `reference state ${cluster.referenceState}` : "reference state not assessed";
  return `<div class="reference-meta"><strong>Reference equivalent</strong>${parts.join("")}<span>${escapeHtml(state)}</span></div>`;
}
function xForNode(node) {
  const labelPad = labelToggle.checked && leafOrder.length <= 1500 ? 150 : 24;
  const usable = Math.max(80, treeWidth - 22 - labelPad);
  if (treeMode === "phylogram" && maxDistance > 0) {
    return 14 + node.distance / maxDistance * usable;
  }
  const tipX = 14 + usable;
  if (!node.children || !node.children.length) {
    return tipX;
  }
  return 14 + node.depth * (usable / Math.max(1, maxDepth));
}
function drawNode(node, parts, xScale) {
  const x = xForNode(node);
  const y = node.row * cellH + cellH / 2;
  if (node.children && node.children.length) {
    node.children.forEach(child => {
      const cx = xForNode(child);
      const cy = child.row * cellH + cellH / 2;
      const branchTitle = treeMode === "phylogram" ? `<title>branch length ${formatBranchLength(child.length)}</title>` : "";
      parts.push(`<line x1="${x}" y1="${cy}" x2="${cx}" y2="${cy}" class="tree-edge">${branchTitle}</line>`);
      parts.push(`<line x1="${x}" y1="${y}" x2="${x}" y2="${cy}" class="tree-edge"/>`);
      drawNode(child, parts, xScale);
    });
  } else {
    const style = tipStyle(node.genome);
    const label = data.genomes[node.genome] || node.label || "";
    const radius = Math.max(0.18, Math.min(5, cellH / 2 - 0.08));
    const strokeWidth = Math.max(0.2, Math.min(1.2, cellH * 0.45));
    parts.push(`<circle cx="${x}" cy="${y}" r="${radius}" fill="${style.fill}" stroke="${style.stroke}" stroke-width="${strokeWidth}"><title>${escapeHtml(label)}: ${escapeHtml(style.title)}</title></circle>`);
    if (labelToggle.checked && leafOrder.length <= 1500 && cellH >= 4) {
      parts.push(`<text x="${x + 9}" y="${y + 3}" class="tip-label">${escapeHtml(label)}</text>`);
    }
  }
}
function formatBranchLength(value) {
  const n = Number(value || 0);
  if (!Number.isFinite(n)) return "0";
  if (n === 0) return "0";
  if (Math.abs(n) < 0.001) return n.toExponential(2);
  return n.toPrecision(4);
}
function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, ch => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[ch]));
}
function clipText(text, maxLen) {
  const value = String(text || "");
  return value.length > maxLen ? `${value.slice(0, Math.max(1, maxLen - 3))}...` : value;
}
function compactAnnotation(text) {
  return String(text || "").split(/[;|]/)[0].trim();
}
function syntenyGeneLabel(name) {
  const cluster = clusterByName.get(name);
  const annotation = compactAnnotation(cluster && cluster.annotation);
  if (annotation && annotation !== name) return `${name} ${annotation}`;
  return name || "NA";
}
function firstExampleLocus(block) {
  const first = String(block.examples || "").split(";").map(item => item.trim()).find(Boolean) || "";
  if (!first) return "example locus unavailable";
  const parts = first.split(":");
  return parts.length > 1 ? parts.slice(1).join(":") : first;
}
function syntenyGenes(block, clusterName) {
  const upstream = [...(block.upstream || [])].reverse().map(name => ({name, target:false}));
  const downstream = (block.downstream || []).map(name => ({name, target:false}));
  return [...upstream, {name:clusterName, target:true}, ...downstream];
}
function arrowPath(x, y, w, h) {
  const head = Math.min(12, Math.max(7, w * 0.24));
  return `M ${x} ${y} H ${x + w - head} L ${x + w} ${y + h / 2} L ${x + w - head} ${y + h} H ${x} Z`;
}
function renderSyntenyPanel() {
  if (syntenyHidden) return;
  if (!syntenyPanel || !syntenySvg || !syntenySummary) return;
  const cluster = data.clusters[selected];
  const blocks = cluster && data.synteny ? (data.synteny[cluster.name] || []) : [];
  const nextKey = cluster ? `${cluster.name}|${blocks.length}` : "none";
  if (nextKey === syntenyKey) return;
  syntenyKey = nextKey;
  if (!cluster) {
    syntenySummary.textContent = "No cluster selected.";
    syntenySvg.setAttribute("width", "620");
    syntenySvg.setAttribute("height", "80");
    syntenySvg.innerHTML = `<text x="14" y="36" class="synteny-label">No cluster selected.</text>`;
    return;
  }
  if (!blocks.length) {
    syntenySummary.textContent = `${cluster.name}: no 5-neighbour synteny blocks in this run.`;
    syntenySvg.setAttribute("width", "620");
    syntenySvg.setAttribute("height", "88");
    syntenySvg.innerHTML = `<text x="14" y="38" class="synteny-label">No unique synteny block data for ${escapeHtml(cluster.name)}.</text>`;
    return;
  }
  const rowH = 78;
  const geneW = 62;
  const geneH = 19;
  const gap = 11;
  const left = 126;
  const top = 52;
  const maxGenes = Math.max(1, ...blocks.map(block => syntenyGenes(block, cluster.name).length));
  const width = Math.max(640, left + maxGenes * (geneW + gap) + 28);
  const height = top + blocks.length * rowH + 20;
  const parts = [];
  let previousPositions = new Map();
  syntenySummary.textContent = `${cluster.name}: ${blocks.length} unique 5-neighbour block${blocks.length === 1 ? "" : "s"} shown; selected cluster in blue.`;
  blocks.forEach((block, rowIdx) => {
    const genes = syntenyGenes(block, cluster.name);
    const y = top + rowIdx * rowH;
    const positions = new Map();
    let targetCenter = null;
    genes.forEach((gene, idx) => {
      const x = left + idx * (geneW + gap);
      const center = {x:x + geneW / 2, y:y + geneH / 2};
      if (!positions.has(gene.name)) positions.set(gene.name, center);
      if (gene.target) targetCenter = center;
    });
    positions.forEach((pos, name) => {
      const previous = previousPositions.get(name);
      if (!previous) return;
      const midY = (previous.y + pos.y) / 2;
      parts.push(`<path d="M ${previous.x} ${previous.y} C ${previous.x} ${midY} ${pos.x} ${midY} ${pos.x} ${pos.y}" class="synteny-link"><title>${escapeHtml(syntenyGeneLabel(name))}</title></path>`);
    });
    parts.push(`<text x="12" y="${y + 13}" class="synteny-row-label">#${block.rank} | n=${block.nGenomes} | f=${escapeHtml(block.fraction || "NA")}</text>`);
    const locus = firstExampleLocus(block);
    const locusY = y - 13;
    parts.push(`<text x="${left}" y="${locusY}" class="synteny-locus">${escapeHtml(clipText(locus, 58))}</text>`);
    if (targetCenter) {
      parts.push(`<line x1="${left}" y1="${locusY + 4}" x2="${targetCenter.x}" y2="${targetCenter.y - geneH / 2}" class="synteny-locus-line"/>`);
    }
    genes.forEach((gene, idx) => {
      const x = left + idx * (geneW + gap);
      const fill = gene.target ? "#155eef" : "#eef2f6";
      const cls = gene.target ? "synteny-gene synteny-target" : "synteny-gene";
      const labelClass = gene.target ? "synteny-target-label" : "synteny-label";
      const label = gene.target ? cluster.name : gene.name;
      parts.push(`<path d="${arrowPath(x, y, geneW, geneH)}" fill="${fill}" class="${cls}"><title>${escapeHtml(syntenyGeneLabel(gene.name))}</title></path>`);
      parts.push(`<text x="${x + geneW / 2}" y="${y + 13}" text-anchor="middle" class="${labelClass}">${escapeHtml(clipText(label, 10))}</text>`);
      parts.push(`<text x="${x + geneW / 2}" y="${y + geneH + 13}" text-anchor="middle" class="synteny-label">${escapeHtml(clipText(compactAnnotation((clusterByName.get(gene.name) || {}).annotation), 11))}</text>`);
    });
    const note = (block.interpretation || "").replace(/_/g, " ");
    if (note) {
      parts.push(`<text x="12" y="${y + geneH + 31}" class="synteny-label">${escapeHtml(clipText(note, 32))}</text>`);
    }
    previousPositions = positions;
  });
  syntenySvg.setAttribute("width", String(width));
  syntenySvg.setAttribute("height", String(height));
  syntenySvg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  syntenySvg.innerHTML = parts.join("");
  const firstTargetIndex = syntenyGenes(blocks[0], cluster.name).findIndex(gene => gene.target);
  if (firstTargetIndex >= 0) {
    const targetX = left + firstTargetIndex * (geneW + gap) + geneW / 2;
    syntenyPanel.scrollLeft = Math.max(0, targetX - syntenyPanel.clientWidth / 2);
  }
  syntenyPanel.scrollTop = 0;
}
function applySyntenyVisibility() {
  if (!viewerSplit || !syntenyToggle) return;
  viewerSplit.classList.toggle("synteny-hidden", syntenyHidden);
  syntenyToggle.classList.toggle("active", syntenyHidden);
  syntenyToggle.textContent = syntenyHidden ? "Show syntenies" : "Hide syntenies";
  syntenyToggle.setAttribute(
    "aria-pressed",
    syntenyHidden ? "true" : "false",
  );
  if (!syntenyHidden) {
    syntenyKey = "";
    renderSyntenyPanel();
  }
}
function renderTree() {
  const xScale = Math.max(24, (treeWidth - 155) / Math.max(1, maxDepth));
  if (treeDirty) {
    const parts = [];
    drawNode(treeRoot, parts, xScale);
    const rootLabel = hasL8Tip(treeRoot) ? "L8-rooted" : "input-rooted";
    if (treeMode === "phylogram" && maxDistance > 0) {
      parts.push(`<text x="14" y="13" class="tree-scale">phylogram, ${rootLabel}, branch length max ${escapeHtml(formatBranchLength(maxDistance))}</text>`);
    } else {
      parts.push(`<text x="14" y="13" class="tree-scale">cladogram, ${rootLabel}, branch lengths ignored</text>`);
    }
    svg.setAttribute("width", String(treeWidth));
    svg.setAttribute("height", String(treeHeight));
    svg.setAttribute("viewBox", `0 0 ${treeWidth} ${treeHeight}`);
    svg.innerHTML = `<line id="treeHoverGuide" x1="0" y1="0" x2="${treeWidth}" y2="0" class="tree-row-guide" style="display:none"/>${parts.join("")}`;
    treeDirty = false;
  }
  const guide = document.getElementById("treeHoverGuide");
  if (!guide) return;
  if (hoverRow >= 0) {
    const hy = hoverRow * cellH + cellH / 2;
    guide.setAttribute("y1", String(hy));
    guide.setAttribute("y2", String(hy));
    guide.style.display = "block";
  } else {
    guide.style.display = "none";
  }
}
function renderColumnHeader() {
  const cols = visibleClusters.length;
  const selectedCol = selectedVisibleColumn();
  headerCanvas.width = Math.max(1, cols * cellW);
  headerCanvas.height = headerH;
  headerCanvas.style.width = `${headerCanvas.width}px`;
  headerCanvas.style.height = `${headerH}px`;
  headerCtx.clearRect(0, 0, headerCanvas.width, headerCanvas.height);
  headerCtx.fillStyle = "#f8fafc";
  headerCtx.fillRect(0, 0, headerCanvas.width, headerH);
  if (matrixColorMode.value === "process") {
    visibleClusters.forEach((clusterIdx, c) => {
      headerCtx.fillStyle = processColor(processKind(data.clusters[clusterIdx]));
      headerCtx.fillRect(c * cellW, 0, Math.max(1, cellW), 7);
    });
  }
  if (hoverCol >= 0) {
    headerCtx.fillStyle = "#fff7ed";
    headerCtx.fillRect(hoverCol * cellW, 0, cellW, headerH);
  }
  if (selectedCol >= 0) {
    headerCtx.fillStyle = "#eef4ff";
    headerCtx.fillRect(selectedCol * cellW, 0, cellW, headerH);
    headerCtx.strokeStyle = "#1d4ed8";
    headerCtx.strokeRect(selectedCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), headerH - 1);
  }
  for (let c = 0; c < cols; c += 1) {
    const clusterIdx = originalClusterAtVisibleColumn(c);
    headerCtx.fillStyle = matrixColorMode.value === "process"
      ? processColor(processKind(data.clusters[clusterIdx]))
      : "#98a2b3";
    headerCtx.fillRect(c * cellW, headerH - 9, Math.max(1, cellW), 6);
  }
}
function rebuildBaseMatrix() {
  const {rows, cols} = updateCellGeometry();
  baseCanvas.width = Math.max(1, cols * cellW);
  baseCanvas.height = Math.max(1, Math.ceil(rows * cellH));
  baseCtx.clearRect(0, 0, baseCanvas.width, baseCanvas.height);
  baseCtx.fillStyle = "#ffffff";
  baseCtx.fillRect(0, 0, baseCanvas.width, baseCanvas.height);
  for (let r = 0; r < rows; r += 1) {
    const genomeIdx = leafOrder[r];
    for (let c = 0; c < cols; c += 1) {
      const clusterIdx = originalClusterAtVisibleColumn(c);
      const state = stateAt(clusterIdx, genomeIdx);
      if (matrixColorMode.value === "process" && state !== "A") {
        baseCtx.fillStyle = processColor(processKind(data.clusters[clusterIdx]));
      } else {
        baseCtx.fillStyle = colors[state] || colors.A;
      }
      const cellDrawH = cellH >= 2 ? Math.max(1, cellH - 1) : Math.max(0.2, cellH);
      baseCtx.fillRect(c * cellW, r * cellH, Math.max(1, cellW - 1), cellDrawH);
    }
  }
  baseCellH = cellH;
  baseCellW = cellW;
  baseValid = true;
}
function renderMatrix(scrollToSelected) {
  const {rows, cols} = updateCellGeometry();
  const selectedCol = selectedVisibleColumn();
  if (!baseValid || baseCellH !== cellH || baseCellW !== cellW) {
    rebuildBaseMatrix();
  }
  canvas.width = baseCanvas.width;
  canvas.height = baseCanvas.height;
  canvas.style.width = `${canvas.width}px`;
  canvas.style.height = `${canvas.height}px`;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.drawImage(baseCanvas, 0, 0);
  if (selectedCol >= 0) {
    ctx.fillStyle = "rgba(37,99,235,.13)";
    ctx.fillRect(selectedCol * cellW, 0, cellW, canvas.height);
  }
  if (hoverRow >= 0 && hoverRow < rows) {
    ctx.fillStyle = "rgba(181,71,8,.11)";
    ctx.fillRect(0, hoverRow * cellH, canvas.width, Math.max(0.25, cellH));
  }
  if (hoverCol >= 0 && hoverCol < cols) {
    ctx.fillStyle = "rgba(181,71,8,.09)";
    ctx.fillRect(hoverCol * cellW, 0, cellW, canvas.height);
  }
  if (hoverCol >= 0) {
    ctx.strokeStyle = "#b54708";
    ctx.lineWidth = 1;
    ctx.strokeRect(hoverCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), Math.max(1, canvas.height - 1));
  }
  if (hoverRow >= 0) {
    ctx.strokeStyle = "#b54708";
    ctx.lineWidth = 1;
    ctx.strokeRect(0.5, hoverRow * cellH + 0.5, Math.max(1, canvas.width - 1), Math.max(0.25, cellH - 0.1));
  }
  ctx.strokeStyle = "#111827";
  ctx.lineWidth = Math.max(1, Math.min(2, cellW));
  if (selectedCol >= 0) {
    ctx.strokeRect(selectedCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), Math.max(1, canvas.height - 1));
  }
  if (scrollToSelected && selectedCol >= 0) {
    matrixPane.scrollLeft = Math.max(0, selectedCol * cellW - 80);
    headerScroller.scrollLeft = matrixPane.scrollLeft;
  }
  const cluster = data.clusters[selected];
  const nextDetailsKey = cluster
    ? `${selected}|${referenceOnly}|${mechanismFilterKey()}|${visibleClusters.length}`
    : "none";
  if (nextDetailsKey !== detailsKey && cluster) {
    const counts = countsForCluster(selected);
    const mechanism = processKind(cluster);
    const support = cluster.support ? `<span>support ${escapeHtml(cluster.support)}</span>` : "";
    const phylo = cluster.phylo ? `<span>phylogeny ${escapeHtml(cluster.phylo)}</span>` : "";
    const detail = cluster.mechanismDetail ? `<span>${escapeHtml(cluster.mechanismDetail)}</span>` : "";
    const wording = cluster.wording ? `<div class="process-wording">${escapeHtml(cluster.wording)}</div>` : "";
    info.textContent = `${cluster.present}/${data.genomes.length} genomes (${Math.round(cluster.prevalence * 1000) / 10}%)`;
    details.innerHTML = `<strong>${escapeHtml(cluster.name)}</strong><small>${escapeHtml(cluster.annotation || "No annotation provided")}</small><div class="state-counts"><span class="state-chip">input ${counts.P}</span><span class="state-chip">rescued ${counts.R}</span><span class="state-chip">NSC ${counts.N}</span><span class="state-chip">absent ${counts.A}</span></div><div class="process-note"><span class="process-chip" style="--process-color:${processColor(mechanism)}"><i></i>${escapeHtml(processLabel(mechanism))}</span><span>${escapeHtml(processDescription(mechanism))}</span></div>${referenceInterpretation(cluster)}<div class="process-meta">${support}${phylo}${detail}</div>${wording}`;
    renderProcessLegend();
    detailsKey = nextDetailsKey;
  } else if (nextDetailsKey !== detailsKey) {
    info.textContent = "";
    details.textContent = "No cluster selected.";
    processLegend.innerHTML = "";
    detailsKey = nextDetailsKey;
  }
}
function render(scrollToSelected = false) {
  renderMatrix(scrollToSelected);
  renderColumnHeader();
  renderTree();
  renderLineageLegend();
  renderSyntenyPanel();
}
function scheduleRender(scrollToSelected = false) {
  if (scrollToSelected) {
    render(true);
    return;
  }
  if (renderQueued) return;
  renderQueued = true;
  requestAnimationFrame(() => {
    renderQueued = false;
    render(false);
  });
}
function setHover(row, col, text) {
  const changed = hoverRow !== row || hoverCol !== col;
  hoverRow = row;
  hoverCol = col;
  if (text) status.textContent = text;
  if (changed) scheduleRender(false);
}
function setTreeMode(mode) {
  treeMode = mode;
  phylogramButton.classList.toggle("active", mode === "phylogram");
  cladogramButton.classList.toggle("active", mode === "cladogram");
  treeDirty = true;
  render(false);
}
function resetMatrixCache() {
  baseValid = false;
  hoverCol = -1;
  treeDirty = true;
}
treeRoot = rootByL8(data.tree);
leafOrder = [];
maxDepth = 0;
maxDistance = 0;
layout(treeRoot, 0, 0);
computeGenomeStats();
refreshVisibleClusters();
if (leafOrder.length > 1500) {
  labelToggle.checked = false;
  sizeControl.value = "2";
  rowHeightControl.value = "2";
} else if (leafOrder.length > 800) {
  labelToggle.checked = false;
  sizeControl.value = "3";
  rowHeightControl.value = "3";
} else if (leafOrder.length > 300) {
  labelToggle.checked = false;
  sizeControl.value = "4";
  rowHeightControl.value = "4";
}
populateOptions();
applySyntenyVisibility();
render(true);
search.addEventListener("input", () => { populateOptions(); treeDirty = true; render(true); });
select.addEventListener("change", () => { selected = Number(select.value || 0); treeDirty = true; render(true); });
colorMode.addEventListener("change", () => { treeDirty = true; render(false); });
labelToggle.addEventListener("change", () => { treeDirty = true; render(false); });
phylogramButton.addEventListener("click", () => setTreeMode("phylogram"));
cladogramButton.addEventListener("click", () => setTreeMode("cladogram"));
referenceButton.addEventListener("click", () => {
  if (!data.referenceGenome) return;
  referenceOnly = !referenceOnly;
  refreshVisibleClusters();
  populateOptions();
  resetMatrixCache();
  render(true);
});
syntenyToggle.addEventListener("click", () => {
  syntenyHidden = !syntenyHidden;
  applySyntenyVisibility();
});
processLegend.addEventListener("click", event => {
  const clear = event.target.closest("[data-clear-mechanism]");
  const card = event.target.closest("[data-mechanism]");
  if (!clear && !card) return;
  if (clear) {
    mechanismFilters.clear();
  } else {
    const mechanism = card.getAttribute("data-mechanism") || "";
    if (mechanismFilters.has(mechanism)) {
      mechanismFilters.delete(mechanism);
    } else {
      mechanismFilters.add(mechanism);
    }
  }
  refreshVisibleClusters();
  populateOptions();
  resetMatrixCache();
  render(true);
});
fullscreenButton.addEventListener("click", async () => {
  if (!document.fullscreenElement && explorerSection.requestFullscreen) {
    await explorerSection.requestFullscreen();
  } else if (document.exitFullscreen) {
    await document.exitFullscreen();
  }
});
document.addEventListener("fullscreenchange", () => {
  fullscreenButton.textContent = document.fullscreenElement === explorerSection
    ? "Exit fullscreen"
    : "Fullscreen explorer";
});
sizeControl.addEventListener("input", () => { resetMatrixCache(); render(false); });
rowHeightControl.addEventListener("input", () => { resetMatrixCache(); render(false); });
matrixColorMode.addEventListener("change", () => { resetMatrixCache(); render(false); });
matrixPane.addEventListener("scroll", () => { headerScroller.scrollLeft = matrixPane.scrollLeft; });
headerScroller.addEventListener("scroll", () => { matrixPane.scrollLeft = headerScroller.scrollLeft; });
function adjustRange(control, delta, min, max) {
  const value = Number(control.value || min);
  const next = Math.max(min, Math.min(max, value + delta));
  if (next === value) return false;
  control.value = String(next);
  return true;
}
function wheelZoom(event, forcedAxis = "") {
  const columnZoom = forcedAxis === "column" || event.shiftKey;
  const rowZoom = forcedAxis === "row" || event.ctrlKey || event.metaKey || event.altKey;
  if (!columnZoom && !rowZoom) return;
  event.preventDefault();
  event.stopPropagation();
  const delta = event.deltaY < 0 ? 1 : -1;
  const changed = columnZoom
    ? adjustRange(sizeControl, delta, 0.25, 18)
    : adjustRange(rowHeightControl, delta, 1, 18);
  if (changed) {
    resetMatrixCache();
    render(false);
  }
}
viewerBody.addEventListener("wheel", wheelZoom, {passive:false});
if (treePane) treePane.addEventListener("wheel", event => wheelZoom(event, "row"), {passive:false});
svg.addEventListener("wheel", event => wheelZoom(event, "row"), {passive:false});
headerScroller.addEventListener("wheel", event => wheelZoom(event, "column"), {passive:false});
headerCanvas.addEventListener("wheel", event => wheelZoom(event, "column"), {passive:false});
matrixPane.addEventListener("wheel", wheelZoom, {passive:false});
svg.addEventListener("mousemove", event => {
  const rect = svg.getBoundingClientRect();
  const r = Math.floor((event.clientY - rect.top) / cellH);
  if (r < 0 || r >= leafOrder.length) return;
  const genomeIdx = leafOrder[r];
  const state = stateAt(selected, genomeIdx);
  const cluster = data.clusters[selected];
  const tip = tipStyle(genomeIdx);
  const tipNode = findLeafByGenome(treeRoot, genomeIdx);
  const branchInfo = tipNode
    ? `; root distance ${formatBranchLength(tipNode.distance)}, terminal branch ${formatBranchLength(tipNode.length)}`
    : "";
  const detail = colorMode.value === "cluster" && cluster
    ? `${cluster.name} / ${data.genomes[genomeIdx]}: ${labels[state]}${branchInfo}`
    : `${data.genomes[genomeIdx]}: ${tip.title}${branchInfo}`;
  tooltip.style.display = "none";
  setHover(r, -1, detail);
});
svg.addEventListener("mouseleave", () => {
  setHover(-1, -1, defaultStatus);
});
headerCanvas.addEventListener("mousemove", event => {
  const rect = headerCanvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  if (c < 0 || c >= visibleClusters.length) return;
  const clusterIdx = originalClusterAtVisibleColumn(c);
  const cluster = data.clusters[clusterIdx];
  setHover(-1, c, `${cluster.name}: ${cluster.annotation || "No annotation provided"}`);
});
headerCanvas.addEventListener("mouseleave", () => {
  setHover(-1, -1, defaultStatus);
});
headerCanvas.addEventListener("click", event => {
  const rect = headerCanvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  if (c < 0 || c >= visibleClusters.length) return;
  selected = originalClusterAtVisibleColumn(c);
  search.value = "";
  populateOptions();
  treeDirty = true;
  render(true);
});
canvas.addEventListener("mousemove", event => {
  const rect = canvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  const r = Math.floor((event.clientY - rect.top) / cellH);
  if (c < 0 || r < 0 || c >= visibleClusters.length || r >= leafOrder.length) return;
  const clusterIdx = originalClusterAtVisibleColumn(c);
  const genomeIdx = leafOrder[r];
  const state = stateAt(clusterIdx, genomeIdx);
  const cluster = data.clusters[clusterIdx];
  const mechanism = processKind(cluster);
  setHover(r, c, `${cluster.name} / ${data.genomes[genomeIdx]}: ${labels[state]}; ${processLabel(mechanism)}`);
  const paneRect = matrixPane.getBoundingClientRect();
  tooltip.style.display = "block";
  tooltip.style.left = `${event.clientX - paneRect.left + matrixPane.scrollLeft + 12}px`;
  tooltip.style.top = `${event.clientY - paneRect.top + 12}px`;
  tooltip.innerHTML = `<strong>${escapeHtml(cluster.name)}</strong><br>${escapeHtml(data.genomes[genomeIdx])}: ${escapeHtml(labels[state])}<br><span>${escapeHtml(processLabel(mechanism))}</span><br><span>${escapeHtml(cluster.annotation || "")}</span>`;
});
canvas.addEventListener("mouseleave", () => {
  tooltip.style.display = "none";
  setHover(-1, -1, defaultStatus);
});
canvas.addEventListener("click", event => {
  const rect = canvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  if (c < 0 || c >= visibleClusters.length) return;
  selected = originalClusterAtVisibleColumn(c);
  search.value = "";
  populateOptions();
  treeDirty = true;
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

    let mut html = String::new();
    html.push_str("<!doctype html><html lang='en'><head><meta charset='utf-8'>");
    html.push_str("<meta name='viewport' content='width=device-width,initial-scale=1'>");
    html.push_str("<title>Recore run report</title>");
    html.push_str(
        "<style>
        :root{color-scheme:light;--ink:#111827;--muted:#667085;--subtle:#475467;--line:#d7dde5;--soft-line:#e6eaf0;--bg:#f3f5f7;--panel:#ffffff;--navy:#182230;--blue:#2563eb;--green:#079455;--purple:#7a5af8;--amber:#b54708;--red:#b42318}
        *{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:13px/1.48 -apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif}
        main{max-width:1200px;margin:0 auto;padding:28px}
        header{background:var(--panel);border:1px solid var(--line);border-top:5px solid var(--navy);border-radius:12px;padding:24px;margin-bottom:14px;box-shadow:0 12px 28px rgba(16,24,40,.06)}
        .header-main{display:flex;justify-content:space-between;gap:24px;align-items:flex-start}.brand{display:flex;align-items:center;gap:18px;min-width:0}.logo{width:118px;height:66px;border-radius:8px;overflow:hidden;background:#fff;border:1px solid var(--soft-line);padding:6px;flex:0 0 auto}.logo svg{width:100%;height:100%;display:block}.hero-text{min-width:0}.path{word-break:break-all;margin-top:8px;color:var(--muted);font-size:12px}
        .kicker{text-transform:uppercase;letter-spacing:.13em;font-size:10px;color:var(--muted);font-weight:800}h1{margin:5px 0 7px;font-size:29px;line-height:1.12;letter-spacing:-.01em}h2{font-size:15px;margin:0 0 12px}p{margin:0}.muted{color:var(--muted)}.report-stamp{text-align:right;min-width:180px}.report-stamp strong{display:block;font-size:12px;color:var(--navy)}.report-stamp span{display:block;color:var(--muted);font-size:12px;margin-top:3px}
        .meta{display:flex;gap:8px;flex-wrap:wrap;margin-top:18px}.pill{border:1px solid var(--line);border-radius:999px;padding:5px 9px;background:#f8fafc;color:var(--subtle);font-size:12px}
        .executive{display:grid;grid-template-columns:minmax(0,1.45fr) minmax(320px,.9fr);gap:18px;margin-bottom:14px}.executive p{font-size:14px;color:#344054}.executive strong{color:var(--navy)}.summary-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}.summary-item{border-left:3px solid var(--navy);background:#f8fafc;padding:10px 12px;border-radius:8px}.summary-item span{display:block;color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.06em}.summary-item strong{display:block;margin-top:3px;font-size:16px}
        .grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}.panel,.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:15px;box-shadow:0 8px 18px rgba(16,24,40,.04)}
        .card .label{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.05em}.card .value{font-size:25px;font-weight:760;margin-top:5px;letter-spacing:-.01em}.card .sub{color:var(--muted);font-size:12px;margin-top:3px}
        .wide{grid-column:span 2}.full{grid-column:1/-1}.section-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px;margin-top:12px}
        .bars{display:grid;gap:10px}.bar-label{display:flex;justify-content:space-between;gap:12px;font-size:12px;margin-bottom:4px}.track{height:8px;background:#eef2f6;border-radius:999px;overflow:hidden}.bar{height:100%;background:linear-gradient(90deg,#1d4ed8,#079455);border-radius:999px}
        table{width:100%;border-collapse:collapse;font-size:12px}th{text-align:left;color:var(--muted);font-weight:750;border-bottom:1px solid var(--line);padding:8px 6px;background:#f8fafc}td{border-bottom:1px solid #eaecf0;padding:8px 6px;vertical-align:top}tr:last-child td{border-bottom:0}
        .badge{display:inline-flex;border-radius:999px;padding:2px 8px;font-size:11px;font-weight:750;background:#eef4ff;color:#1849a9;border:1px solid #d1e0ff}.badge.review{background:#fff6ed;color:var(--amber);border-color:#fedf89}.badge.filter{background:#fef3f2;color:var(--red);border-color:#fecdca}.badge.pass{background:#ecfdf3;color:var(--green);border-color:#abefc6}
        .files{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}.file{text-decoration:none;color:var(--ink);border:1px solid var(--line);border-radius:9px;padding:12px;background:#fff;display:block}.file span{font-weight:750;display:block}.file small{color:var(--muted)}
        .explorer-head{display:flex;justify-content:space-between;gap:16px;align-items:flex-start;border-bottom:1px solid var(--soft-line);padding-bottom:12px}.legend{display:flex;gap:10px;flex-wrap:wrap;font-size:12px;color:var(--muted);justify-content:flex-end}.sw{display:inline-block;width:10px;height:10px;border-radius:2px;border:1px solid var(--line);margin-right:5px;vertical-align:-1px}.sw.p{background:#2563eb}.sw.r{background:#079455}.sw.n{background:#7a5af8}.sw.a{background:#fff}
        .explorer-controls,.viewer-controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px;align-items:center;margin-top:14px}.cluster-controls{grid-template-columns:minmax(240px,1fr) minmax(280px,1.2fr) auto}.viewer-controls{position:sticky;top:8px;z-index:6;margin-top:12px;padding:10px;border:1px solid var(--soft-line);border-radius:10px;background:rgba(255,255,255,.97);box-shadow:0 8px 20px rgba(16,24,40,.08)}.explorer input,.explorer select,.explorer button{border:1px solid var(--line);border-radius:8px;padding:8px 10px;background:#fff;color:var(--ink);font:inherit;min-width:0}.explorer button{cursor:pointer;font-weight:700}.explorer button:disabled{cursor:not-allowed;color:#98a2b3;background:#f8fafc}.action-button.active,.segmented-control button.active{background:#eef4ff;border-color:#84adff;color:#1849a9}.segmented-control{display:grid;grid-template-columns:1fr 1fr;gap:4px;border:1px solid var(--line);border-radius:9px;padding:3px;background:#f8fafc}.segmented-control button{border:0;border-radius:6px;background:transparent;padding:6px 8px}.scale-control,.toggle-control{display:flex;align-items:center;gap:8px;color:var(--muted);font-size:12px}.scale-control input{padding:0}.toggle-control input{width:14px;height:14px}#clusterInfo{white-space:nowrap}.cluster-details{margin-top:10px;border:1px solid var(--soft-line);background:#f8fafc;border-radius:9px;padding:10px 12px}.cluster-details strong{font-size:13px}.cluster-details small{display:block;color:var(--muted);margin-top:2px}.state-counts{display:flex;gap:8px;flex-wrap:wrap;margin-top:8px}.state-chip{font-size:11px;border:1px solid var(--line);border-radius:999px;background:#fff;padding:3px 8px}.process-note{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-top:9px;color:#344054;font-size:12px}.process-chip{display:inline-flex;align-items:center;gap:6px;border:1px solid var(--line);border-radius:999px;background:#fff;padding:3px 8px;font-size:11px;font-weight:750;color:#344054}.process-chip i{display:inline-block;width:9px;height:9px;border-radius:2px;background:var(--process-color)}.reference-meta{display:flex;gap:7px;align-items:center;flex-wrap:wrap;margin-top:8px;color:#344054;font-size:11px}.reference-meta strong{font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}.reference-meta span{border:1px solid var(--soft-line);border-radius:999px;background:#fff;padding:2px 7px}.process-meta{display:flex;gap:8px;flex-wrap:wrap;margin-top:6px;color:var(--muted);font-size:11px}.process-meta span{border:1px solid var(--soft-line);border-radius:999px;background:#fff;padding:2px 7px}.process-wording{margin-top:7px;color:#344054;font-size:12px}.process-legend{display:grid;gap:8px;margin-top:10px}.process-selected{display:flex;gap:8px;align-items:center;flex-wrap:wrap;color:#344054;font-size:12px}.process-clear{border-radius:999px!important;padding:4px 9px!important;background:#eef4ff!important;border-color:#84adff!important;color:#1849a9!important;font-size:11px!important}.process-list{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:8px}.process-item{display:grid!important;grid-template-columns:12px minmax(0,1fr);gap:8px;align-items:start;text-align:left;border:1px solid var(--soft-line);border-radius:9px;background:#fff;padding:8px}.process-item:hover{border-color:#84adff;background:#f8fbff}.process-item.active{border-color:#84adff;background:#eef4ff}.process-item.filtered{outline:2px solid #1849a9;outline-offset:0}.process-swatch{width:11px;height:11px;border-radius:3px;margin-top:3px}.process-item strong{display:block;font-size:12px}.process-item small{display:block;margin-top:2px;color:var(--muted);font-size:11px;line-height:1.35;font-weight:400}.lineage-legend{display:none;align-items:center;gap:8px;flex-wrap:wrap;margin-top:8px;border:1px solid var(--soft-line);border-radius:9px;background:#fff;padding:8px 10px}.lineage-title{font-size:11px;font-weight:800;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}.lineage-chip{display:inline-flex;align-items:center;gap:5px;border:1px solid var(--soft-line);border-radius:999px;padding:3px 7px;background:#f8fafc}.lineage-chip i{width:9px;height:9px;border-radius:2px;display:inline-block}.lineage-chip strong{font-size:11px}.lineage-chip small{font-size:11px;color:var(--muted)}
        .viewer-split{display:grid;grid-template-columns:minmax(0,1fr) minmax(390px,1fr);gap:12px;align-items:start;margin-top:14px}.viewer-split>*{min-width:0}.viewer-split.synteny-hidden{grid-template-columns:minmax(0,1fr)}.viewer-split.synteny-hidden .synteny-panel{display:none}.viewer{border:1px solid var(--line);border-radius:10px;background:#fff;overflow:hidden}.viewer-top{display:grid;grid-template-columns:410px minmax(0,1fr);border-bottom:1px solid var(--soft-line);background:#f8fafc}.tree-title,.matrix-title{padding:8px 10px;font-size:11px;font-weight:800;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}.matrix-title{border-left:1px solid var(--soft-line);position:sticky;left:0;background:#f8fafc;z-index:2}.matrix-scroller{overflow-x:auto;overflow-y:hidden}.column-head{display:block;height:24px;image-rendering:auto;cursor:crosshair}.viewer-body{display:grid;grid-template-columns:410px minmax(0,1fr);max-height:680px;overflow-y:auto}.tree-pane{border-right:1px solid var(--soft-line);overflow:hidden;background:#fff}.matrix-pane{position:relative;overflow-x:auto;overflow-y:hidden}.tree-pane svg,.matrix-pane canvas{display:block}.tree-pane svg,.matrix-pane canvas{cursor:crosshair}.matrix-pane canvas{image-rendering:pixelated}.tree-edge{stroke:#98a2b3;stroke-width:1.1;fill:none}.tree-row-guide{stroke:#f97316;stroke-width:.9;stroke-opacity:.32}.tree-scale{font-size:9px;fill:#667085}.tip-label{font-size:10px;fill:#344054}.matrix-tooltip{position:absolute;display:none;pointer-events:none;z-index:4;background:#111827;color:white;border-radius:7px;padding:7px 9px;font-size:12px;box-shadow:0 8px 20px rgba(16,24,40,.22);max-width:280px}.matrix-tooltip span{color:#d0d5dd}.synteny-panel{border:1px solid var(--line);border-radius:10px;background:#fff;min-height:260px;max-height:760px;overflow:auto}.synteny-head{position:sticky;top:0;z-index:3;background:#f8fafc;border-bottom:1px solid var(--soft-line);padding:9px 10px}.synteny-head span{display:block;font-size:11px;font-weight:800;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}.synteny-head small{display:block;margin-top:3px;color:#344054;font-size:12px}.synteny-svg{display:block;min-width:620px}.synteny-link{stroke:#84adff;stroke-width:7;stroke-linecap:round;stroke-opacity:.18;fill:none}.synteny-gene{stroke:#98a2b3;stroke-width:1;rx:3;ry:3}.synteny-target{stroke:#1849a9;stroke-width:1.3}.synteny-label{font-size:9px;fill:#344054}.synteny-target-label{font-size:9px;fill:#fff;font-weight:700}.synteny-row-label{font-size:10px;fill:#475467;font-weight:700}.synteny-locus{font-size:10px;fill:#111827;font-weight:750}.synteny-locus-line{stroke:#475467;stroke-width:.8;stroke-dasharray:2 2}.synteny-empty{padding:14px;color:var(--muted);font-size:12px}.matrix-status{margin-top:8px;font-size:12px}.explorer:fullscreen{overflow:auto;background:var(--bg);padding:18px}.explorer:fullscreen .viewer-split{grid-template-columns:minmax(0,1fr) minmax(430px,1fr)}.explorer:fullscreen .viewer-split.synteny-hidden{grid-template-columns:minmax(0,1fr)}.explorer:fullscreen .viewer-body{max-height:calc(100vh - 310px)}
        footer{color:var(--muted);font-size:12px;margin:18px 4px;border-top:1px solid var(--soft-line);padding-top:12px}
        @media(max-width:1100px){.explorer-controls,.viewer-controls{grid-template-columns:minmax(0,1fr) minmax(0,1fr)}}
        @media(max-width:820px){.viewer-split{grid-template-columns:1fr}}
        @media(max-width:900px){.header-main,.executive{grid-template-columns:1fr;display:block}.report-stamp{text-align:left;margin-top:14px}.grid,.section-grid,.files,.explorer-grid,.explorer-controls,.viewer-controls{grid-template-columns:1fr}.wide{grid-column:auto}main{padding:16px}h1{font-size:24px}.brand{align-items:flex-start}.logo{width:92px;height:54px}.summary-grid{grid-template-columns:1fr}}
        </style>",
    );
    html.push_str("</head><body><main>");

    html.push_str("<header><div class='header-main'><div class='brand'><div class='logo'>");
    html.push_str(REPORT_LOGO);
    write!(
        html,
        "</div><div class='hero-text'><div class='kicker'>Recore audit report</div><h1>Corrected pangenome matrix review</h1><p class='muted'>A reproducible audit of rescued false absences, guarded matrix corrections, and downstream interpretation signals.</p></div></div><div class='report-stamp'><strong>Recore v{}</strong><span>schema recore.outputs.v6</span><span>{} clusters x {} genomes</span></div></div><div class='meta'><span class='pill'>policy {:?}</span><span class='pill'>synteny {}</span><span class='pill'>NSC {}</span><span class='pill'>reference {}</span><span class='pill'>HTML explorer {}</span></div><p class='path'>Output: {}</p></header>",
        html_escape(env!("CARGO_PKG_VERSION")),
        format_count(pa.n_clusters()),
        format_count(pa.n_genomes()),
        stats.matrix_policy,
        if args.synteny { "on" } else { "off" },
        if args.nscluster { "on" } else { "off" },
        html_escape(args.reference_genome.as_deref().unwrap_or("not assessed")),
        if phylo_tree.is_some() { "available" } else { "static" },
        html_escape(&output_dir.display().to_string()),
    )
    .expect("write to String is infallible");

    write!(
        html,
        "<section class='panel executive'><div><div class='kicker'>Executive summary</div><p><strong>{}</strong> binary matrix cells were corrected under the <strong>{:?}</strong> policy from <strong>{}</strong> searchable gaps. Recore retained <strong>{}</strong> accepted evidence rows for audit and kept <strong>{}</strong> candidates out of the binary matrix through synteny/paralog guardrails. Core changed from <strong>{}</strong> to <strong>{}</strong> clusters ({:+}); soft-core changed from <strong>{}</strong> to <strong>{}</strong> ({:+}).</p></div><div class='summary-grid'><div class='summary-item'><span>Correction rate</span><strong>{}</strong></div><div class='summary-item'><span>Evidence retained</span><strong>{}</strong></div><div class='summary-item'><span>Guardrail rejections</span><strong>{}</strong></div><div class='summary-item'><span>Review mode</span><strong>{}</strong></div></div></section>",
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
        format_percent(matrix_corrected_total, stats.searchable_gaps),
        format_count(rescued),
        format_count(rejected_total),
        if phylo_tree.is_some() { "tree-linked" } else { "table-only" },
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
