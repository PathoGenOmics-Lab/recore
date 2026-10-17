use super::types::RescueEvent;
use crate::types::GenomeIdx;
use std::collections::HashMap;

pub fn resolve_conflicts(rescues: &mut Vec<RescueEvent>) {
    if rescues.is_empty() {
        return;
    }

    // Strongest first; at equal confidence the earlier locus wins, as before.
    // A rescue is kept unless it substantially overlaps one already kept. The
    // pairwise sweep this replaces let a rescue that later lost still remove
    // another whose only conflict was with it.
    let mut order: Vec<usize> = (0..rescues.len()).collect();
    order.sort_by(|&a, &b| {
        let (ra, rb) = (&rescues[a], &rescues[b]);
        rb.confidence
            .total_cmp(&ra.confidence)
            .then(ra.genome_idx.cmp(&rb.genome_idx))
            .then(ra.contig.cmp(&rb.contig))
            .then(ra.start.cmp(&rb.start))
            .then(ra.end.cmp(&rb.end))
            .then(ra.cluster_idx.cmp(&rb.cluster_idx))
    });
    let mut kept_by_locus: HashMap<(GenomeIdx, &str), Vec<usize>> = HashMap::new();
    let mut keep = vec![false; rescues.len()];
    for idx in order {
        let rescue = &rescues[idx];
        let kept = kept_by_locus
            .entry((rescue.genome_idx, rescue.contig.as_str()))
            .or_default();
        if kept
            .iter()
            .any(|&other| is_substantial_overlap(rescue, &rescues[other]))
        {
            continue;
        }
        kept.push(idx);
        keep[idx] = true;
    }
    drop(kept_by_locus);

    let mut idx = 0;
    rescues.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
    rescues.sort_by(|a, b| {
        a.genome_idx
            .cmp(&b.genome_idx)
            .then(a.contig.cmp(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.cluster_idx.cmp(&b.cluster_idx))
    });
}

fn is_substantial_overlap(a: &RescueEvent, b: &RescueEvent) -> bool {
    let overlap = a.end.min(b.end).saturating_sub(a.start.max(b.start));
    if overlap == 0 {
        return false;
    }
    let a_len = a.end.saturating_sub(a.start);
    let b_len = b.end.saturating_sub(b.start);
    let shorter = a_len.min(b_len);
    shorter > 0 && (overlap as f64 / shorter as f64) >= 0.50
}
