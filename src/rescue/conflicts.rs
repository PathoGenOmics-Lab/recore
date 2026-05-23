use super::types::RescueEvent;

pub fn resolve_conflicts(rescues: &mut Vec<RescueEvent>) {
    if rescues.is_empty() {
        return;
    }

    rescues.sort_by(|a, b| {
        a.genome_idx
            .cmp(&b.genome_idx)
            .then(a.contig.cmp(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.cluster_idx.cmp(&b.cluster_idx))
    });

    let mut to_remove: Vec<bool> = vec![false; rescues.len()];

    for i in 0..rescues.len() {
        if to_remove[i] {
            continue;
        }
        for j in (i + 1)..rescues.len() {
            if to_remove[j] {
                continue;
            }

            let ri = &rescues[i];
            let rj = &rescues[j];

            if ri.genome_idx != rj.genome_idx {
                break;
            }
            if ri.contig != rj.contig {
                break;
            }
            if rj.start >= ri.end {
                break;
            }

            if !is_substantial_overlap(ri, rj) {
                continue;
            }

            if ri.confidence >= rj.confidence {
                to_remove[j] = true;
            } else {
                to_remove[i] = true;
                break;
            }
        }
    }

    let mut idx = 0;
    rescues.retain(|_| {
        let keep = !to_remove[idx];
        idx += 1;
        keep
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
