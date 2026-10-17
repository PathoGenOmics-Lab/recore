use crate::search::Hit;
use ahash::HashMap;

/// Scan a nucleotide sequence for internal stop codons (TAA, TAG, TGA) in all three frames.
///
/// The minimum stop count across frames is used because the expected coding
/// frame should be the cleanest one when the sequence is intact.
pub fn count_internal_stops(seq: &[u8]) -> usize {
    if seq.len() < 6 {
        return 0;
    }

    let mut min_stops = usize::MAX;

    for frame in 0..3usize {
        let mut stops = 0usize;
        let mut pos = frame;
        while pos + 3 <= seq.len() {
            let codon = &seq[pos..pos + 3];
            if pos + 3 < seq.len() {
                let upper = [
                    codon[0].to_ascii_uppercase(),
                    codon[1].to_ascii_uppercase(),
                    codon[2].to_ascii_uppercase(),
                ];
                if matches!(
                    upper,
                    [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
                ) {
                    stops += 1;
                }
            }
            pos += 3;
        }
        min_stops = min_stops.min(stops);
    }

    min_stops
}

/// Count internal stop codons in a single reading frame starting at the first
/// base of `seq`, ignoring the final codon (the gene's own stop).
pub fn count_stops_in_frame(seq: &[u8]) -> usize {
    let mut stops = 0usize;
    let mut pos = 0usize;
    while pos + 3 < seq.len() {
        let codon = [
            seq[pos].to_ascii_uppercase(),
            seq[pos + 1].to_ascii_uppercase(),
            seq[pos + 2].to_ascii_uppercase(),
        ];
        if matches!(
            codon,
            [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
        ) {
            stops += 1;
        }
        pos += 3;
    }
    stops
}

/// Internal stops in an annotated CDS, read in its own frame from its own start.
///
/// For a rescue that sits on an annotated CDS this is the right question, and
/// the alignment span is the wrong one. The minimum over three frames of the
/// span let a premature stop hide whenever another frame happened to be
/// stop-free, which GC-rich third positions make common. And the span follows
/// the cluster's representative: when that allele reads through its stop, the
/// span runs past this genome's own, canonical stop, which then counted as
/// internal and an intact gene was called a pseudogene.
pub fn check_internal_stops_in_cds(
    contig: &str,
    start: u64,
    end: u64,
    strand: char,
    genome_seqs: &HashMap<String, Vec<u8>>,
) -> Option<usize> {
    let contig_seq = genome_seqs.get(contig)?;
    let start = start as usize;
    let end = (end as usize).min(contig_seq.len());
    if start >= end || end - start < 6 {
        return None;
    }
    let region = &contig_seq[start..end];
    let n_stops = if strand == '-' {
        count_stops_in_frame(&reverse_complement(region))
    } else {
        count_stops_in_frame(region)
    };
    (n_stops > 0).then_some(n_stops)
}

/// Reverse complement a DNA sequence.
pub(super) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Check if a rescued hit is a pseudogene by scanning the genome sequence for internal stops.
///
/// For hits on the minus strand, the sequence is reverse-complemented before
/// scanning. Returns `Some(n_stops)` if stops are found, otherwise `None`.
pub fn check_internal_stops_in_hit(
    hit: &Hit,
    genome_seqs: Option<&HashMap<String, Vec<u8>>>,
) -> Option<usize> {
    let seqs = genome_seqs?;
    let contig_seq = seqs.get(&hit.target)?;

    let start = hit.target_start as usize;
    let end = (hit.target_end as usize).min(contig_seq.len());
    if start >= end || end - start < 6 {
        return None;
    }

    let region = &contig_seq[start..end];
    let n_stops = if hit.strand == '-' {
        let rc = reverse_complement(region);
        count_internal_stops(&rc)
    } else {
        count_internal_stops(region)
    };

    if n_stops > 0 {
        Some(n_stops)
    } else {
        None
    }
}
