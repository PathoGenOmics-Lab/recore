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

/// Internal stops in an alignment span, read in the frame the alignment fixes.
///
/// `query_start` is where the alignment starts on the cluster's
/// representative, which begins at its start codon; so the first whole codon
/// of the span, read in the gene's direction, starts `(3 - query_start % 3) %
/// 3` bases in. A rescue with no annotation used to be scanned in all three
/// frames and judged by the cleanest, so a premature stop hid whenever another
/// frame happened to be stop-free, which G/C-ending codons make common, and a
/// pseudogene went into the matrix as present.
pub fn check_internal_stops_in_span(
    contig: &str,
    start: u64,
    end: u64,
    strand: char,
    query_start: u64,
    genome_seqs: &HashMap<String, Vec<u8>>,
) -> Option<usize> {
    let contig_seq = genome_seqs.get(contig)?;
    let start = start as usize;
    let end = (end as usize).min(contig_seq.len());
    if start >= end || end - start < 6 {
        return None;
    }
    let region = &contig_seq[start..end];
    let oriented = if strand == '-' {
        reverse_complement(region)
    } else {
        region.to_vec()
    };
    let offset = ((3 - query_start % 3) % 3) as usize;
    if oriented.len() < offset + 6 {
        return None;
    }
    let n_stops = count_stops_in_frame(&oriented[offset..]);
    (n_stops > 0).then_some(n_stops)
}

/// Translate a coding sequence with the bacterial code (NCBI table 11), from
/// its first base, one amino acid per whole codon, with `*` for stops and `X`
/// for a codon with an ambiguous base. A final stop codon is left off. When
/// `from_start_codon`, the first codon is the gene's start and reads as `M`
/// whatever table 11 start it is (GTG and TTG included), as a CDS does.
pub fn translate_bacterial(seq: &[u8], from_start_codon: bool) -> String {
    const STARTS: [&[u8; 3]; 7] = [b"ATG", b"GTG", b"TTG", b"CTG", b"ATT", b"ATC", b"ATA"];
    let mut protein = String::with_capacity(seq.len() / 3);
    for (i, codon) in seq.chunks_exact(3).enumerate() {
        let upper = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        if i == 0 && from_start_codon && STARTS.contains(&&upper) {
            protein.push('M');
            continue;
        }
        protein.push(codon_amino_acid(upper));
    }
    if protein.ends_with('*') {
        protein.pop();
    }
    protein
}

/// The standard genetic code, which table 11 shares outside start codons.
fn codon_amino_acid(codon: [u8; 3]) -> char {
    let index = |base: u8| match base {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    // Codons in TCAG order: first base slowest, third fastest.
    const TABLE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    match (index(codon[0]), index(codon[1]), index(codon[2])) {
        (Some(a), Some(b), Some(c)) => TABLE[a * 16 + b * 4 + c] as char,
        _ => 'X',
    }
}

/// Reverse complement a DNA sequence.
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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
