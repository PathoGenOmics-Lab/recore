//! Compact encoding of the corrected matrix for the HTML explorer.
//!
//! The explorer used to embed one character per (cluster, genome) cell. That is
//! fine for a toy run and impossible for a real one: 10,000 genomes by 30,000
//! clusters is 300 million characters, so the report never finishes parsing.
//!
//! A cell is one of four states, and `Absent` is overwhelmingly the common one,
//! so instead of a dense grid each cluster stores the genome indices that are
//! *not* absent, split by state. Each set picks whichever of three encodings is
//! smallest for its shape:
//!
//! * `Sparse` - delta-encoded varint indices. Wins for accessory clusters.
//! * `Complement` - the same, over the genomes *missing* from the set. Wins for
//!   core clusters, where a set covering every genome costs a single byte.
//! * `Bitset` - one bit per genome. Caps the pathological middle case, where a
//!   cluster sits near 50% prevalence with no run structure to exploit.
//!
//! Sets are written in genome order, which the caller supplies as matrix order;
//! the browser maps to tree-tip order on load.

const MODE_SPARSE: u8 = 0;
const MODE_COMPLEMENT: u8 = 1;
const MODE_BITSET: u8 = 2;

/// Number of state sets stored per cluster: present, rescued, NSC, and
/// evidence kept out of the binary matrix.
pub(super) const STATE_SETS: usize = 4;

pub(super) fn push_varint(out: &mut Vec<u8>, mut value: u64) {
    loop {
        let mut byte = (value & 0x7f) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80;
        }
        out.push(byte);
        if value == 0 {
            break;
        }
    }
}

fn varint_len(value: u64) -> usize {
    let mut len = 1;
    let mut rest = value >> 7;
    while rest != 0 {
        len += 1;
        rest >>= 7;
    }
    len
}

/// Byte cost of delta-encoding a sorted index list, including its length prefix.
fn sparse_len(indices: &[u32]) -> usize {
    let mut len = varint_len(indices.len() as u64);
    let mut previous = 0u64;
    for &idx in indices {
        let idx = u64::from(idx);
        len += varint_len(idx - previous);
        previous = idx;
    }
    len
}

fn push_sparse(out: &mut Vec<u8>, indices: &[u32]) {
    push_varint(out, indices.len() as u64);
    let mut previous = 0u64;
    for &idx in indices {
        let idx = u64::from(idx);
        push_varint(out, idx - previous);
        previous = idx;
    }
}

/// Encode one set of genome indices, choosing the cheapest of the three modes.
fn encode_set(out: &mut Vec<u8>, indices: &[u32], n_genomes: usize) {
    let complement: Vec<u32> = if indices.is_empty() {
        (0..n_genomes as u32).collect()
    } else {
        let mut present = vec![false; n_genomes];
        for &idx in indices {
            if (idx as usize) < n_genomes {
                present[idx as usize] = true;
            }
        }
        (0..n_genomes as u32)
            .filter(|idx| !present[*idx as usize])
            .collect()
    };

    let bitset_len = n_genomes.div_ceil(8);
    let sparse = sparse_len(indices);
    let complement_cost = sparse_len(&complement);

    if sparse <= complement_cost && sparse <= bitset_len {
        out.push(MODE_SPARSE);
        push_sparse(out, indices);
    } else if complement_cost <= bitset_len {
        out.push(MODE_COMPLEMENT);
        push_sparse(out, &complement);
    } else {
        out.push(MODE_BITSET);
        let mut bits = vec![0u8; bitset_len];
        for &idx in indices {
            if (idx as usize) < n_genomes {
                bits[idx as usize / 8] |= 1 << (idx % 8);
            }
        }
        out.extend_from_slice(&bits);
    }
}

/// Encode a single genome set, such as the carriers of a synteny block.
pub(super) fn encode_one_set(indices: &[u32], n_genomes: usize) -> Vec<u8> {
    let mut out = Vec::new();
    encode_set(&mut out, indices, n_genomes);
    out
}

/// Encode one cluster row: the present, rescued and NSC index sets.
pub(super) fn encode_cluster(sets: &[Vec<u32>; STATE_SETS], n_genomes: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for set in sets {
        encode_set(&mut out, set, n_genomes);
    }
    out
}

const B64: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/// Standard base64 with padding. The payload travels inside a JSON string, so it
/// has to survive the HTML and JSON layers without escaping.
pub(super) fn base64(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let b0 = chunk[0] as u32;
        let b1 = chunk.get(1).copied().unwrap_or(0) as u32;
        let b2 = chunk.get(2).copied().unwrap_or(0) as u32;
        let triple = (b0 << 16) | (b1 << 8) | b2;
        out.push(B64[(triple >> 18) as usize & 0x3f] as char);
        out.push(B64[(triple >> 12) as usize & 0x3f] as char);
        if chunk.len() > 1 {
            out.push(B64[(triple >> 6) as usize & 0x3f] as char);
        } else {
            out.push('=');
        }
        if chunk.len() > 2 {
            out.push(B64[triple as usize & 0x3f] as char);
        } else {
            out.push('=');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cell states in the order the sets are written, mirroring the browser.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum CellState {
        PresentInput,
        Rescued,
        NscMerged,
        KeptOut,
        Absent,
    }

    /// Mirror of the browser decoder, so the round trip is checked in one place.
    fn decode_cluster(bytes: &[u8], n_genomes: usize) -> Vec<CellState> {
        let mut states = vec![CellState::Absent; n_genomes];
        let order = [
            CellState::PresentInput,
            CellState::Rescued,
            CellState::NscMerged,
            CellState::KeptOut,
        ];
        let mut pos = 0usize;
        for state in order {
            let mode = bytes[pos];
            pos += 1;
            match mode {
                MODE_SPARSE | MODE_COMPLEMENT => {
                    let (count, used) = read_varint(&bytes[pos..]);
                    pos += used;
                    let mut previous = 0u64;
                    let mut indices = Vec::with_capacity(count as usize);
                    for _ in 0..count {
                        let (delta, used) = read_varint(&bytes[pos..]);
                        pos += used;
                        previous += delta;
                        indices.push(previous as usize);
                    }
                    if mode == MODE_SPARSE {
                        for idx in indices {
                            states[idx] = state;
                        }
                    } else {
                        let mut absent = vec![false; n_genomes];
                        for idx in indices {
                            absent[idx] = true;
                        }
                        for (idx, missing) in absent.iter().enumerate() {
                            if !missing {
                                states[idx] = state;
                            }
                        }
                    }
                }
                MODE_BITSET => {
                    let len = n_genomes.div_ceil(8);
                    for (idx, slot) in states.iter_mut().enumerate() {
                        if bytes[pos + idx / 8] & (1 << (idx % 8)) != 0 {
                            *slot = state;
                        }
                    }
                    pos += len;
                }
                other => panic!("unknown mode {other}"),
            }
        }
        states
    }

    fn read_varint(bytes: &[u8]) -> (u64, usize) {
        let mut value = 0u64;
        let mut shift = 0;
        let mut used = 0;
        loop {
            let byte = bytes[used];
            value |= u64::from(byte & 0x7f) << shift;
            used += 1;
            if byte & 0x80 == 0 {
                break;
            }
            shift += 7;
        }
        (value, used)
    }

    fn round_trip(present: &[u32], rescued: &[u32], nsc: &[u32], n: usize) -> Vec<CellState> {
        let sets = [present.to_vec(), rescued.to_vec(), nsc.to_vec(), Vec::new()];
        let encoded = encode_cluster(&sets, n);
        decode_cluster(&encoded, n)
    }

    #[test]
    fn empty_cluster_is_all_absent() {
        let states = round_trip(&[], &[], &[], 8);
        assert!(states.iter().all(|s| *s == CellState::Absent));
    }

    #[test]
    fn core_cluster_round_trips_through_complement_mode() {
        let all: Vec<u32> = (0..1000).collect();
        let sets = [all.clone(), Vec::new(), Vec::new(), Vec::new()];
        let encoded = encode_cluster(&sets, 1000);
        // A set covering every genome must not cost a byte per genome.
        assert!(
            encoded.len() < 16,
            "core cluster took {} bytes",
            encoded.len()
        );
        let states = decode_cluster(&encoded, 1000);
        assert!(states.iter().all(|s| *s == CellState::PresentInput));
    }

    #[test]
    fn sparse_cluster_round_trips() {
        let states = round_trip(&[0, 5, 999], &[7], &[123], 1000);
        assert_eq!(states[0], CellState::PresentInput);
        assert_eq!(states[5], CellState::PresentInput);
        assert_eq!(states[999], CellState::PresentInput);
        assert_eq!(states[7], CellState::Rescued);
        assert_eq!(states[123], CellState::NscMerged);
        assert_eq!(states[1], CellState::Absent);
    }

    #[test]
    fn alternating_cluster_falls_back_to_bitset() {
        let half: Vec<u32> = (0..2000).filter(|i| i % 2 == 0).collect();
        let sets = [half.clone(), Vec::new(), Vec::new(), Vec::new()];
        let encoded = encode_cluster(&sets, 2000);
        // 1000 sparse indices would cost ~1000 bytes; the bitset is 250.
        assert!(
            encoded.len() < 400,
            "bitset fallback took {} bytes",
            encoded.len()
        );
        let states = decode_cluster(&encoded, 2000);
        for (idx, state) in states.iter().enumerate() {
            let expected = if idx % 2 == 0 {
                CellState::PresentInput
            } else {
                CellState::Absent
            };
            assert_eq!(*state, expected, "genome {idx}");
        }
    }

    #[test]
    fn later_sets_win_over_earlier_ones() {
        // A genome listed in both present and rescued should read as rescued,
        // matching the browser decoder, which applies the sets in order.
        let states = round_trip(&[4], &[4], &[], 8);
        assert_eq!(states[4], CellState::Rescued);
    }

    #[test]
    fn base64_matches_known_vectors() {
        assert_eq!(base64(b""), "");
        assert_eq!(base64(b"f"), "Zg==");
        assert_eq!(base64(b"fo"), "Zm8=");
        assert_eq!(base64(b"foo"), "Zm9v");
        assert_eq!(base64(b"foob"), "Zm9vYg==");
        assert_eq!(base64(b"fooba"), "Zm9vYmE=");
        assert_eq!(base64(b"foobar"), "Zm9vYmFy");
        assert_eq!(base64(&[0xff, 0xfe, 0xfd]), "//79");
    }

    #[test]
    fn payload_stays_small_for_a_realistic_pangenome() {
        // 10,000 genomes; 3,000 core clusters and 27,000 accessory clusters at
        // low prevalence, which is the shape the explorer has to survive.
        let n = 10_000usize;
        let mut total = 0usize;
        for _ in 0..3_000 {
            let all: Vec<u32> = (0..n as u32).collect();
            total += base64(&encode_cluster(
                &[all, Vec::new(), Vec::new(), Vec::new()],
                n,
            ))
            .len();
        }
        for cluster in 0..27_000u32 {
            let present: Vec<u32> = (0..40)
                .map(|k| (cluster.wrapping_mul(97).wrapping_add(k * 211)) % n as u32)
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect();
            total += base64(&encode_cluster(
                &[present, Vec::new(), Vec::new(), Vec::new()],
                n,
            ))
            .len();
        }
        assert!(
            total < 12_000_000,
            "10k x 30k payload would be {total} bytes of base64"
        );
    }
}
