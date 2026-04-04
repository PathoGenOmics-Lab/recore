//! Core type-safe index wrappers.
//!
//! Using raw `usize` for both cluster and genome indices is error-prone:
//! nothing prevents passing a genome index where a cluster index is expected.
//! These newtypes make such errors a compile-time failure.

/// Type-safe index into the cluster (gene) dimension.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ClusterIdx(pub usize);

/// Type-safe index into the genome (sample) dimension.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct GenomeIdx(pub usize);

impl ClusterIdx {
    #[inline]
    pub fn val(self) -> usize {
        self.0
    }
}

impl GenomeIdx {
    #[inline]
    pub fn val(self) -> usize {
        self.0
    }
}

impl std::fmt::Display for ClusterIdx {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "C{}", self.0)
    }
}

impl std::fmt::Display for GenomeIdx {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "G{}", self.0)
    }
}
