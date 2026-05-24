#![forbid(unsafe_code)]

//! Library entry point for Recore.
//!
//! The command-line binary is intentionally thin; all parsing, rescue,
//! reporting, and audit logic lives in this library so it can be tested and
//! reused without going through process-level CLI execution.

pub mod cli;
mod nscluster;
mod parser;
pub mod pipeline;
mod relationships;
mod rescue;
mod search;
mod types;
