mod cli;
mod types;
mod parser;
mod search;
mod rescue;
mod nscluster;
mod pipeline;

use anyhow::Result;
use clap::Parser;
use env_logger::Env;
use log::info;

use cli::Args;
use pipeline::run;

fn main() -> Result<()> {
    let args = Args::parse();

    let default_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(Env::default().default_filter_or(default_level))
        .format_timestamp_secs()
        .init();
    info!("Recore v{}", env!("CARGO_PKG_VERSION"));

    run(&args)?;
    Ok(())
}
