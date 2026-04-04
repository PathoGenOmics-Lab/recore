mod cli;
mod parser;
mod search;
mod rescue;
mod pipeline;

use anyhow::Result;
use clap::Parser;
use env_logger::Env;
use log::info;

use cli::Args;
use pipeline::run;

fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info"))
        .format_timestamp_secs()
        .init();

    let args = Args::parse();
    info!("Recore v{}", env!("CARGO_PKG_VERSION"));

    run(&args)?;
    Ok(())
}
