use std::env;

use upl_finder_rust::{find_upl_matches, load_upl_probes};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <probes.(tsv|json|txt)> <sequence>", args[0]);
        std::process::exit(2);
    }

    let probes_path = &args[1];
    let sequence = &args[2];

    let probes = load_upl_probes(probes_path)?;
    let matches = find_upl_matches(sequence, &probes);
    let out = serde_json::to_string_pretty(&matches)?;
    println!("{out}");
    Ok(())
}
