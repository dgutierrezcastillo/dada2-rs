use clap::Parser;
use dada2_rs::{
    parse_sequence_file, dereplicate, Clustering
};

/// dada2-rs: High-performance amplicon denoising and sample inference in pure Rust.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the FASTQ/FASTA file (can be gzipped)
    #[arg(required = true)]
    input_file: String,

    /// Abundance p-value threshold (omega_a)
    #[arg(short = 'a', long, default_value_t = 1e-40)]
    omega_a: f64,

    /// Verbose output
    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    if args.verbose {
        println!("Parsing sequence file: {}", args.input_file);
    }

    // Phase 1: Parse sequence file
    let records = parse_sequence_file(&args.input_file)?;
    
    if args.verbose {
        println!("Successfully parsed {} records.", records.len());
        println!("Performing sequence dereplication...");
    }

    // Phase 1: Dereplicate sequences
    let derep = dereplicate(&records)?;

    // Output stats
    println!("=== DADA2-rs Dereplication Report ===");
    println!("Total Reads:        {}", derep.total_reads);
    println!("Unique Sequences:   {}", derep.uniques.len());
    
    if derep.uniques.is_empty() {
        println!("No sequences found.");
        return Ok(());
    }

    println!("\n=== Running DADA2 Divisive Clustering Core (Phase 3) ===");
    
    // Generate default transition error probability matrix based on max Phred score 40
    let err_mat = dada2_rs::cluster::generate_default_error_matrix(40);
    
    // Run full EM-like clustering inference
    let clustering = Clustering::run_dada(
        &derep.uniques,
        &err_mat,
        args.omega_a,
        1e-40,
        true, // use qualities
        5,    // Match score
        -4,   // Mismatch score
        -8,   // Gap penalty
        -2,   // Homopolymer gap discount
        16,   // Band size
        false, // Greedy locking
        args.verbose,
    );

    println!("\n=== DADA2-rs Denoising (ASV) Inference Report ===");
    println!("Inferred Amplicon Sequence Variants (ASVs): {}", clustering.clusters.len());
    println!("Total Sample Reads:                            {}", clustering.total_reads);

    println!("\nASV Clusters Summary Table:");
    println!(
        "{:<6} {:<15} {:<12} {:<12} {:<10} {:<12}",
        "ASV_ID", "Center_Unique", "Center_Abund", "Total_Reads", "Unique_Seqs", "Birth_Type"
    );
    for (i, cluster) in clustering.clusters.iter().enumerate() {
        println!(
            "{:<6} Unique_{:<9} {:<12} {:<12} {:<10} {:<12}",
            i + 1,
            cluster.center_index + 1,
            clustering.records[cluster.center_index].abundance,
            cluster.total_reads,
            cluster.raw_indices.len(),
            cluster.birth_type
        );
    }

    if args.verbose && !clustering.clusters.is_empty() {
        println!("\n=== Inferred ASV Representative Sequences ===");
        for (i, cluster) in clustering.clusters.iter().enumerate() {
            let seq_str = String::from_utf8_lossy(&derep.uniques[cluster.center_index].sequence);
            let display_seq = if seq_str.len() > 60 {
                format!("{}...", &seq_str[..57])
            } else {
                seq_str.into_owned()
            };
            println!("ASV {} (Abundance = {}): {}", i + 1, clustering.records[cluster.center_index].abundance, display_seq);
        }
    }

    Ok(())
}
