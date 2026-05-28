use clap::Parser;
use dada2_rs::{
    parse_sequence_file, dereplicate, assign_kmer, kmer_dist, nwalign_endsfree
};

/// dada2-rs: High-performance amplicon denoising and sample inference in pure Rust.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the FASTQ/FASTA file (can be gzipped)
    #[arg(required = true)]
    input_file: String,

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
    
    if !derep.uniques.is_empty() {
        println!("\nTop 5 Unique Sequences by Abundance:");
        println!("{:<4} {:<10} {:<7} {:<50}", "Rank", "Abundance", "Length", "Sequence (Prefix)");
        for (i, unique) in derep.uniques.iter().take(5).enumerate() {
            let seq_str = String::from_utf8_lossy(&unique.sequence);
            let display_seq = if seq_str.len() > 50 {
                format!("{}...", &seq_str[..47])
            } else {
                seq_str.into_owned()
            };
            println!(
                "{:<4} {:<10} {:<7} {:<50}",
                i + 1,
                unique.abundance,
                unique.sequence.len(),
                display_seq
            );
        }

        // Phase 2 Showcase: Sequence Comparison and Alignment
        if derep.uniques.len() >= 2 {
            println!("\n=== DADA2-rs Sequence Comparison (Phase 2 Showcase) ===");
            println!("Comparing top unique sequence (Rank 1) against other top sequences:");
            println!(
                "{:<8} {:<10} {:<15} {:<10} {:<30}",
                "Target", "Abundance", "K-mer Distance", "NW Score", "Alignment (Query vs Ref)"
            );

            let rank1_seq = &derep.uniques[0].sequence;
            let rank1_kvec = assign_kmer(rank1_seq, 5);

            for (i, unique) in derep.uniques.iter().skip(1).take(4).enumerate() {
                let target_seq = &unique.sequence;
                let target_kvec = assign_kmer(target_seq, 5);
                
                // Calculate 5-mer distance
                let dist = kmer_dist(&rank1_kvec, rank1_seq.len(), &target_kvec, target_seq.len(), 5);
                
                // Calculate ends-free banded alignment (band = 16)
                let align_res = nwalign_endsfree(
                    rank1_seq,
                    target_seq,
                    5,    // Match
                    -4,   // Mismatch
                    -8,   // Gap penalty
                    16,   // Band
                );

                // Create a brief alignment preview
                let q_align_str = if align_res.query_align.len() > 15 {
                    format!("{}...", &align_res.query_align[..12])
                } else {
                    align_res.query_align.clone()
                };
                let r_align_str = if align_res.ref_align.len() > 15 {
                    format!("{}...", &align_res.ref_align[..12])
                } else {
                    align_res.ref_align.clone()
                };

                println!(
                    "Rank {:<3} {:<10} {:<15.4} {:<10} {:<15} vs {:<15}",
                    i + 2,
                    unique.abundance,
                    dist,
                    align_res.score,
                    q_align_str,
                    r_align_str
                );
            }
        }
    } else {
        println!("No sequences found.");
    }

    Ok(())
}

