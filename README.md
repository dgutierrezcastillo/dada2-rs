# dada2-rs 🧬

> [!WARNING]
> **UNDER ACTIVE DEVELOPMENT (WORK IN PROGRESS)**
> This crate is an early-stage high-performance rewrite of the DADA2 amplicon denoising pipeline in pure Rust. It is **not yet deployed** to crates.io, and it is **not yet production-ready**. API stability is not guaranteed.

`dada2-rs` is a high-performance, parallelized bioinformatic toolkit designed to recreate the core sequence denoising, error-correction, and amplicon processing capabilities of the popular R/C++ package **DADA2** (Divisive Amplicon Denoising Algorithm 2) in pure, safe Rust.

By leveraging Rust's zero-cost abstractions, memory safety guarantees, and robust parallel programming ecosystem (`rayon`), `dada2-rs` aims to deliver substantial speedups and memory footprint reductions for massive-scale amplicon sequencing studies (e.g., PacBio, Illumina, and Oxford Nanopore datasets).

---

## 🚀 Key Features (Phases 1 & 2 Completed)

- [x] **Zero-Copy FASTQ/FASTA Parsing**: Seamless low-overhead ingestion of uncompressed or gzip-compressed files via `needletail`.
- [x] **DADA2-Compliant Dereplication**: Accurately aggregates duplicate sequences into unique reads, tracking abundances and positional consensus quality scores matching R DADA2's `derep-class`.
- [x] **Vectorized 5-Mer Distance Screening**: Computes intersection-based k-mer distance tables for rapid pre-filtering. The nucleotide indexing is optimized at the bit-level to allow automatic compiler SIMD vectorization.
- [x] **Banded Ends-Free Needleman-Wunsch Alignment**: Custom global-local alignment dynamic programming that supports:
  - Banded optimization ($O(L \times W)$ time complexity).
  - Customized **homopolymer gap penalty discounts**, crucial for correcting indel-prone sequencing runs (e.g., PacBio CCS or Oxford Nanopore reads).
- [x] **Multithreaded Execution**: Native parallel processing using `rayon` for scalable, high-throughput analysis.

---

## 📦 Directory and Architecture Overview

The crate is structured around highly optimized modular components located in the `src/` directory:

| Module | File | Description |
| :--- | :--- | :--- |
| **Models** | [`models.rs`](src/models.rs) | Defines basic structs (`DadaRecord`, `UniqueRecord`, and `Derep`) representing dereplicated profiles and reads. |
| **I/O Parser** | [`io.rs`](src/io.rs) | Ingests `.fasta`, `.fastq`, `.fasta.gz`, and `.fastq.gz` sequences with zero-copy stream processing. |
| **Dereplication** | [`derep.rs`](src/derep.rs) | Accumulates exact sequence abundances and maps consensus quality scores using mathematical average mapping. |
| **K-Mer Metrics** | [`kmers.rs`](src/kmers.rs) | Computes count vectors for 5-mers mapped into $4^5 = 1024$ dimensions, generating fast distance matrix screenings. |
| **Sequence Alignment** | [`align.rs`](src/align.rs) | Implements banded, ends-free Needleman-Wunsch dynamic programming with targeted homopolymer run discounts. |

---

## 🛠️ Installation & Building

Since the crate is in active development and **not yet deployed** to crates.io, you can build it locally from the source.

### Prerequisites

You need the Rust toolchain installed. If you do not have it, install it via [rustup.rs](https://rustup.rs):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Compiling from Source

Clone the repository and compile in release mode to enable maximum compiler optimizations:

```bash
git clone https://github.com/dgutierrezcastillo/dada2-rs.git
cd dada2-rs
cargo build --release
```

The compiled binary will be available at `./target/release/dada2-rs`.

---

## 🧪 Verification & Testing

Verify that the algorithms match the DADA2 standards by running our comprehensive suite of unit tests. These test the parser, abundance accumulation, quality consensus averages, k-mer metrics, ends-free bounds, and the homopolymer discount logic:

```bash
cargo test
```

### Running the CLI Showcase

You can run the CLI directly on your FASTX files (compressed or uncompressed). For example:

```bash
cargo run --release -- /path/to/sequences.fasta
```

This will:
1. Parse the input dataset.
2. Accumulate duplicate reads into unique abundances.
3. Print a summary of dereplication abundances (displaying the top unique reads).
4. Perform an ends-free banded pairwise alignment and calculate the 5-mer distance matrix comparing the Rank 1 sequence against other top sequences.

---

## 📄 License

This project is licensed under the **GNU Lesser General Public License version 3.0 (LGPL-3.0)** - see the [LICENSE](LICENSE) file for details. This license choice ensures compatibility with the original upstream R DADA2 project, which is also distributed under the LGPL-3.0.

---

## 📚 Citations & References

If you use these algorithms or ideas in your research, please cite the original DADA2 literature:

1. **DADA2 Method Paper**:
   > Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP (2016). **DADA2: High-resolution sample inference from Illumina amplicon data.** *Nature Methods*, 13(7), 581-583. [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)

2. **Denoising Core Algorithm**:
   > Rosen MJ, Callahan BJ, Fisher DS, Holmes SP (2012). **Denoising PCR-amplified metagenomic data.** *BMC Bioinformatics*, 13(1), 283. [doi:10.1186/1471-2105-13-283](https://doi.org/10.1186/1471-2105-13-283)
