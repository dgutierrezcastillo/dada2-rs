pub mod models;
pub mod io;
pub mod derep;
pub mod kmers;
pub mod align;
pub mod cluster;

// Re-export key functions and types for convenience
pub use models::{DadaRecord, UniqueRecord, Derep};
pub use io::parse_sequence_file;
pub use derep::dereplicate;
pub use kmers::{assign_kmer, kmer_dist};
pub use align::{nwalign_endsfree, nwalign_endsfree_homo, AlignResult};
pub use cluster::{Clustering, Cluster, Comparison};


