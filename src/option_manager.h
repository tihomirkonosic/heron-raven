
#ifndef RAVEN_OPTION_MANAGER_H
#define RAVEN_OPTION_MANAGER_H

#include <cstdint>
#include <string>

struct Program_Parameters {
  unsigned split = 0;

  std::uint8_t ploidy = 2;

  std::uint8_t kmer_len = 15;
  std::uint8_t window_len = 7;
  std::uint16_t bandwidth = 500;
  std::uint16_t chain_n = 3;
  std::uint16_t match_n = 75;
  std::uint16_t gap_size = 10000;
  double freq = 0.001;
  bool hpc = false;

  std::string error_corrected_reads = "";
  std::string herro_snps_path = "";
  std::string load_paf = "";

  std::int32_t num_polishing_rounds = 0;

  std::string ul_read_path;

  double disagreement = 0.01;
  std::string gfa_path {"final.gfa"};
  std::string input_gfa_path = "";
  bool skip_contruction = false;
  bool skip_loading_fasta = false;
  bool resume = false;
  bool checkpoints = true;

  std::string root_path{"out/raven"};

  bool load_input_gfa = false;
  bool load_input_paf = false;
  std::string input_paf_path = "";
  std::string sequence_path;

  std::string gfa_after_construction_filename {"after_construction.gfa"};
  std::string gfa_after_transitive_filename {"after_transitive.gfa"};
  std::string gfa_after_bubble_filename {"after_bubble.gfa"};
  std::string gfa_after_force_filename {"after_force.gfa"};
  std::string gfa_post_construction_filename {"post_construction.gfa"};
  std::string gfa_post_cleaning_filename {"post_cleaning.gfa"};
  std::string gfa_unitig_graph_filename {"unitig_graph.gfa"};

  std::string paf_after_snp_filename {"afterSNP.paf"};
  std::string paf_after_contained_filename {"afterContained.paf"};
  std::string paf_after_chimeric_filename {"afterChimeric.paf"};

  std::string cereal_filename {"backup.cereal"};

  std::uint32_t num_threads = 1;

  std::size_t max_overlaps = 16;
  std::uint32_t min_unitig_size = 9999;
  std::uint16_t valid_region_length_threshold = 1260;
  std::uint16_t valid_region_coverage_threshold = 4;
  bool paf = false;
  bool print_gfa_seq = false;

  bool compress_homopolymers = false;
};

int ProcessParameters(int argc, char **argv, Program_Parameters& param);


#endif //RAVEN_OPTION_MANAGER_H
