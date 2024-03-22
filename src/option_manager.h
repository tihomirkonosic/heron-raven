
#ifndef RAVEN_OPTION_MANAGER_H
#define RAVEN_OPTION_MANAGER_H

#include <cstdint>
#include <string>

struct Program_Parameters {
  unsigned split = 0;

  std::uint8_t ploidy = 2;

  std::uint8_t kmer_len = 31;
  std::uint8_t window_len = 17;
  std::uint16_t bandwidth = 500;
  std::uint16_t chain_n = 4;
  std::uint16_t match_n = 100;
  std::uint16_t gap_size = 10000;
  double freq = 0.001;
  bool hpc = false;

  std::string error_corrected_reads = "";

  std::int32_t num_polishing_rounds = 0;

  std::string ul_read_path;

  double disagreement = 0.1;
  std::string gfa_path = "";
  std::string input_gfa_path = "";
  bool skip_contruction = false;
  bool resume = false;
  bool checkpoints = true;

  bool stdoutput = true;
  std::string output_path = "";

  std::uint32_t num_threads = 1;

  std::size_t max_overlaps = 16;
  std::uint32_t min_unitig_size = 9999;
  std::uint16_t valid_region_size = 4;
  bool paf = false;
  bool print_gfa_seq = false;

  std::string sequence_path;
};

int ProcessParameters(int argc, char **argv, Program_Parameters& param);


#endif //RAVEN_OPTION_MANAGER_H
