
#ifndef RAVEN_GRAPH_CONSTRUCTOR_H
#define RAVEN_GRAPH_CONSTRUCTOR_H

#include <memory>
#include "biosoup/nucleic_acid.hpp"
#include "pile.hpp"

namespace raven {

  class Graph_Constructor {
  public:
    Graph_Constructor() = default;

    Graph_Constructor(Graph &graph,
                      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

    // break chimeric sequences, remove contained sequences and overlaps not
    // spanning bridged repeats at sequence ends
    void Construct(
        std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,  // NOLINT
        double disagreement = 0.1,
        unsigned split = 0,
        std::size_t kMaxNumOverlaps = 16,
				std::uint8_t ploidy = 2,
        std::uint8_t kmer_len = 41,
        std::uint8_t window_len = 41,
        std::uint16_t bandwidth = 500,
        std::uint16_t chain_n = 4,
        std::uint16_t match_n = 100,
        std::uint16_t gap_size = 10000,
        double freq = 0.001,
				bool hpc = false,
        bool paf = false,
        std::uint16_t valid_region_length_threshold = 1260,
        std::uint16_t valid_region_coverage_threshold = 4,
        std::string herro_snps_path = "",
        std::string load_paf = "");

    void LoadFromGfa(const std::string &gfa_path);
    void LoadHerroSNPs(const std::string &herro_snps_path, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequencess);
    void LoadOverlaps(const std::string &overlaps_path, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, std::vector<std::vector<extended_overlap>> &extended_overlaps);
  private:
    Graph &graph_;
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_;

    double disagreement_;
  };

} // raven

#endif //RAVEN_GRAPH_CONSTRUCTOR_H
