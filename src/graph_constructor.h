
#ifndef RAVEN_GRAPH_CONSTRUCTOR_H
#define RAVEN_GRAPH_CONSTRUCTOR_H

#include <memory>
#include "biosoup/nucleic_acid.hpp"
#include "pile.hpp"
#include "option_manager.h"

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
        Program_Parameters &param);

    void LoadFromGfa(const std::string &gfa_path);
    void LoadHerroSNPs(const std::string &herro_snps_path, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequencess);
    void LoadOverlaps(const std::string &overlaps_path, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, std::vector<std::vector<extended_overlap>> &extended_overlaps);
    void LoadFromPaf(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, const std::string &paf_path);
  private:
    Graph &graph_;
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_;

    double disagreement_;
  };

} // raven

#endif //RAVEN_GRAPH_CONSTRUCTOR_H
