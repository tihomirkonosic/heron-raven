#ifndef RAVEN_UNITIG_H
#define RAVEN_UNITIG_H

#include <vector>
#include <memory>
#include "biosoup/nucleic_acid.hpp"
#include "graph.hpp"

namespace raven {

  class Unitig {
  public:
    Unitig() = delete;

    Unitig(std::vector<std::shared_ptr<Node>> nodes, std::vector<std::shared_ptr<Edge>> edges);

    // ignore nodes that are less than epsilon away from any junction node
    std::uint32_t CreateUnitigs(std::uint32_t epsilon = 0);

    std::uint32_t CreateUnitigsAlternate(std::uint32_t epsilon = 0);

    void CreateUnitigGraph();

    void ResolveGraphWithUl(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_reads);

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetUnitigPairs(bool drop_unpolished = false);

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetAssembledData(bool primary = true);

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetUnitigs(bool drop_unpolished = false);

    // draw unitig graph with Bandage
    void PrintUnitigGfa(const std::string& path) const;

    std::vector<std::shared_ptr<Node>> nodes_;
    std::vector<std::shared_ptr<Edge>> edges_;

    std::vector<std::shared_ptr<Node>> unitig_nodes_;
    std::vector<std::shared_ptr<Edge>> unitig_edges_;
  };

}

#endif //RAVEN_UNITIG_H
