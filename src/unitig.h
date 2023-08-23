//
// Created by Tiho on 9.8.2023..
//

#ifndef RAVEN_UNITIG_H
#define RAVEN_UNITIG_H

#include <vector>
#include <memory>
#include "biosoup/nucleic_acid.hpp"
#include "graph.hpp"

namespace raven {

    class Unitig {
    public:
        // ignore nodes that are less than epsilon away from any junction node
        std::uint32_t CreateUnitigs(std::uint32_t epsilon = 0);
        std::uint32_t CreateUnitigsAlternate(std::uint32_t epsilon = 0);

        void CreateUnitigGraph();
        void ResolveGraphWithUl(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_reads);

        std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetUnitigPairs(bool drop_unpolished = false);

        std::vector<std::unique_ptr<biosoup::NucleicAcid>> GetAssembledData(bool primary = true);
        std::vector <std::unique_ptr<biosoup::NucleicAcid>> GetUnitigs(bool drop_unpolished = false);

        std::vector<std::shared_ptr<Node>> nodes_;
        std::vector<std::shared_ptr<Edge>> edges_;

        std::vector<std::shared_ptr<Node>> unitig_nodes_;
        std::vector<std::shared_ptr<Edge>> unitig_edges_;

        std::vector<std::shared_ptr<Node>> nodes_alternate_;
        std::vector<std::shared_ptr<Edge>> edges_alternate_;
    };

}

#endif //RAVEN_UNITIG_H
