
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
                std::uint8_t kmer_len = 15,
                std::uint8_t window_len = 5,
                double freq = 0.001,
                std::uint16_t valid_region_treshold = 4);

    private:
        Graph &graph_;
        std::shared_ptr<thread_pool::ThreadPool> thread_pool_;

        double disagreement_;
    };

} // raven

#endif //RAVEN_GRAPH_CONSTRUCTOR_H
