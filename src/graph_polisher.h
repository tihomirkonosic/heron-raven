
#ifndef RAVEN_GRAPH_POLISHER_H
#define RAVEN_GRAPH_POLISHER_H

#include <vector>
#include <memory>
#include "biosoup/nucleic_acid.hpp"
#include "graph.hpp"

namespace raven {

	class Graph_Polisher {
	public:
		Graph_Polisher() = default;
		Graph_Polisher(Graph &graph,
					   std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

		// Racon wrapper
		void Polish(
				const std::vector<std::unique_ptr<biosoup::NucleicAcid>>& sequences,
				std::uint8_t match,
				std::uint8_t mismatch,
				std::uint8_t gap,
				std::uint32_t cuda_poa_batches,
				bool cuda_banded_alignment,
				std::uint32_t cuda_alignment_batches,
				std::uint32_t num_rounds);

	private:
		int stage_;
		Graph &graph_;
		std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
	};

} // raven

#endif //RAVEN_GRAPH_POLISHER_H
