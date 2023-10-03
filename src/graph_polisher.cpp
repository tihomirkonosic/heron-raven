
#include "graph_polisher.h"
#include "biosoup/timer.hpp"
#include "racon/polisher.hpp"

namespace raven {
  Graph_Polisher::Graph_Polisher(Graph &graph, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
      : stage_(0),
        graph_(graph),
        thread_pool_(thread_pool ?
                     thread_pool :
                     std::make_shared<thread_pool::ThreadPool>(1)) {
  }

  void Graph_Polisher::Polish(
      const std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
      std::uint8_t match,
      std::uint8_t mismatch,
      std::uint8_t gap,
      std::uint32_t cuda_poa_batches,
      bool cuda_banded_alignment,
      std::uint32_t cuda_alignment_batches,
      std::uint32_t num_rounds) {
    if (sequences.empty() || num_rounds == 0) {
      return;
    }

    auto unitigs = graph_.GetUnitigs();
    if (unitigs.empty()) {
      return;
    }

    graph_.piles_.clear();  // save memory

    double avg_q = 0.;
    for (const auto &it: sequences) {
      if (it->block_quality.empty()) {
        continue;
      }
      double q = std::accumulate(
          it->block_quality.begin(),
          it->block_quality.end(),
          0.);
      avg_q += q / it->block_quality.size();
    }
    if (avg_q == 0.) {  // when all values equal to '!'
      for (const auto &it: sequences) {
        it->block_quality.clear();
      }
    } else {
      avg_q /= sequences.size();
    }

    auto polisher = racon::Polisher::Create(
        thread_pool_,
        avg_q,
        0.3,
        500,
        true,
        match, mismatch, gap,
        cuda_poa_batches,
        cuda_banded_alignment,
        cuda_alignment_batches);

    while (stage_ < static_cast<std::int32_t>(num_rounds)) {
      auto polished = polisher->Polish(unitigs, sequences, false);
      unitigs.swap(polished);

      for (const auto &it: unitigs) {  // store unitigs
        const auto &node = graph_.nodes_[std::atoi(&it->name[3])];
        std::size_t tag;
        if ((tag = it->name.rfind(':')) != std::string::npos) {
          if (std::atof(&it->name[tag + 1]) > 0) {
            if (node->is_circular) {  // rotate
              auto s = it->InflateData();
              std::size_t b = 0.42 * s.size();
              s = s.substr(b) + s.substr(0, b);
              it->deflated_data = biosoup::NucleicAcid{"", s}.deflated_data;
            }

            node->is_polished = node->pair->is_polished = true;
            node->sequence.deflated_data = node->pair->sequence.deflated_data = it->deflated_data;  // NOLINT
            node->sequence.inflated_len = node->pair->sequence.inflated_len = it->inflated_len;  // NOLINT
          }
        }
      }

      ++stage_;
      if (graph_.use_checkpoints()) {
        biosoup::Timer timer{};
        timer.Start();
        //Store();
        std::cerr << "[raven::Graph::Polish] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }
  }

} // raven