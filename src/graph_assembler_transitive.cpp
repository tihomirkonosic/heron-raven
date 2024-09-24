
#include "graph_assembler_transitive.h"

namespace raven {

GraphAssemblerTransitive::GraphAssemblerTransitive(raven::Graph &graph)
    : graph_{graph} {}

namespace {

bool is_comparable(double a, double b) {
  double eps = 0.12;
  return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
      (b >= a * (1 - eps) && b <= a * (1 + eps));
};

}

std::uint32_t GraphAssemblerTransitive::RemoveTransitiveEdges() {
  std::vector<Edge *> candidate(graph_.nodes_.size(), nullptr);
  std::unordered_set<std::uint32_t> marked_edges;
  for (const auto &it : graph_.nodes_) {
    if (it == nullptr) {
      continue;
    }
    for (auto jt : it->outedges) {
      candidate[jt->head->id] = jt;
    }
    for (auto jt : it->outedges) {
      for (auto kt : jt->head->outedges) {
        if (candidate[kt->head->id] &&
            is_comparable(jt->length + kt->length, candidate[kt->head->id]->length)) {  // NOLINT
          marked_edges.emplace(candidate[kt->head->id]->id);
          marked_edges.emplace(candidate[kt->head->id]->pair->id);
        }
      }
    }
    for (auto jt : it->outedges) {
      candidate[jt->head->id] = nullptr;
    }
  }

  for (auto i : marked_edges) {  // store for force directed layout
    if (i & 1) {
      auto lhs = graph_.edges_[i]->tail->id & ~1UL;
      auto rhs = graph_.edges_[i]->head->id & ~1UL;
      graph_.nodes_[lhs]->transitive.emplace(rhs);
      graph_.nodes_[rhs]->transitive.emplace(lhs);
    }
  }

  graph_.RemoveEdges(marked_edges);
  return marked_edges.size() / 2;
}

}