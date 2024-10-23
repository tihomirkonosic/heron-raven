
#include "graph_assembler_tips.h"

namespace raven {

GraphAssemblerTips::GraphAssemblerTips(raven::Graph &graph)
    : graph_{graph} {}

std::uint32_t GraphAssemblerTips::RemoveTips() {
  std::uint32_t num_tips = 0;
  std::vector<char> is_visited(graph_.nodes_.size(), 0);

  for (const auto &it: graph_.nodes_) {
    if (it == nullptr || is_visited[it->id] || !it->is_tip()) {
      continue;
    }
    bool is_circular = false;
    std::uint32_t num_sequences = 0;

    auto end = it.get();
    while (!end->is_junction()) {
      num_sequences += end->count;
      is_visited[end->id] = 1;
      is_visited[end->pair->id] = 1;
      if (end->outdegree() == 0 ||
          end->outedges.front()->head->is_junction()) {
        break;
      }
      end = end->outedges.front()->head;
      if (end == it.get()) {
        is_circular = true;
        break;
      }
    }

    if (is_circular || end->outdegree() == 0 || num_sequences > 5) {
      continue;
    }

    std::unordered_set<std::uint32_t> marked_edges;
    for (auto jt: end->outedges) {
      if (jt->head->indegree() > 1) {
        marked_edges.emplace(jt->id);
        marked_edges.emplace(jt->pair->id);
      }
    }
    if (marked_edges.size() / 2 == end->outedges.size()) {  // delete whole
      auto begin = it.get();
      while (begin != end) {
        marked_edges.emplace(begin->outedges.front()->id);
        marked_edges.emplace(begin->outedges.front()->pair->id);
        begin = begin->outedges.front()->head;
      }
      ++num_tips;
    }
    graph_.RemoveEdges(marked_edges, true);
  }

  return num_tips;
}

}