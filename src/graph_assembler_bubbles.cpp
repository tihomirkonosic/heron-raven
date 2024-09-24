
#include "graph_assembler_bubbles.h"
#include "edlib.h"

namespace {

std::unordered_set<std::uint32_t> annotation_extract(
    raven::Graph &graph,
    std::uint32_t i,
    std::uint32_t begin,
    std::uint32_t end,
    std::uint32_t len,
    bool strand) {

  std::unordered_set<std::uint32_t> dst;
  if (graph.annotations_[i].empty()) {
    return dst;
  }

  for (const auto &it : graph.annotations_[i]) {
    if (begin <= it && it <= end) {
      dst.emplace(strand ? it : len - 1 - it);
    }
  }

  return dst;
}

std::vector<raven::Node *> path_extract(raven::Node *begin, raven::Node *end, std::vector<raven::Node *> predecessor) {
  std::vector<raven::Node *> dst;

  while (end != begin) {
    dst.emplace_back(end);
    end = predecessor[end->id];
  }

  dst.emplace_back(begin);
  std::reverse(dst.begin(), dst.end());

  return dst;
}

bool path_is_simple(const std::vector<raven::Node *> &path) {
  if (path.empty()) {
    return false;
  }

  for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
    if (path[i]->is_junction()) {
      return false;  // complex
    }
  }

  return true;  // without branches
}

std::string path_sequence(const std::vector<raven::Node *> &path) {
  std::string data{};
  for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
    for (auto it: path[i]->outedges) {
      if (it->head == path[i + 1]) {
        data += it->Label();
        break;
      }
    }
  }

  data += path.back()->sequence.InflateData();
  return data;
}

std::unordered_set<std::uint32_t> path_annotation(raven::Graph& graph, const std::vector<raven::Node *> &path) {
  std::unordered_set<std::uint32_t> dst;
  std::uint32_t offset = 0;

  for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
    for (auto it: path[i]->outedges) {
      if (it->head == path[i + 1]) {
        const auto &annotations = annotation_extract(
            graph,
            it->tail->sequence.id,
            0,
            it->length,
            it->tail->sequence.inflated_len,
            !it->tail->is_rc());
        for (const auto &jt: annotations) {
          dst.emplace(offset + jt);
        }
        offset += it->length;
      }
    }
  }

  const auto &annotations = annotation_extract(
      graph,
      path.back()->sequence.id,
      0,
      path.back()->sequence.inflated_len,
      path.back()->sequence.inflated_len,
      !path.back()->is_rc());

  for (const auto &jt: annotations) {
    dst.emplace(offset + jt);
  }

  return dst;
}

std::unordered_set<std::uint32_t> bubble_pop_int(raven::Graph& graph, const std::vector<raven::Node *> &lhs, const std::vector<raven::Node *> &rhs) {
  if (lhs.empty() || rhs.empty()) {
    return std::unordered_set<std::uint32_t>{};
  }

  // check BFS result
  std::unordered_set<raven::Node *> bubble;
  bubble.insert(lhs.begin(), lhs.end());
  bubble.insert(rhs.begin(), rhs.end());
  if (lhs.size() + rhs.size() - 2 != bubble.size()) {
    return std::unordered_set<std::uint32_t>{};
  }

  for (const auto &it: lhs) {
    if (bubble.count(it->pair) != 0) {
      return std::unordered_set<std::uint32_t>{};
    }
  }

  if (!path_is_simple(lhs) || !path_is_simple(rhs)) {  // complex path(s)
    // check poppability
    if (graph.FindRemovableEdges(lhs).empty() && graph.FindRemovableEdges(rhs).empty()) {
      return std::unordered_set<std::uint32_t>{};
    }

    // check similarity
    auto l = path_sequence(lhs);
    auto r = path_sequence(rhs);
    if (std::min(l.size(), r.size()) < std::max(l.size(), r.size()) * 0.8) {
      return std::unordered_set<std::uint32_t>{};
    }

    EdlibAlignResult result = edlibAlign(
        l.c_str(), l.size(),
        r.c_str(), r.size(),
        edlibDefaultAlignConfig());
    double score = 0;
    if (result.status == EDLIB_STATUS_OK) {
      score = 1 - result.editDistance /
          static_cast<double>(std::max(l.size(), r.size()));
      edlibFreeAlignResult(result);
    }
    if (score < 0.8) {
      return std::unordered_set<std::uint32_t>{};
    }
    if (!graph.annotations_.empty()) {
      std::unordered_set<std::uint32_t> marked_edges;
      if (!path_is_simple(lhs)) {
        marked_edges = graph.FindRemovableEdges(lhs);
        if (marked_edges.size() > 2) {
          marked_edges.clear();
        }
      }
      if (marked_edges.empty() && !path_is_simple(rhs)) {
        marked_edges = graph.FindRemovableEdges(rhs);
        if (marked_edges.size() > 2) {
          marked_edges.clear();
        }
      }
      return marked_edges;
    }
  }

  if (!graph.annotations_.empty()) {
    auto la = path_annotation(graph, lhs);
    auto ra = path_annotation(graph, rhs);

    if (!la.empty() && !ra.empty()) {
      auto l = path_sequence(lhs);
      auto r = path_sequence(rhs);

      EdlibAlignResult result = edlibAlign(
          l.c_str(), l.size(),
          r.c_str(), r.size(),
          edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

      if (result.status == EDLIB_STATUS_OK) {
        std::uint32_t lhs_pos = 0;
        std::uint32_t rhs_pos = 0;

        std::uint32_t mismatches = 0;
        std::uint32_t snps = 0;

        for (int a = 0; a < result.alignmentLength; ++a) {
          if (la.find(lhs_pos) != la.end() ||
              ra.find(rhs_pos) != ra.end()) {
            ++snps;
            if (result.alignment[a] == 3) {
              ++mismatches;
            }
          }
          switch (result.alignment[a]) {
            case 0:
            case 3: {
              ++lhs_pos;
              ++rhs_pos;
              break;
            }
            case 1: {
              ++lhs_pos;
              break;
            }
            case 2: {
              ++rhs_pos;
              break;
            }
            default:
              break;
          }
        }

        edlibFreeAlignResult(result);

        if (mismatches / static_cast<double>(snps) > 0.1) {  // disagreement_) {
          return std::unordered_set<std::uint32_t>{};
        }
      }
    }
  }

  std::uint32_t lhs_count = 0;
  for (auto jt: lhs) {
    lhs_count += jt->count;
  }

  std::uint32_t rhs_count = 0;
  for (auto jt: rhs) {
    rhs_count += jt->count;
  }

  auto marked_edges = graph.FindRemovableEdges(lhs_count > rhs_count ? rhs : lhs);
  if (marked_edges.empty()) {
    marked_edges = graph.FindRemovableEdges(lhs_count > rhs_count ? lhs : rhs);
  }

  return marked_edges;
}

std::vector<raven::MarkedEdge> bubble_pop_snp(raven::Graph& graph, const std::vector<raven::Node *> &lhs, const std::vector<raven::Node *> &rhs) {
  if (lhs.empty() || rhs.empty()) {
    return std::vector<raven::MarkedEdge>{};
  }

  // check BFS result
  std::unordered_set<raven::Node *> bubble;
  bubble.insert(lhs.begin(), lhs.end());
  bubble.insert(rhs.begin(), rhs.end());
  if (lhs.size() + rhs.size() - 2 != bubble.size()) {
    return std::vector<raven::MarkedEdge>{};
  }
  for (const auto &it: lhs) {
    if (bubble.count(it->pair) != 0) {
      return std::vector<raven::MarkedEdge>{};
    }
  }

  if (!path_is_simple(lhs) || !path_is_simple(rhs)) {  // complex path(s)
    // check poppability
    if (graph.FindRemovableEdges(lhs).empty() && graph.FindRemovableEdges(rhs).empty()) {
      return std::vector<raven::MarkedEdge>{};
    }

    // check similarity
    auto l = path_sequence(lhs);
    auto r = path_sequence(rhs);
    if (std::min(l.size(), r.size()) < std::max(l.size(), r.size()) * 0.8) {
      return std::vector<raven::MarkedEdge>{};
    }

    EdlibAlignResult result = edlibAlign(
        l.c_str(), l.size(),
        r.c_str(), r.size(),
        edlibDefaultAlignConfig());
    double score = 0;
    if (result.status == EDLIB_STATUS_OK) {
      score = 1 - result.editDistance /
          static_cast<double>(std::max(l.size(), r.size()));
      edlibFreeAlignResult(result);
    }
    if (score < 0.8) {
      return std::vector<raven::MarkedEdge>{};
    }
    if (!graph.annotations_.empty()) {
      std::unordered_set<std::uint32_t> marked_edges;
      if (!path_is_simple(lhs)) {
        marked_edges = graph.FindRemovableEdges(lhs);
        if (marked_edges.size() > 2) {
          marked_edges.clear();
        }
      }
      if (marked_edges.empty() && !path_is_simple(rhs)) {
        marked_edges = graph.FindRemovableEdges(rhs);
        if (marked_edges.size() > 2) {
          marked_edges.clear();
        }
      }
      std::vector<raven::MarkedEdge> marked_edge_list;
      for (const auto &it: marked_edges)
        marked_edge_list.emplace_back(it);
      return marked_edge_list;
    }
  }

  if (!graph.annotations_.empty()) {
    auto la = path_annotation(graph, lhs);
    auto ra = path_annotation(graph, rhs);

    if (!la.empty() && !ra.empty()) {
      auto l = path_sequence(lhs);
      auto r = path_sequence(rhs);

      EdlibAlignResult result = edlibAlign(
          l.c_str(), l.size(),
          r.c_str(), r.size(),
          edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

      if (result.status == EDLIB_STATUS_OK) {
        std::uint32_t lhs_pos = 0;
        std::uint32_t rhs_pos = 0;

        std::uint32_t mismatches = 0;
        std::uint32_t snps = 0;

        for (int a = 0; a < result.alignmentLength; ++a) {
          if (la.find(lhs_pos) != la.end() ||
              ra.find(rhs_pos) != ra.end()) {
            ++snps;
            if (result.alignment[a] == 3) {
              ++mismatches;
            }
          }
          switch (result.alignment[a]) {
            case 0:
            case 3: {
              ++lhs_pos;
              ++rhs_pos;
              break;
            }
            case 1: {
              ++lhs_pos;
              break;
            }
            case 2: {
              ++rhs_pos;
              break;
            }
            default:
              break;
          }
        }

        edlibFreeAlignResult(result);

        if (mismatches / static_cast<double>(snps) > 0.1) {  // disagreement_) {
          std::vector<raven::MarkedEdge> marked_edge_list;

          auto marked_edges = graph.FindRemovableEdges(lhs);
          for (const auto &it: marked_edges)
            marked_edge_list.emplace_back(it, 1);

          marked_edges = graph.FindRemovableEdges(rhs);
          for (const auto &it: marked_edges)
            marked_edge_list.emplace_back(it, 2);

          return marked_edge_list;
        }
      }
    }
  }

  std::uint32_t lhs_count = 0;
  for (auto jt: lhs) {
    lhs_count += jt->count;
  }
  std::uint32_t rhs_count = 0;
  for (auto jt: rhs) {
    rhs_count += jt->count;
  }
  auto marked_edges = graph.FindRemovableEdges(lhs_count > rhs_count ? rhs : lhs);
  if (marked_edges.empty()) {
    marked_edges = graph.FindRemovableEdges(lhs_count > rhs_count ? lhs : rhs);
  }
  std::vector<raven::MarkedEdge> marked_edge_list;
  for (const auto &it: marked_edges)
    marked_edge_list.emplace_back(it);
  return marked_edge_list;
}

}

namespace raven {

GraphAssemblerBubbles::GraphAssemblerBubbles(raven::Graph &graph)
    : graph_{graph} {}

std::uint32_t GraphAssemblerBubbles::RemoveBubbles() {
  std::vector<std::uint32_t> distance(graph_.nodes_.size(), 0);
  std::vector<Node *> predecessor(graph_.nodes_.size(), nullptr);

  std::uint32_t num_bubbles = 0;
  for (const auto &it: graph_.nodes_) {
    if (it == nullptr || it->outdegree() < 2) {
      continue;
    }

    // BFS
    Node *begin = it.get();
    Node *end = nullptr;
    Node *other_end = nullptr;
    std::deque<Node *> que{begin};
    std::vector<Node *> visited{1, begin};
    while (!que.empty() && !end) {
      auto jt = que.front();
      que.pop_front();

      for (auto kt: jt->outedges) {
        if (kt->head == begin) {  // cycle
          continue;
        }
        if (distance[jt->id] + kt->length > 5000000) {  // out of reach
          continue;
        }
        if (kt->head->id > distance.size())
          continue;
        distance[kt->head->id] = distance[jt->id] + kt->length;
        visited.emplace_back(kt->head);
        que.emplace_back(kt->head);

        if (predecessor[kt->head->id]) {  // found bubble
          end = kt->head;
          other_end = jt;
          break;
        }

        predecessor[kt->head->id] = jt;
      }
    }

    std::unordered_set<std::uint32_t> marked_edges;
    if (end) {
      auto lhs = path_extract(begin, end, predecessor);
      auto rhs = path_extract(begin, other_end, predecessor);
      rhs.emplace_back(end);
      marked_edges = bubble_pop_int(graph_, lhs, rhs);
    }

    for (auto jt: visited) {
      distance[jt->id] = 0;
      predecessor[jt->id] = nullptr;
    }

    graph_.RemoveEdges(marked_edges, true);
    num_bubbles += 1 - marked_edges.empty();
  }

  return num_bubbles;
}

std::uint32_t GraphAssemblerBubbles::RemoveSnpBubbles() {
  std::vector<std::uint32_t> distance(graph_.nodes_.size(), 0);
  std::vector<Node *> predecessor(graph_.nodes_.size(), nullptr);

  std::uint32_t num_bubbles = 0;
  for (const auto &it: graph_.nodes_) {
    if (it == nullptr || it->outdegree() < 2) {
      continue;
    }

    // BFS
    Node *begin = it.get();
    Node *end = nullptr;
    Node *other_end = nullptr;
    std::deque<Node *> que{begin};
    std::vector<Node *> visited{1, begin};
    while (!que.empty() && !end) {
      auto jt = que.front();
      que.pop_front();

      for (auto kt: jt->outedges) {
        if (kt->head == begin) {  // cycle
          continue;
        }
        if (distance[jt->id] + kt->length > 5000000) {  // out of reach
          continue;
        }
        if (kt->head->id > distance.size())
          continue;
        distance[kt->head->id] = distance[jt->id] + kt->length;
        visited.emplace_back(kt->head);
        que.emplace_back(kt->head);

        if (predecessor[kt->head->id]) {  // found bubble
          end = kt->head;
          other_end = jt;
          break;
        }

        predecessor[kt->head->id] = jt;
      }
    }

    std::vector<MarkedEdge> marked_edges;
    if (end) {
      auto lhs = path_extract(begin, end, predecessor);
      auto rhs = path_extract(begin, other_end, predecessor);
      rhs.emplace_back(end);
      marked_edges = bubble_pop_snp(graph_, lhs, rhs);
    }

    for (auto jt: visited) {
      distance[jt->id] = 0;
      predecessor[jt->id] = nullptr;
    }

    graph_.RemoveDiploidEdges(marked_edges, true);
    num_bubbles += 1 - marked_edges.empty();
    std::cerr << "[raven::Graph::Assemble] RemoveBubbles: num_bubbles=" << num_bubbles << std::endl;
  }

  return num_bubbles;
}

}