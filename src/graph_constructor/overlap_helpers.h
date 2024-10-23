
#ifndef OVERLAP_HELPERS_H
#define OVERLAP_HELPERS_H

#include <string>

inline std::string edlib_alignment_reverse(const std::string &s) {
  std::string rs;
  if (s.empty()) {
    return rs;
  }
  for (char c : s) {
    switch (c) {
      case '2':rs += '\001';
        break;
      case '1':rs += '\002';
        break;
      default:rs += c;
        break;
    }
  }
  return rs;
}

inline std::string cigar_alignment_reverse(const std::string &s) {
  std::string rs;
  if (s.empty()) {
    return rs;
  }
  for (char c : s) {
    switch (c) {
      case 'I':rs += 'D';
        break;
      case 'D':rs += 'I';
        break;
      default:rs += c;
        break;
    }
  }
  return rs;
}

inline std::string cigar_to_edlib_alignment(const std::string &s) {
  std::string rs = "";
  std::uint64_t pos = 0;
  std::uint64_t start_pos = 0;
  std::uint64_t total_num = 0;

  if (s.empty()) {
    return rs;
  }
  for (int i = 0; i < s.length(); i++) {
    if (std::isdigit(s[i])) {
      if (pos == 0) {
        start_pos = i;
      }
      ++pos;
    } else {
      total_num = 0;
      for (int j = start_pos; j < start_pos + pos; j++) {
        total_num += (s[j] - '0') * std::pow(10, (start_pos + pos) - j - 1);
      }
      pos = 0;
      switch (s[i]) {
        case '=':
          for (int j = 0; j < total_num; j++) {
            rs += '\000';
          };
          break;
        case 'X':
          for (int j = 0; j < total_num; j++) {
            rs += '\003';
          };
          break;
        case 'I':
          for (int j = 0; j < total_num; j++) {
            rs += '\001';
          };
          break;
        case 'D':
          for (int j = 0; j < total_num; j++) {
            rs += '\002';
          };
          break;
        default:
          //rs += '\000';
          break;
      }
    }
  }
  return rs;
}

inline biosoup::Overlap overlap_reverse(const biosoup::Overlap &o) {
  return biosoup::Overlap(
      o.rhs_id, o.rhs_begin, o.rhs_end,
      o.lhs_id, o.lhs_begin, o.lhs_end,
      o.score, edlib_alignment_reverse(o.alignment),
      o.strand);
}

inline biosoup::Overlap cigar_overlap_reverse(const biosoup::Overlap &o) {
  return biosoup::Overlap(
      o.rhs_id, o.rhs_begin, o.rhs_end,
      o.lhs_id, o.lhs_begin, o.lhs_end,
      o.score, cigar_alignment_reverse(o.alignment),
      o.strand);
}

inline extended_overlap cigar_extended_overlap_reverse(const extended_overlap &eo) {
  return extended_overlap{
      biosoup::Overlap(
          eo.overlap.rhs_id, eo.overlap.rhs_begin, eo.overlap.rhs_end,
          eo.overlap.lhs_id, eo.overlap.lhs_begin, eo.overlap.lhs_end,
          eo.overlap.score,
          eo.overlap.strand),
      edlib_align{
          eo.edlib_alignment.matches,
          eo.edlib_alignment.block_length,
          cigar_alignment_reverse(eo.edlib_alignment.cigar),
          eo.edlib_alignment.edit_distance
      }
  };
}

inline std::uint32_t overlap_length(const biosoup::Overlap &o) {
  return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
}

inline bool overlap_update(biosoup::Overlap &o, raven::Graph &graph) {
  if (graph.piles_[o.lhs_id]->is_invalid() ||
      graph.piles_[o.rhs_id]->is_invalid()) {
    return false;
  }
  if (o.lhs_begin >= graph.piles_[o.lhs_id]->end() ||
      o.lhs_end <= graph.piles_[o.lhs_id]->begin() ||
      o.rhs_begin >= graph.piles_[o.rhs_id]->end() ||
      o.rhs_end <= graph.piles_[o.rhs_id]->begin()) {
    return false;
  }

  std::uint32_t lhs_begin = o.lhs_begin + (o.strand ?
                                           (o.rhs_begin < graph.piles_[o.rhs_id]->begin() ?
                                            graph.piles_[o.rhs_id]->begin() - o.rhs_begin : 0)
                                                    :
                                           (o.rhs_end > graph.piles_[o.rhs_id]->end() ?
                                            o.rhs_end - graph.piles_[o.rhs_id]->end() : 0));
  std::uint32_t lhs_end = o.lhs_end - (o.strand ?
                                       (o.rhs_end > graph.piles_[o.rhs_id]->end() ?
                                        o.rhs_end - graph.piles_[o.rhs_id]->end() : 0)
                                                :
                                       (o.rhs_begin < graph.piles_[o.rhs_id]->begin() ?
                                        graph.piles_[o.rhs_id]->begin() - o.rhs_begin : 0));

  std::uint32_t rhs_begin = o.rhs_begin + (o.strand ?
                                           (o.lhs_begin < graph.piles_[o.lhs_id]->begin() ?
                                            graph.piles_[o.lhs_id]->begin() - o.lhs_begin : 0)
                                                    :
                                           (o.lhs_end > graph.piles_[o.lhs_id]->end() ?
                                            o.lhs_end - graph.piles_[o.lhs_id]->end() : 0));
  std::uint32_t rhs_end = o.rhs_end - (o.strand ?
                                       (o.lhs_end > graph.piles_[o.lhs_id]->end() ?
                                        o.lhs_end - graph.piles_[o.lhs_id]->end() : 0)
                                                :
                                       (o.lhs_begin < graph.piles_[o.lhs_id]->begin() ?
                                        graph.piles_[o.lhs_id]->begin() - o.lhs_begin : 0));

  if (lhs_begin >= graph.piles_[o.lhs_id]->end() ||
      lhs_end <= graph.piles_[o.lhs_id]->begin() ||
      rhs_begin >= graph.piles_[o.rhs_id]->end() ||
      rhs_end <= graph.piles_[o.rhs_id]->begin()) {
    return false;
  }

  lhs_begin = std::max(lhs_begin, graph.piles_[o.lhs_id]->begin());
  lhs_end = std::min(lhs_end, graph.piles_[o.lhs_id]->end());
  rhs_begin = std::max(rhs_begin, graph.piles_[o.rhs_id]->begin());
  rhs_end = std::min(rhs_end, graph.piles_[o.rhs_id]->end());

  if (lhs_begin >= lhs_end || lhs_end - lhs_begin < 84 ||
      rhs_begin >= rhs_end || rhs_end - rhs_begin < 84) {
    return false;
  }

  o.lhs_begin = lhs_begin;
  o.lhs_end = lhs_end;
  o.rhs_begin = rhs_begin;
  o.rhs_end = rhs_end;

  return true;
}

inline std::uint32_t extend_overlap(biosoup::Overlap &o,
                                    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences) {
  std::uint32_t start_offset = std::min(o.lhs_begin, o.rhs_begin);
  std::uint32_t end_offset = std::min(sequences[o.lhs_id]->inflated_len - o.lhs_end,
                                      sequences[o.rhs_id]->inflated_len - o.rhs_end);

  o.lhs_begin = o.lhs_begin - start_offset;
  o.lhs_end = o.lhs_end + end_offset;
  o.rhs_begin = o.rhs_begin - start_offset;
  o.rhs_end = o.rhs_end + end_offset;
}

inline std::uint32_t overlap_type(const biosoup::Overlap &o, raven::Graph &graph) {
  std::uint32_t lhs_length =
      graph.piles_[o.lhs_id]->end() - graph.piles_[o.lhs_id]->begin();
  std::uint32_t lhs_begin = o.lhs_begin - graph.piles_[o.lhs_id]->begin();
  std::uint32_t lhs_end = o.lhs_end - graph.piles_[o.lhs_id]->begin();

  std::uint32_t rhs_length =
      graph.piles_[o.rhs_id]->end() - graph.piles_[o.rhs_id]->begin();
  std::uint32_t rhs_begin = o.strand ?
                            o.rhs_begin - graph.piles_[o.rhs_id]->begin() :
                            rhs_length - (o.rhs_end - graph.piles_[o.rhs_id]->begin());
  std::uint32_t rhs_end = o.strand ?
                          o.rhs_end - graph.piles_[o.rhs_id]->begin() :
                          rhs_length - (o.rhs_begin - graph.piles_[o.rhs_id]->begin());

  std::uint32_t overhang =
      std::min(lhs_begin, rhs_begin) +
          std::min(lhs_length - lhs_end, rhs_length - rhs_end);

  if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
      rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
    return 0;  // internal
  }
  if (lhs_begin <= rhs_begin &&
      lhs_length - lhs_end <= rhs_length - rhs_end) {
    return 1;  // lhs contained
  }
  if (rhs_begin <= lhs_begin &&
      rhs_length - rhs_end <= lhs_length - lhs_end) {
    return 2;  // rhs contained
  }
  if (lhs_begin > rhs_begin) {
    return 3;  // lhs -> rhs
  }
  return 4;  // rhs -> lhs
}

inline bool overlap_finalize(biosoup::Overlap &o, raven::Graph &graph) {
  o.score = overlap_type(o, graph);
  if (o.score < 3) {
    return false;
  }

  o.lhs_begin -= graph.piles_[o.lhs_id]->begin();
  o.lhs_end -= graph.piles_[o.lhs_id]->begin();

  o.rhs_begin -= graph.piles_[o.rhs_id]->begin();
  o.rhs_end -= graph.piles_[o.rhs_id]->begin();
  if (!o.strand) {
    auto rhs_begin = o.rhs_begin;
    o.rhs_begin = graph.piles_[o.rhs_id]->length() - o.rhs_end;
    o.rhs_end = graph.piles_[o.rhs_id]->length() - rhs_begin;
  }
  return true;
}

std::vector<std::vector<std::uint32_t>> connected_components(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                                             std::vector<std::vector<biosoup::Overlap>> overlaps,
                                                             raven::Graph &graph) {  // NOLINT
  std::vector<std::vector<std::uint32_t>> connections(sequences.size());
  for (const auto &it : overlaps) {
    for (const auto &jt : it) {
      if (overlap_type(jt, graph) > 2) {
        connections[jt.lhs_id].emplace_back(jt.rhs_id);
        connections[jt.rhs_id].emplace_back(jt.lhs_id);
      }
    }
  }

  std::vector<std::vector<std::uint32_t>> dst;
  std::vector<char> is_visited(sequences.size(), false);
  for (std::uint32_t i = 0; i < connections.size(); ++i) {
    if (graph.piles_[i]->is_invalid() || is_visited[i]) {
      continue;
    }

    dst.resize(dst.size() + 1);
    std::deque<std::uint32_t> que = {i};
    while (!que.empty()) {
      std::uint32_t j = que.front();
      que.pop_front();

      if (is_visited[j]) {
        continue;
      }
      is_visited[j] = true;
      dst.back().emplace_back(j);

      for (const auto &it : connections[j]) {
        que.emplace_back(it);
      }
    }
  }

  return dst;
}

#endif //OVERLAP_HELPERS_H
