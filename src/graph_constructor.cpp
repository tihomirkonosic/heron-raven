
#include <cmath>
#include <deque>
#include <fstream>
#include "graph.hpp"
#include "graph_constructor.h"
#include "biosoup/overlap.hpp"
#include "edlib.h"
#include "biosoup/timer.hpp"
#include "extended_overlap.h"

#define STOP_AFTER_PILES 1

namespace raven {
  Graph_Constructor::Graph_Constructor(Graph &graph, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
      : graph_(graph), thread_pool_(thread_pool ?
                                    thread_pool :
                                    std::make_shared<thread_pool::ThreadPool>(1)) {
  }


  void Graph_Constructor::Construct(
      std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,  // NOLINT
      double disagreement,
      unsigned split,
      std::size_t kMaxNumOverlaps,
			std::uint8_t ploidy,
      std::uint8_t kmer_len,
      std::uint8_t window_len,
      std::uint16_t bandwidth,
      std::uint16_t chain_n,
      std::uint16_t match_n,
      std::uint16_t gap_size,
      double freq,
			bool hpc,
      bool paf,
      std::uint16_t valid_region_length_threshold,
      std::uint16_t valid_region_coverage_threshold,
      std::string herro_snps_path,
      std::string load_paf) {

    disagreement_ = disagreement;

    if (sequences.empty()) {
      return;
    }

    std::vector<std::vector<biosoup::Overlap>> overlaps(sequences.size());
    std::vector<std::vector<extended_overlap>> extended_overlaps(sequences.size());
//region biosoup::Overlap helper functions

    auto edlib_alignment_reverse = [](const std::string &s) -> std::string {
      std::string rs;
      if (s.empty()) {
        return rs;
      }
      for (char c: s) {
        switch (c) {
          case '2':
            rs += '\001';
            break;
          case '1':
            rs += '\002';
            break;
          default:
            rs += c;
            break;
        }
      }
      return rs;
    };

    auto cigar_alignment_reverse = [](const std::string &s) -> std::string {
      std::string rs;
      if (s.empty()) {
        return rs;
      }
      for (char c: s) {
        switch (c) {
          case 'I':
            rs += 'D';
            break;
          case 'D':
            rs += 'I';
            break;
          default:
            rs += c;
            break;
        }
      }
      return rs;
    };

    auto cigar_to_edlib_alignment = [](const std::string &s) -> std::string {
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
    };

    auto overlap_reverse = [&edlib_alignment_reverse](const biosoup::Overlap &o) -> biosoup::Overlap {
      return biosoup::Overlap(
          o.rhs_id, o.rhs_begin, o.rhs_end,
          o.lhs_id, o.lhs_begin, o.lhs_end,
          o.score, edlib_alignment_reverse(o.alignment),
          o.strand);
    };

    auto cigar_overlap_reverse = [&cigar_alignment_reverse](const biosoup::Overlap &o) -> biosoup::Overlap {
      return biosoup::Overlap(
          o.rhs_id, o.rhs_begin, o.rhs_end,
          o.lhs_id, o.lhs_begin, o.lhs_end,
          o.score, cigar_alignment_reverse(o.alignment),
          o.strand);
    };

    auto cigar_extended_overlap_reverse = [&cigar_alignment_reverse](const extended_overlap &eo) -> extended_overlap {
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
    };

    auto overlap_length = [](const biosoup::Overlap &o) -> std::uint32_t {
      return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
    };

    auto overlap_update = [&](biosoup::Overlap &o) -> bool {
      if (graph_.piles_[o.lhs_id]->is_invalid() ||
          graph_.piles_[o.rhs_id]->is_invalid()) {
        return false;
      }
      if (o.lhs_begin >= graph_.piles_[o.lhs_id]->end() ||
          o.lhs_end <= graph_.piles_[o.lhs_id]->begin() ||
          o.rhs_begin >= graph_.piles_[o.rhs_id]->end() ||
          o.rhs_end <= graph_.piles_[o.rhs_id]->begin()) {
        return false;
      }

      std::uint32_t lhs_begin = o.lhs_begin + (o.strand ?
                                               (o.rhs_begin < graph_.piles_[o.rhs_id]->begin() ?
                                                graph_.piles_[o.rhs_id]->begin() - o.rhs_begin : 0)
                                                        :
                                               (o.rhs_end > graph_.piles_[o.rhs_id]->end() ?
                                                o.rhs_end - graph_.piles_[o.rhs_id]->end() : 0));
      std::uint32_t lhs_end = o.lhs_end - (o.strand ?
                                           (o.rhs_end > graph_.piles_[o.rhs_id]->end() ?
                                            o.rhs_end - graph_.piles_[o.rhs_id]->end() : 0)
                                                    :
                                           (o.rhs_begin < graph_.piles_[o.rhs_id]->begin() ?
                                            graph_.piles_[o.rhs_id]->begin() - o.rhs_begin : 0));

      std::uint32_t rhs_begin = o.rhs_begin + (o.strand ?
                                               (o.lhs_begin < graph_.piles_[o.lhs_id]->begin() ?
                                                graph_.piles_[o.lhs_id]->begin() - o.lhs_begin : 0)
                                                        :
                                               (o.lhs_end > graph_.piles_[o.lhs_id]->end() ?
                                                o.lhs_end - graph_.piles_[o.lhs_id]->end() : 0));
      std::uint32_t rhs_end = o.rhs_end - (o.strand ?
                                           (o.lhs_end > graph_.piles_[o.lhs_id]->end() ?
                                            o.lhs_end - graph_.piles_[o.lhs_id]->end() : 0)
                                                    :
                                           (o.lhs_begin < graph_.piles_[o.lhs_id]->begin() ?
                                            graph_.piles_[o.lhs_id]->begin() - o.lhs_begin : 0));

      if (lhs_begin >= graph_.piles_[o.lhs_id]->end() ||
          lhs_end <= graph_.piles_[o.lhs_id]->begin() ||
          rhs_begin >= graph_.piles_[o.rhs_id]->end() ||
          rhs_end <= graph_.piles_[o.rhs_id]->begin()) {
        return false;
      }

      lhs_begin = std::max(lhs_begin, graph_.piles_[o.lhs_id]->begin());
      lhs_end = std::min(lhs_end, graph_.piles_[o.lhs_id]->end());
      rhs_begin = std::max(rhs_begin, graph_.piles_[o.rhs_id]->begin());
      rhs_end = std::min(rhs_end, graph_.piles_[o.rhs_id]->end());

      if (lhs_begin >= lhs_end || lhs_end - lhs_begin < 84 ||
          rhs_begin >= rhs_end || rhs_end - rhs_begin < 84) {
        return false;
      }

      o.lhs_begin = lhs_begin;
      o.lhs_end = lhs_end;
      o.rhs_begin = rhs_begin;
      o.rhs_end = rhs_end;

      return true;
    };

    auto extend_overlap = [&] (biosoup::Overlap &o) -> std::uint32_t {
      std::uint32_t start_offset = std::min(o.lhs_begin, o.rhs_begin);
      std::uint32_t end_offset = std::min(sequences[o.lhs_id]->inflated_len - o.lhs_end,
                                          sequences[o.rhs_id]->inflated_len - o.rhs_end);

      
      o.lhs_begin = o.lhs_begin - start_offset;
      o.lhs_end = o.lhs_end + end_offset;
      o.rhs_begin = o.rhs_begin - start_offset;
      o.rhs_end = o.rhs_end + end_offset;
    };

    auto overlap_type = [&](const biosoup::Overlap &o) -> std::uint32_t {
      std::uint32_t lhs_length =
          graph_.piles_[o.lhs_id]->end() - graph_.piles_[o.lhs_id]->begin();
      std::uint32_t lhs_begin = o.lhs_begin - graph_.piles_[o.lhs_id]->begin();
      std::uint32_t lhs_end = o.lhs_end - graph_.piles_[o.lhs_id]->begin();

      std::uint32_t rhs_length =
          graph_.piles_[o.rhs_id]->end() - graph_.piles_[o.rhs_id]->begin();
      std::uint32_t rhs_begin = o.strand ?
                                o.rhs_begin - graph_.piles_[o.rhs_id]->begin() :
                                rhs_length - (o.rhs_end - graph_.piles_[o.rhs_id]->begin());
      std::uint32_t rhs_end = o.strand ?
                              o.rhs_end - graph_.piles_[o.rhs_id]->begin() :
                              rhs_length - (o.rhs_begin - graph_.piles_[o.rhs_id]->begin());

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
    };

    auto overlap_finalize = [&](biosoup::Overlap &o) -> bool {
      o.score = overlap_type(o);
      if (o.score < 3) {
        return false;
      }

      o.lhs_begin -= graph_.piles_[o.lhs_id]->begin();
      o.lhs_end -= graph_.piles_[o.lhs_id]->begin();

      o.rhs_begin -= graph_.piles_[o.rhs_id]->begin();
      o.rhs_end -= graph_.piles_[o.rhs_id]->begin();
      if (!o.strand) {
        auto rhs_begin = o.rhs_begin;
        o.rhs_begin = graph_.piles_[o.rhs_id]->length() - o.rhs_end;
        o.rhs_end = graph_.piles_[o.rhs_id]->length() - rhs_begin;
      }
      return true;
    };

    auto connected_components = [&]() -> std::vector<std::vector<std::uint32_t>> {  // NOLINT
      std::vector<std::vector<std::uint32_t>> connections(sequences.size());
      for (const auto &it: overlaps) {
        for (const auto &jt: it) {
          if (overlap_type(jt) > 2) {
            connections[jt.lhs_id].emplace_back(jt.rhs_id);
            connections[jt.rhs_id].emplace_back(jt.lhs_id);
          }
        }
      }

      std::vector<std::vector<std::uint32_t>> dst;
      std::vector<char> is_visited(sequences.size(), false);
      for (std::uint32_t i = 0; i < connections.size(); ++i) {
        if (graph_.piles_[i]->is_invalid() || is_visited[i]) {
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

          for (const auto &it: connections[j]) {
            que.emplace_back(it);
          }
        }
      }

      return dst;
    };

//endregion

    graph_.annotations_.resize(sequences.size());
    struct base_pile {
      std::uint32_t a;
      std::uint32_t c;
      std::uint32_t g;
      std::uint32_t t;
      std::uint32_t i;
      std::uint32_t d;
    };

//region annotations_ helper functions

    auto edlib_wrapper = [&](
        std::uint32_t i,
        const biosoup::Overlap &it,
        const std::string &lhs,
        const std::string &rhs) -> edlib_align {
      edlib_align alignment_result;
      EdlibAlignResult result = edlibAlign(
          lhs.c_str(), lhs.size(),
          rhs.c_str(), rhs.size(),
          edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0)); // align lhs and rhs
      if (result.status == EDLIB_STATUS_OK) {
        alignment_result.cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
        alignment_result.edit_distance = result.editDistance;
        alignment_result.block_length= result.alignmentLength;
        alignment_result.matches = 0;
        for (int i = 0; i < result.alignmentLength; i++) {
          if (result.alignment[i] == 0) {
            ++alignment_result.matches;
          }
        }
        return alignment_result;
      } else {
        edlibFreeAlignResult(result);
        alignment_result.cigar = "";
        alignment_result.matches = 0;
        alignment_result.edit_distance = -1;
        alignment_result.block_length = 0;
        return alignment_result;
      }
    };

    auto find_pairwise_alignment = [&](std::uint32_t i, std::vector<extended_overlap> &ovlps) -> void {
      for(auto &it : ovlps){

        std::string lhs = sequences[i]->InflateData(it.overlap.lhs_begin, it.overlap.lhs_end - it.overlap.lhs_begin);
        biosoup::NucleicAcid rhs_{"",
                          sequences[it.overlap.rhs_id]->InflateData(it.overlap.rhs_begin, it.overlap.rhs_end - it.overlap.rhs_begin)};
                                                              
        if (!it.overlap.strand) rhs_.ReverseAndComplement();

        auto rhs = rhs_.InflateData();

        it.edlib_alignment = edlib_wrapper(i, it.overlap, lhs, rhs);
      }
    };

    auto call_snps = [&](std::uint32_t i, std::vector<extended_overlap> ovlps_final) -> void {
      std::uint32_t seq_inflated_len = sequences[i]->inflated_len;
      std::vector<base_pile> base_pile_tmp(seq_inflated_len);
      std::vector<std::uint32_t> cov;

      for (auto &ovlp: ovlps_final) {
        if (!(ovlp.edlib_alignment.cigar.empty())) {
          std::uint32_t lhs_pos = ovlp.overlap.lhs_begin;
          std::uint32_t rhs_pos = 0;
          biosoup::NucleicAcid rhs{"",
                                   sequences[ovlp.overlap.rhs_id]->InflateData(ovlp.overlap.rhs_begin,
                                                                       ovlp.overlap.rhs_end - ovlp.overlap.rhs_begin)};
          if (!ovlp.overlap.strand) rhs.ReverseAndComplement();
          std::string rhs_tmp = rhs.InflateData();

          std::string edlib_alignment = cigar_to_edlib_alignment(ovlp.edlib_alignment.cigar);
          for (auto &edlib_align: edlib_alignment) {
            switch (edlib_align) {
              case 0:
              case 3: {
                switch (rhs_tmp[rhs_pos]) {
                  case 'A':
                    ++base_pile_tmp[lhs_pos].a;
                    break;
                  case 'C':
                    ++base_pile_tmp[lhs_pos].c;
                    break;
                  case 'G':
                    ++base_pile_tmp[lhs_pos].g;
                    break;
                  case 'T':
                    ++base_pile_tmp[lhs_pos].t;
                    break;
                  default:
                    break; // if they align
                }
                ++lhs_pos;
                ++rhs_pos;
                break;
              }
              case 1: {
                ++base_pile_tmp[lhs_pos].i;
                ++lhs_pos;
                break; // insertion on the left hand side
              }
              case 2: {
                if (!(lhs_pos >= base_pile_tmp.size())) {
                  ++base_pile_tmp[lhs_pos].d;
                }
                ++rhs_pos;
                break; // deletion on the left hand side
              }
              default:
                break;
            }
          }
        }
      }


      for (const auto &jt: base_pile_tmp) {
        //cov.emplace_back(jt.a + jt.c + jt.g + jt.t);
        cov.emplace_back(jt.a + jt.c + jt.g + jt.t + jt.d + jt.i);
      }

      std::nth_element(cov.begin(), cov.begin() + cov.size() / 2, cov.end());
      double m = cov[cov.size() / 2] * 2. / 3.;

      std::size_t j = 0;
      for (const auto &jt: base_pile_tmp) {
        std::vector<double> counts = {
            static_cast<double>(jt.a),
            static_cast<double>(jt.c),
            static_cast<double>(jt.g),
            static_cast<double>(jt.t),
            static_cast<double>(jt.d),
            static_cast<double>(jt.i)
        };

        double sum = std::accumulate(counts.begin(), counts.end(), 0);

        if (use_frequencies) {
          for (auto &kt: counts) {
            kt /= sum;
          };
        };

        if (sum > m) {
          std::size_t variants = 0;
          for (const auto &it: counts) {
            if (use_frequencies) {
              if (freq_low_th < it && it < freq_high_th) {
                ++variants;
              }
            } else {
              if (it > variant_call_th) {
                ++variants;
              }
            }
          }
          if (variants > 1) graph_.annotations_[i].emplace(j);


        };
        ++j;
      };
    };

    auto annotation_extract = [&](
        std::uint32_t i,
        std::uint32_t begin,
        std::uint32_t end,
        std::uint32_t len,
        bool strand) -> std::unordered_set<std::uint32_t> {
      std::unordered_set<std::uint32_t> dst;
      if (graph_.annotations_[i].empty()) {
        return dst;
      }
      for (const auto &it: graph_.annotations_[i]) {
        if (begin <= it && it <= end) {
          dst.emplace(strand ? it : len - 1 - it);
        }
      }
      return dst;
    };

//endregion

    // checkpoint test
    if (graph_.stage() == Graph_Stage::Construct_Graph && graph_.use_checkpoints()) {
      graph_.Store();
    }

    biosoup::Timer timer{};

    ram::MinimizerEngine minimizer_engine{
        thread_pool_,
        kmer_len,
        window_len,
        bandwidth,
        chain_n,
        match_n,
        gap_size
    };

    // find overlaps and create piles
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      for (const auto &it: sequences) {
        graph_.piles_.emplace_back(new Pile(it->id, it->inflated_len));
      }
      
      if(!load_paf.empty()){
        LoadOverlaps(load_paf, sequences, extended_overlaps);
        std::vector<std::future<void>> void_futures;
        for(int i = 0; i < sequences.size(); i++){
          void_futures.emplace_back(thread_pool_->Submit(
              [&](std::uint32_t i) -> void {
                find_pairwise_alignment(i, extended_overlaps[i]);
              },
              i));
        };
        for (const auto &it: void_futures) {
          it.wait();
        }

        void_futures.clear();

        
        std::vector<std::future<void>> extended_layers_futures;
        std::uint16_t counter = 0;
        for (const auto &it: graph_.piles_) {
          counter += 1;
          //std::cerr << counter << std::endl;
          if (extended_overlaps[it->id()].empty()){
                    continue;
                  }
          extended_layers_futures.emplace_back(thread_pool_->Submit(
           [&](std::uint32_t i) -> void {
              it->AddExtendedLayers(
              extended_overlaps[it->id()].begin(),
              extended_overlaps[it->id()].end());

              // if (extended_overlaps[it->id()].size() < kMaxNumOverlaps) {
              //   return;
              // }

              // std::sort(extended_overlaps[i].begin(), extended_overlaps[i].end(),
              //           [&](const extended_overlap &lhs,
              //               const extended_overlap &rhs) -> bool {
              //             return overlap_length(lhs.overlap) > overlap_length(rhs.overlap);
              //           });

              // std::vector<extended_overlap> tmp;
              // tmp.insert(tmp.end(), extended_overlaps[i].begin(),
              //             extended_overlaps[i].begin() + kMaxNumOverlaps);  // NOLINT
              // tmp.swap(extended_overlaps[i]);
           },
           it->id()));
        }
        for (const auto &it: extended_layers_futures) {
          it.wait();
        }
      }
       else {
        std::size_t bytes = 0;
        for (std::uint32_t i = 0, j = 0; i < sequences.size(); ++i) {
          bytes += sequences[i]->inflated_len;
          if (i != sequences.size() - 1 && bytes < (1ULL << 32)) {
            continue;
          }
          bytes = 0;

          timer.Start();

          minimizer_engine.Minimize(
              sequences.begin() + j,
              sequences.begin() + i + 1,
              true);
          minimizer_engine.Filter(freq);

          std::cerr << "[raven::Graph::Construct] minimized "
                    << j << " - " << i + 1 << " / " << sequences.size() << " "
                    << std::fixed << timer.Stop() << "s"
                    << std::endl;

          timer.Start();

          std::vector<std::uint32_t> num_overlaps(overlaps.size());
          for (std::uint32_t k = 0; k < overlaps.size(); ++k) {
            num_overlaps[k] = overlaps[k].size();
          }

          std::vector<std::future<std::vector<extended_overlap>>> thread_futures;

          for (std::uint32_t k = 0; k < i + 1; ++k) {
            thread_futures.emplace_back(thread_pool_->Submit(
                [&](std::uint32_t i) -> std::vector<extended_overlap> { // map sequences
                  std::vector<biosoup::Overlap> ovlps = minimizer_engine.Map(sequences[i], true, true,
                                                                            true, hpc);
                  if (!ovlps.empty()) {
                    std::vector<extended_overlap> ovlps_final;

                    std::sort(ovlps.begin(), ovlps.end(),
                              [&](const biosoup::Overlap &lhs,
                                  const biosoup::Overlap &rhs) -> bool {
                                return overlap_length(lhs) > overlap_length(rhs);
                              });


                    std::vector<biosoup::Overlap> tmp;
                    // tmp.insert(tmp.end(), ovlps.begin(),
                    //            ovlps.begin() + (ovlps.size() > 120 ? 120 : ovlps.size()));  // NOLINT // factor in coverage
                    // tmp.swap(ovlps);

                    for (auto &ovlp: ovlps) {
                    if(overlap_length(ovlp) > 500){
                        ovlp.lhs_begin = ovlp.lhs_begin - (window_len+kmer_len-1) ? ovlp.lhs_begin - (window_len+kmer_len-1) : 0;
                        ovlp.lhs_end = ovlp.lhs_end + (window_len+kmer_len-1) < sequences[ovlp.lhs_id]->inflated_len ? ovlp.lhs_end + (window_len+kmer_len-1) : sequences[ovlp.lhs_id]->inflated_len;

                        ovlp.rhs_begin = ovlp.rhs_begin - (window_len+kmer_len-1) ? ovlp.rhs_begin - (window_len+kmer_len-1) : 0;
                        ovlp.rhs_end = ovlp.rhs_end + (window_len+kmer_len-1) < sequences[ovlp.rhs_id]->inflated_len ? ovlp.rhs_end + (window_len+kmer_len-1) : sequences[ovlp.rhs_id]->inflated_len;

                        auto lhs = sequences[i]->InflateData(ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);
                                                            
                        biosoup::NucleicAcid rhs_{"",
                                                  sequences[ovlp.rhs_id]->InflateData(ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin)};
                                                                                      

                        if (!ovlp.strand) rhs_.ReverseAndComplement();

                        auto rhs = rhs_.InflateData();

                        
                        edlib_align tmp = edlib_wrapper(i, ovlp, lhs, rhs);
                        if(static_cast<float>(tmp.matches) / tmp.block_length > 0.9){
                         // edlib_align tmp;
                          biosoup::Overlap ovlp_tmp{ovlp.lhs_id, ovlp.lhs_begin, ovlp.lhs_end,
                                                    ovlp.rhs_id, ovlp.rhs_begin, ovlp.rhs_end, 
                                                    ovlp.score, ovlp.strand};

                          extended_overlap total_ovlp{ovlp_tmp, tmp};
                          ovlps_final.emplace_back(total_ovlp);
                    }
                  }

                    };
                    return ovlps_final;
                  }
                  std::vector<extended_overlap> total_ovlps;
                  return total_ovlps;
                },
                k));

            bytes += sequences[k]->inflated_len;
            if (k != i && bytes < (1U << 30)) {
              continue;
            }
            bytes = 0;

            for (auto &it: thread_futures) {
              for (const auto &jt: it.get()) {
                extended_overlaps[jt.overlap.lhs_id].emplace_back(jt);
                //overlaps.emplace_back(jt.overlap);
                extended_overlaps[jt.overlap.rhs_id].emplace_back(cigar_extended_overlap_reverse(jt));
                //overlaps.emplace_back(overlap_reverse(jt.overlap));
              }
            }
            thread_futures.clear();
          }

        std::vector<std::future<void>> void_futures;
        for (const auto &it: graph_.piles_) {
          if (extended_overlaps[it->id()].empty()            
              || extended_overlaps[it->id()].size() == num_overlaps[it->id()]
              ){
            continue;
          }

          void_futures.emplace_back(thread_pool_->Submit(
              [&](std::uint32_t i) -> void {

                graph_.piles_[i]->AddExtendedLayers(
                    extended_overlaps[i].begin(),
                    extended_overlaps[i].end());

                // num_overlaps[i] = std::min(
                //     extended_overlaps[i].size(),
                //     kMaxNumOverlaps);

                // if (extended_overlaps[i].size() < kMaxNumOverlaps) {
                //   return;
                // }

                // std::sort(extended_overlaps[i].begin(), extended_overlaps[i].end(),
                //           [&](const extended_overlap &lhs,
                //               const extended_overlap &rhs) -> bool {
                //             return overlap_length(lhs.overlap) > overlap_length(rhs.overlap);
                //           });

                // std::vector<extended_overlap> tmp;
                // tmp.insert(tmp.end(), extended_overlaps[i].begin(),
                //             extended_overlaps[i].begin() + kMaxNumOverlaps);  // NOLINT
                // tmp.swap(extended_overlaps[i]);
              },
              it->id()));
        }
        for (const auto &it: void_futures) {
          it.wait();
        }
      


      std::cerr << "[raven::Graph::Construct] mapped sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

        j = i + 1;
      }
    }
    }

      graph_.PrintOverlaps(extended_overlaps, sequences, true, "initial.paf");
      std::cerr << "[raven::Graph::Construct] initial overlaps printed"
               << std::endl;

    if(herro_snps_path == ""){
      std::vector<std::future<void>> void_futures;
      for(int i = 0; i < sequences.size(); i++){
        void_futures.emplace_back(thread_pool_->Submit(
            [&](std::uint32_t i) -> void {
              call_snps(i, extended_overlaps[i]);
            },
            i));
      };

      for (const auto &it: void_futures) {
        it.wait();
      }
    
      std::cerr << "[raven::Graph::Construct] snps called"
                << std::endl;

      if (print_snp_data && ploidy >= 2) {
        std::ofstream outdata;
        outdata.open("snp_annotations.anno");

        for (std::uint32_t i = 0; i < graph_.annotations_.size(); ++i) {
          if (graph_.annotations_[i].empty()) {
            continue;
          }
          outdata << sequences[i]->name << " ";
          for (const auto &jt: graph_.annotations_[i]) {
            outdata << " " << jt;
          }
          outdata << std::endl;
        }

      };
    } else {
      std::cerr << "Loading snps" << std::endl;
      LoadHerroSNPs(herro_snps_path, sequences);
      std::ofstream outdata;
      outdata.open("snp_annotations_check.anno");
      for (std::uint32_t i = 0; i < graph_.annotations_.size(); ++i) {
        if (graph_.annotations_[i].empty()) {
          continue;
        }
        outdata << sequences[i]->name << " ";
        for (const auto &jt: graph_.annotations_[i]) {
          outdata << " " << jt;
        }
        outdata << std::endl;
      }
    };
    
    std::ofstream outdata;
    outdata.open("piles.csv");
    for(int i = 0; i < graph_.piles_.size(); i++){
      outdata << graph_.piles_[i]->id() << ",";
      std::vector<std::uint16_t> coverages = graph_.piles_[i]->get_data();
      for(auto &element: coverages){
        outdata << element << ";";
      }
      outdata << std::endl;

    }

    if(STOP_AFTER_PILES == 1){
      exit(0);
    }
    
    graph_.PrintOverlaps(extended_overlaps, sequences, true, "beforeValidRegion.paf");

    // trim and annotate piles
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      timer.Start();
      std::vector<std::future<void>> thread_futures;
      for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->Submit(
        [&](std::uint32_t i) -> void {
              graph_.piles_[i]->FindValidRegion(valid_region_coverage_threshold, valid_region_length_threshold);
              if (graph_.piles_[i]->is_invalid()) {
                std::vector<extended_overlap>().swap(extended_overlaps[i]);
                
              } else {
                graph_.piles_[i]->FindMedian();
                graph_.piles_[i]->FindChimericRegions();
              }
           },
           i));
      }
      for (const auto &it: thread_futures) {
        it.wait();
      }
      thread_futures.clear();

      std::cerr << "[raven::Graph::Construct] annotated piles "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }

    // std::ofstream outdata;
    // outdata.open("piles.csv");
    // for(int i = 0; i < graph_.piles_.size(); i++){
    //   outdata << graph_.piles_[i]->id() << ",";
    //   std::vector<std::uint16_t> coverages = graph_.piles_[i]->get_data();
    //   for(auto &element: coverages){
    //     outdata << element << ";";
    //   }
    //   outdata << std::endl;

    // }

    // resolve contained reads
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      timer.Start();
      std::vector<std::future<void>> futures;
      for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
         futures.emplace_back(thread_pool_->Submit(
             [&](std::uint32_t i) -> void {
              std::uint32_t k = 0;
              for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {
                if (!overlap_update(extended_overlaps[i][j].overlap)) {
                  continue;
                }

                const auto &it = extended_overlaps[i][j];

                auto lhs_anno = annotation_extract(
                    it.overlap.lhs_id,
                    it.overlap.lhs_begin,
                    it.overlap.lhs_end,
                    sequences[it.overlap.lhs_id]->inflated_len,
                    true);

                auto rhs_anno = annotation_extract(
                    it.overlap.rhs_id,
                    it.overlap.rhs_begin,
                    it.overlap.rhs_end,
                    sequences[it.overlap.rhs_id]->inflated_len,
                    it.overlap.strand);

                if (!lhs_anno.empty() || !rhs_anno.empty()) {
                  //std::vector<std::pair<char, int>> cigar = parse_cigar_string(it.alignment);
                  std::string edlib_alignment = cigar_to_edlib_alignment(it.edlib_alignment.cigar);
                  std::uint32_t lhs_pos = it.overlap.lhs_begin;
                  std::uint32_t rhs_pos = it.overlap.strand ?
                                          it.overlap.rhs_begin :
                                          sequences[it.overlap.rhs_id]->inflated_len - it.overlap.rhs_end;

                  std::uint32_t mismatches = 0;
                  std::uint32_t snps = 0;

                  for (int t = 0; t < static_cast<std::int32_t>(edlib_alignment.length()); t++) {
                    if (lhs_anno.find(lhs_pos) != lhs_anno.end() ||
                        rhs_anno.find(rhs_pos) != rhs_anno.end()) {
                      ++snps;
                      if (edlib_alignment[t] == 3) {
                        ++mismatches;
                      }
                    }
                    switch (edlib_alignment[t]) {
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

                  if (mismatches / static_cast<double>(snps) > disagreement_) {
                    continue;
                  }
                }
                
                extended_overlaps[i][k++] = extended_overlaps[i][j];
              }
              extended_overlaps[i].resize(k);
             },
             i));
      }
       for (const auto &it: futures) {
         it.wait();
       }
      futures.clear();
      graph_.PrintOverlaps(extended_overlaps, sequences, true, "afterSnps.paf");

      for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
        std::uint32_t k = 0;
        for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
          if (!overlap_update(overlaps[i][j])) {
            continue;
          }
          std::uint32_t type = overlap_type(overlaps[i][j]);
          if (type == 1 &&
              !graph_.piles_[overlaps[i][j].rhs_id]->is_maybe_chimeric()) {
            graph_.piles_[i]->set_is_contained();
          } else if (type == 2 &&
                     !graph_.piles_[i]->is_maybe_chimeric()) {
            graph_.piles_[overlaps[i][j].rhs_id]->set_is_contained();
          } else {
            overlaps[i][k++] = overlaps[i][j];
          }
        }
        overlaps[i].resize(k);
      }
      for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
        if (graph_.piles_[i]->is_contained()) {
          graph_.piles_[i]->set_is_invalid();
          std::vector<biosoup::Overlap>().swap(overlaps[i]);
        }
      }

      std::cerr << "[raven::Graph::Construct] removed contained sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }

    graph_.PrintOverlaps(extended_overlaps, sequences, true, "afterContained.paf");

    // resolve chimeric sequences
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      timer.Start();

      while (true) {
        auto components = connected_components();
        for (const auto &it: components) {
          std::vector<std::uint16_t> medians;
          for (const auto &jt: it) {
            medians.emplace_back(graph_.piles_[jt]->median());
          }
          std::nth_element(
              medians.begin(),
              medians.begin() + medians.size() / 2,
              medians.end());
          std::uint16_t median = medians[medians.size() / 2];

          std::vector<std::future<void>> thread_futures;
          for (const auto &jt: it) {
            thread_futures.emplace_back(thread_pool_->Submit(
                [&](std::uint32_t i) -> void {
                  graph_.piles_[i]->ClearChimericRegions(median);
                  if (graph_.piles_[i]->is_invalid()) {
                    std::vector<biosoup::Overlap>().swap(overlaps[i]);
                  }
                },
                jt));
          }
          for (const auto &it: thread_futures) {
            it.wait();
          }
          thread_futures.clear();
        }

        bool is_changed = false;
        for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
          std::uint32_t k = 0;
          for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
            if (overlap_update(overlaps[i][j])) {
              overlaps[i][k++] = overlaps[i][j];
            } else {
              is_changed = true;
            }
          }
          overlaps[i].resize(k);
        }

        if (!is_changed) {
          for (const auto &it: overlaps) {
            for (const auto &jt: it) {
              std::uint32_t type = overlap_type(jt);
              if (type == 1) {
                graph_.piles_[jt.lhs_id]->set_is_contained();
                graph_.piles_[jt.lhs_id]->set_is_invalid();
              } else if (type == 2) {
                graph_.piles_[jt.rhs_id]->set_is_contained();
                graph_.piles_[jt.rhs_id]->set_is_invalid();
              }
            }
          }
          overlaps.clear();
          break;
        }
      }

      std::cerr << "[raven::Graph::Construct] removed chimeric sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }


    // checkpoint
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      graph_.advance_stage();
      if (graph_.use_checkpoints()) {
        timer.Start();
        graph_.Store();
        std::cerr << "[raven::Graph::Construct] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }

    // find overlaps and repetitive regions
    if (graph_.stage() == Graph_Stage::Construct_Graph_2) {
      std::sort(sequences.begin(), sequences.end(),
                [&](const std::unique_ptr<biosoup::NucleicAcid> &lhs,
                    const std::unique_ptr<biosoup::NucleicAcid> &rhs) -> bool {
                  return graph_.piles_[lhs->id]->is_invalid() < graph_.piles_[rhs->id]->is_invalid() ||
                         // NOLINT
                         (graph_.piles_[lhs->id]->is_invalid() == graph_.piles_[rhs->id]->is_invalid() &&
                          lhs->id < rhs->id);  // NOLINT
                });

      std::vector<std::uint32_t> sequences_map(sequences.size());
      for (std::uint32_t i = 0; i < sequences.size(); ++i) {
        sequences_map[sequences[i]->id] = i;
      }

      std::uint32_t s = 0;
      for (std::uint32_t i = 0; i < sequences.size(); ++i) {
        if (graph_.piles_[sequences[i]->id]->is_invalid()) {
          s = i;
          break;
        }
      }

      // map valid reads to each other
      overlaps.resize(sequences.size() + 1);
      std::size_t bytes = 0;
      for (std::uint32_t i = 0, j = 0; i < s; ++i) {
        bytes += sequences[i]->inflated_len;
        if (i != s - 1 && bytes < (1U << 30)) {
          continue;
        }
        bytes = 0;

        timer.Start();

        minimizer_engine.Minimize(
            sequences.begin() + j,
            sequences.begin() + i + 1);

        std::cerr << "[raven::Graph::Construct] minimized "
                  << j << " - " << i + 1 << " / " << s << " "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;

        timer.Start();

        std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;
        minimizer_engine.Filter(0.001);
        for (std::uint32_t k = 0; k < i + 1; ++k) {
          thread_futures.emplace_back(thread_pool_->Submit(
              [&](std::uint32_t i) -> std::vector<biosoup::Overlap> {
                std::vector<std::uint32_t> filtered;
                auto dst = minimizer_engine.Map(
                    sequences[i],
                    true,  // avoid equal
                    true,  // avoid symmetric
                    false,  // minhash
                    hpc,
                    &filtered);
                graph_.piles_[sequences[i]->id]->AddKmers(filtered, kmer_len, sequences[i]); // NOLINT

                std::uint32_t k = 0;
                for (std::uint32_t j = 0; j < dst.size(); ++j) {
                  if (!overlap_update(dst[j]) || overlap_length(dst[j]) < 1000) {
                    continue;
                  }

                  const auto &it = dst[j];

                  auto lhs_anno = annotation_extract(
                      it.lhs_id,
                      it.lhs_begin,
                      it.lhs_end,
                      sequences[sequences_map[it.lhs_id]]->inflated_len,
                      true);

                  auto rhs_anno = annotation_extract(
                      it.rhs_id,
                      it.rhs_begin,
                      it.rhs_end,
                      sequences[sequences_map[it.rhs_id]]->inflated_len,
                      it.strand);

                  if (!lhs_anno.empty() || !rhs_anno.empty()) {
                    auto lhs = sequences[sequences_map[it.lhs_id]]->InflateData(
                        it.lhs_begin,
                        it.lhs_end - it.lhs_begin);

                    auto rhs = sequences[sequences_map[it.rhs_id]]->InflateData(
                        it.rhs_begin,
                        it.rhs_end - it.rhs_begin);
                    if (!it.strand) {
                      biosoup::NucleicAcid rhs_{"", rhs};
                      rhs_.ReverseAndComplement();
                      rhs = rhs_.InflateData();
                    }

                    EdlibAlignResult result = edlibAlign(
                        lhs.c_str(), lhs.size(),
                        rhs.c_str(), rhs.size(),
                        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr,
                                            0));  // NOLINT

                    if (result.status == EDLIB_STATUS_OK) {
                      std::string edlib_alignment = it.alignment;
                      std::uint32_t lhs_pos = it.lhs_begin;
                      std::uint32_t rhs_pos = it.strand ?
                                              it.rhs_begin :
                                              sequences[sequences_map[it.rhs_id]]->inflated_len -
                                              it.rhs_end;  // NOLINT

                      std::uint32_t mismatches = 0;
                      std::uint32_t snps = 0;

                      for (int a = 0; a < result.alignmentLength; ++a) {
                        //    if(result.alignment[a] != edlib_alignment[a]) std::cerr << "disagreement on pos " << a << std::endl;
                        if (lhs_anno.find(lhs_pos) != lhs_anno.end() ||
                            rhs_anno.find(rhs_pos) != rhs_anno.end()) {
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

                      if (mismatches / static_cast<double>(snps) > disagreement_) {
                        continue;
                      }
                    }
                  }

                  dst[k++] = dst[j];
                }
                dst.resize(k);

                return dst;
              },
              k));
        }
        for (auto &it: thread_futures) {
          for (auto &jt: it.get()) {
            if (!overlap_update(jt)) {
              continue;
            }
            std::uint32_t type = overlap_type(jt);
            if (type == 0) {
              continue;
            } else if (type == 1) {
              graph_.piles_[jt.lhs_id]->set_is_contained();
            } else if (type == 2) {
              graph_.piles_[jt.rhs_id]->set_is_contained();
            } else {
              if (overlaps.back().size() &&
                  overlaps.back().back().lhs_id == jt.lhs_id &&
                  overlaps.back().back().rhs_id == jt.rhs_id) {
                if (overlap_length(overlaps.back().back()) < overlap_length(jt)) {
                  overlaps.back().back() = jt;
                }
              } else {
                overlaps.back().emplace_back(jt);
              }
            }
          }
        }
        thread_futures.clear();

        std::cerr << "[raven::Graph::Construct] mapped valid sequences "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;

        j = i + 1;
      }

      timer.Start();

      std::vector<std::future<void>> thread_futures;
      for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
        if (graph_.piles_[i]->is_contained()) {
          graph_.piles_[i]->set_is_invalid();
        }
      }

      {
        std::uint32_t k = 0;
        for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
          if (overlap_update(overlaps.back()[i])) {
            overlaps.back()[k++] = overlaps.back()[i];
          }
        }
        overlaps.back().resize(k);
      }

      std::cerr << "[raven::Graph::Construct] updated overlaps "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      std::sort(sequences.begin(), sequences.end(),
                [&](const std::unique_ptr<biosoup::NucleicAcid> &lhs,
                    const std::unique_ptr<biosoup::NucleicAcid> &rhs) -> bool {
                  return lhs->id < rhs->id;
                });
    }

    // resolve repeat induced overlaps
    if (graph_.stage() == Graph_Stage::Construct_Graph_2) {
      timer.Start();

      while (true) {
        auto components = connected_components();
        for (const auto &it: components) {
          std::vector<std::uint16_t> medians;
          for (const auto &jt: it) {
            medians.emplace_back(graph_.piles_[jt]->median());
          }
          std::nth_element(
              medians.begin(),
              medians.begin() + medians.size() / 2,
              medians.end());
          std::uint16_t median = medians[medians.size() / 2];

          std::vector<std::future<void>> futures;
          for (const auto &jt: it) {
            futures.emplace_back(thread_pool_->Submit(
                [&](std::uint32_t i) -> void {
                  graph_.piles_[i]->FindRepetitiveRegions(median);
                },
                jt));
          }
          for (const auto &it: futures) {
            it.wait();
          }
        }

        for (const auto &it: overlaps.back()) {
          graph_.piles_[it.lhs_id]->UpdateRepetitiveRegions(it);
          graph_.piles_[it.rhs_id]->UpdateRepetitiveRegions(it);
        }

        bool is_changed = false;
        std::uint32_t j = 0;
        for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
          const auto &it = overlaps.back()[i];
          if (graph_.piles_[it.lhs_id]->CheckRepetitiveRegions(it) ||
              graph_.piles_[it.rhs_id]->CheckRepetitiveRegions(it)) {
            is_changed = true;
          } else {
            overlaps.back()[j++] = it;
          }
        }
        overlaps.back().resize(j);

        if (!is_changed) {
          break;
        }

        for (const auto &it: components) {
          for (const auto &jt: it) {
            graph_.piles_[jt]->ClearRepetitiveRegions();
          }
        }
      }

      std::cerr << "[raven::Graph::Construct] removed false overlaps "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();
    }

    // construct assembly graph
    if (graph_.stage() == Graph_Stage::Construct_Graph_2) {
      Node::num_objects = 0;
      Edge::num_objects = 0;

      std::vector<std::int32_t> sequence_to_node(graph_.piles_.size(), -1);
      for (const auto &it: graph_.piles_) {  // create nodes
        if (it->is_invalid()) {
          continue;
        }

        std::unordered_set<std::uint32_t> annotations;
        for (const auto &jt: graph_.annotations_[it->id()]) {
          if (it->begin() <= jt && jt < it->end()) {
            annotations.emplace(jt - it->begin());
          }
        }
        graph_.annotations_[it->id()].swap(annotations);

        auto sequence = biosoup::NucleicAcid{
            sequences[it->id()]->name,
            sequences[it->id()]->InflateData(it->begin(), it->end() - it->begin())};  // NOLINT
        sequence.id = it->id();

        sequence_to_node[it->id()] = Node::num_objects;

        auto node = std::make_shared<Node>(sequence);
        sequence.ReverseAndComplement();
        graph_.nodes_.emplace_back(node);
        graph_.nodes_.emplace_back(std::make_shared<Node>(sequence));
        node->pair = graph_.nodes_.back().get();
        node->pair->pair = node.get();

        if (it->id() < split) {
          node->color = 1;
          node->pair->color = 1;
        }
      }

      std::cerr << "[raven::Graph::Construct] stored " << graph_.nodes_.size() << " nodes "  // NOLINT
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();

      for (auto &it: overlaps.back()) {  // create edges
        if (!overlap_finalize(it)) {
          continue;
        }

        auto tail = graph_.nodes_[sequence_to_node[it.lhs_id]].get();
        auto head = graph_.nodes_[sequence_to_node[it.rhs_id] + 1 - it.strand].get();

        auto length = it.lhs_begin - it.rhs_begin;
        auto length_pair =
            (graph_.piles_[it.rhs_id]->length() - it.rhs_end) -
            (graph_.piles_[it.lhs_id]->length() - it.lhs_end);

        if (it.score == 4) {
          std::swap(head, tail);
          length *= -1;
          length_pair *= -1;
        }

        auto edge = std::make_shared<Edge>(tail, head, length);
        graph_.edges_.emplace_back(edge);
        graph_.edges_.emplace_back(std::make_shared<Edge>(head->pair, tail->pair, length_pair));  // NOLINT
        edge->pair = graph_.edges_.back().get();
        edge->pair->pair = edge.get();
      }

      std::cerr << "[raven::Graph::Construct] stored " << graph_.edges_.size() << " edges "  // NOLINT
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa("after_construction.gfa", false);
    }

    // checkpoint
    if (graph_.stage() == Graph_Stage::Construct_Graph_2) {
      graph_.advance_stage();
      if (graph_.use_checkpoints()) {
        timer.Start();
        graph_.Store();
        std::cerr << "[raven::Graph::Construct] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }

    std::cerr << "[raven::Graph::Construct] "
              << std::fixed << timer.elapsed_time() << "s"
              << std::endl;
  }

  void Graph_Constructor::LoadOverlaps(const std::string &overlaps_path, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, std::vector<std::vector<extended_overlap>> &extended_overlaps){
    std::ifstream file(overlaps_path);
    if (!file.is_open()) {
      throw std::runtime_error("Error opening file: " + overlaps_path);
    }
    
    auto reverse_hifiasm_overlap = [](const extended_overlap &eo) -> extended_overlap {
      return extended_overlap{
        biosoup::Overlap(
              eo.overlap.rhs_id, eo.overlap.rhs_begin, eo.overlap.rhs_end,
              eo.overlap.lhs_id, eo.overlap.lhs_begin, eo.overlap.lhs_end,
              eo.overlap.score,
              eo.overlap.strand),
        edlib_align{}
      };
    };
    
    auto get_read_id = [](const std::string &read_name, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences) -> std::uint32_t {
      for (const auto &it: sequences) {
        if (it->name == read_name) {
          return it->id;
        }
      }
      return -1;
    };

    std::cerr << "[raven::Graph::LoadHerroSNPs] loading overlaps from: " << overlaps_path << std::endl;
    std::string line;
    while(std::getline(file, line)){
      std::istringstream iss(line);
      std::string item;
      std::vector<std::string> items;
      std::uint32_t lhs_seq_id;
      std::uint32_t rhs_seq_id;

      while(std::getline(iss, item, '\t')){
        items.push_back(item);
      };

      lhs_seq_id = get_read_id(items[0], sequences);
      rhs_seq_id = get_read_id(items[5], sequences);

      if(lhs_seq_id == -1 || rhs_seq_id == -1){
        continue;
      } else {
        biosoup::Overlap overlap{lhs_seq_id, std::stoi(items[2]), std::stoi(items[3]), 
                                  rhs_seq_id, std::stoi(items[7]), std::stoi(items[8]), 
                                  255, items[4] == "+"};
        edlib_align tmp = {};
        extended_overlap total_ovlp{overlap, tmp};
        extended_overlaps[lhs_seq_id].emplace_back(total_ovlp);
        //extended_overlaps[rhs_seq_id].emplace_back(reverse_hifiasm_overlap(total_ovlp));
      }
    };
    std::cerr << "[raven::Graph::LoadHerroSNPs] loaded overlaps from: " << overlaps_path << std::endl;
  };

  void Graph_Constructor::LoadHerroSNPs(const std::string &herro_snps_path, std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences){
    std::ifstream file(herro_snps_path);
    if (!file.is_open()) {
      throw std::runtime_error("Error opening file: " + herro_snps_path);
    }

    std::cerr << "[raven::Graph::LoadHerroSNPs] loading snps from: " << herro_snps_path << std::endl;
    std::string line;
    while(std::getline(file, line)){
          std::string single_line = line;
          std::istringstream iss(single_line);

          std::string item;
          std::string first_item;
          bool is_first = true;
          std::uint32_t seq_id;
          bool found = false;


          while (std::getline(iss, item, '\t')) {
            if (is_first) {
              first_item = item;
              is_first = false;
              for(auto &read : sequences){

                if(read.get()->name == first_item){
                  seq_id = read.get()->id;
                  found = true;
                  break;
                }
              }
              continue;
              } else if(found){
                std::uint32_t pos = std::stoi(item);
                graph_.annotations_[seq_id].emplace(pos);
              }

          }
    };


  };


  void Graph_Constructor::LoadFromGfa(const std::string &gfa_path){
    try {
        std::string gfa_path_without_leading_whitespace;
        if (!gfa_path.empty()) {
            gfa_path_without_leading_whitespace = gfa_path.substr(1);
        }
        std::ifstream file(gfa_path_without_leading_whitespace);

        if (!file.is_open()) {
            throw std::runtime_error("Error opening file: " + gfa_path_without_leading_whitespace);
        }

        std::string line;
        Node::num_objects = 0;
        Edge::num_objects = 0;
        std::map<std::string, std::shared_ptr<Node>> sequence_to_node;
        
        while (std::getline(file, line)) {
          // Process each line here
          std::string single_line = line;
          std::istringstream iss(single_line);

          std::string item;
          std::string first_item;
          bool is_first = true;
          std::uint8_t counter = 0;
          std::uint32_t sequence_counter = 0;
          std::string seq_name;
          std::string nuc_sequence;

          while (std::getline(iss, item, '\t')) {
              if (counter == 0) {
                  first_item = item;
                  if(first_item == "S"){
                    sequence_counter++;
                    counter++;
                    continue;
                  }else{
                    break;
                  }
              }

              if(first_item == "S"){
                if(counter == 1){
                  seq_name = item;
                } else if(counter == 2){
                  nuc_sequence = item;
              
                }
              }
              counter++;
          }
          if(first_item == "S"){
            auto sequence = biosoup::NucleicAcid{
              seq_name,
              nuc_sequence
            };
            sequence.id = sequence_counter;

            auto node = std::make_shared<Node>(sequence);
            sequence.ReverseAndComplement();
            graph_.nodes_.emplace_back(node);
            graph_.nodes_.emplace_back(std::make_shared<Node>(sequence));
            node->pair = graph_.nodes_.back().get();
            node->pair->pair = node.get();        
            sequence_to_node.emplace(seq_name, node);
          }
        }
        if(file.eof()){
          std::cerr << "[raven::Graph::LoadFromGfa] loaded sequences from: " << gfa_path_without_leading_whitespace << std::endl; 
        }
        //file.seekg(0, file.beg);

        std::ifstream file2(gfa_path_without_leading_whitespace);
        
        while (std::getline(file2, line)) {
          // Process each line here
          std::string single_line = line;
          std::istringstream iss(single_line);

          std::string item;
          std::string first_item;
          bool is_first = true;
          std::uint8_t counter = 0;
          std::uint32_t sequence_counter = 0;
          std::string seq_name;
          std::string nuc_sequence;
          std::string tail_node_name;
          std::string head_node_name;

          bool tail_node_strand;
          bool head_node_strand;

          std::string edge_length;
          std::string item2;
          std::string ol_length;

          while (std::getline(iss, item, '\t')) {
            if (counter == 0) {
                first_item = item;
                if(first_item == "L"){
                  sequence_counter++;
                  counter++;
                  continue;
                } else {
                  break;
                }
            }
            if(first_item == "L"){
              if(counter == 1){
                tail_node_name = item;
              }else if(counter == 2){
                tail_node_strand = item == "+" ? true : false;
              }else if(counter == 3){
                head_node_name = item;
              }else if(counter == 4){
                head_node_strand = item == "+" ? true : false;
              }else if(counter == 5){
                std::stringstream ss(item);
                while(std::getline(ss, item2 , 'M')){
                  ol_length = item2;
                }
              }else if(counter == 6){
                std::uint8_t mini_counter = 0;
                std::stringstream ss(item);
                while(std::getline(ss, item2 , ':')){
                  if(mini_counter == 2) edge_length = item2;
                  mini_counter++;
                }
              }
            }
            counter++;
          }
          if(first_item == "L"){
            
            auto tail_node = tail_node_strand ? sequence_to_node[tail_node_name].get() : sequence_to_node[tail_node_name]->pair;
            auto head_node = head_node_strand ? sequence_to_node[head_node_name].get() : sequence_to_node[head_node_name]->pair;

            auto length = std::stoi(edge_length);
            auto length_pair = head_node->sequence.inflated_len - std::stoi(ol_length);

            auto edge = std::make_shared<Edge>(tail_node, head_node, length);
            graph_.edges_.emplace_back(edge);
            graph_.edges_.emplace_back(std::make_shared<Edge>(head_node->pair, tail_node->pair, length_pair));  // NOLINT
            edge->pair = graph_.edges_.back().get();
            edge->pair->pair = edge.get();
            
          }
        }
          //std::cout << line << std::endl;
        file2.close();
      if (file.eof()) {
          // File has been read successfully
          file.close();
          std::cerr << "[raven::Graph::LoadFromGfa] successfully loaded graph from: " << gfa_path_without_leading_whitespace << std::endl;
      } else {
          throw std::runtime_error("Error reading file: " + gfa_path_without_leading_whitespace);
      }
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    while(graph_.stage() != Graph_Stage::Assemble_Transitive_Edges){
      graph_.advance_stage();
    }
  };
  // NOLINT
} // raven