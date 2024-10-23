
#include <cmath>
#include <deque>
#include <fstream>
#include "graph.hpp"
#include "graph_constructor.h"
#include "overlap.h"
#include "overlap_parser.h"
#include "biosoup/overlap.hpp"
#include "edlib.h"
#include "biosoup/timer.hpp"
#include "extended_overlap.h"
#include "overlap_helpers.h"
#include "annotation_helpers.h"

namespace raven {
  Graph_Constructor::Graph_Constructor(Graph &graph, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
      : graph_(graph), thread_pool_(thread_pool ?
                                    thread_pool :
                                    std::make_shared<thread_pool::ThreadPool>(1)) {
  }


  void Graph_Constructor::Construct(
      std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,  // NOLINT
      Program_Parameters &param) {

    disagreement_ = param.disagreement;

    if (sequences.empty()) {
      return;
    }

    std::vector<std::vector<biosoup::Overlap>> overlaps(sequences.size());
    std::vector<std::vector<extended_overlap>> extended_overlaps(sequences.size());

    graph_.annotations_.resize(sequences.size());

    // checkpoint test
    if (graph_.stage() == Graph_Stage::Construct_Graph && graph_.use_checkpoints()) {
      graph_.Store(param.cereal_filename);
    }

    biosoup::Timer timer{};

    ram::MinimizerEngine minimizer_engine{
        thread_pool_,
        param.kmer_len,
        param.window_len,
        param.bandwidth,
        param.chain_n,
        param.match_n,
        param.gap_size
    };

    // find overlaps and create piles
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      for (const auto &it: sequences) {
        graph_.piles_.emplace_back(new Pile(it->id, it->inflated_len));
      }
      
      if(!param.load_paf.empty()) {
        LoadOverlaps(param.load_paf, sequences, extended_overlaps);
        std::vector<std::future<void>> void_futures;
        for(int i = 0; i < sequences.size(); i++){
          void_futures.emplace_back(thread_pool_->Submit(
              [&](std::uint32_t i) -> void {
                find_pairwise_alignment(i, extended_overlaps[i], sequences);
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
        minimizer_engine.Filter(param.freq);

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
                                                                           true, param.hpc);
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
                        ovlp.lhs_begin = ovlp.lhs_begin - (param.window_len+param.kmer_len-1) ? ovlp.lhs_begin - (param.window_len+param.kmer_len-1) : 0;
                        ovlp.lhs_end = ovlp.lhs_end + (param.window_len+param.kmer_len-1) < sequences[ovlp.lhs_id]->inflated_len ? ovlp.lhs_end + (param.window_len+param.kmer_len-1) : sequences[ovlp.lhs_id]->inflated_len;

                        ovlp.rhs_begin = ovlp.rhs_begin - (param.window_len+param.kmer_len-1) ? ovlp.rhs_begin - (param.window_len+param.kmer_len-1) : 0;
                        ovlp.rhs_end = ovlp.rhs_end + (param.window_len+param.kmer_len-1) < sequences[ovlp.rhs_id]->inflated_len ? ovlp.rhs_end + (param.window_len+param.kmer_len-1) : sequences[ovlp.rhs_id]->inflated_len;

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

    if(param.herro_snps_path == ""){
      std::vector<std::future<void>> void_futures;
      for(int i = 0; i < sequences.size(); i++){
        void_futures.emplace_back(thread_pool_->Submit(
            [&](std::uint32_t i) -> void {
              call_snps(i, extended_overlaps[i], sequences, graph_);
            },
            i));
      };

      for (const auto &it: void_futures) {
        it.wait();
      }
    
      std::cerr << "[raven::Graph::Construct] snps called"
                << std::endl;

      if (print_snp_data && param.ploidy >= 2) {
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
      LoadHerroSNPs(param.herro_snps_path, sequences);
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

    graph_.PrintOverlaps(extended_overlaps, sequences, true, "beforeValidRegion.paf");

    // trim and annotate piles
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      timer.Start();
      std::vector<std::future<void>> thread_futures;
      for (std::uint32_t i = 0; i < graph_.piles_.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->Submit(
            [&](std::uint32_t i) -> void {
              graph_.piles_[i]->FindValidRegion(param.valid_region_coverage_threshold, param.valid_region_length_threshold);
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

    // resolve contained reads
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      timer.Start();
      std::vector<std::future<void>> futures;
      for (std::uint32_t i = 0; i < extended_overlaps.size(); ++i) {
        futures.emplace_back(thread_pool_->Submit(
            [&](std::uint32_t i) -> void {
              std::uint32_t k = 0;
              for (std::uint32_t j = 0; j < extended_overlaps[i].size(); ++j) {
                if (!overlap_update(extended_overlaps[i][j].overlap, graph_)) {
                  continue;
                }

                const auto &it = extended_overlaps[i][j];

                auto lhs_anno = annotation_extract(
                    it.overlap.lhs_id,
                    it.overlap.lhs_begin,
                    it.overlap.lhs_end,
                    sequences[it.overlap.lhs_id]->inflated_len,
                    true,
                    graph_);

                auto rhs_anno = annotation_extract(
                    it.overlap.rhs_id,
                    it.overlap.rhs_begin,
                    it.overlap.rhs_end,
                    sequences[it.overlap.rhs_id]->inflated_len,
                    it.overlap.strand,
                    graph_);

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
          if (!overlap_update(overlaps[i][j], graph_)) {
            continue;
          }
          std::uint32_t type = overlap_type(overlaps[i][j], graph_);
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

    graph_.PrintOverlaps(extended_overlaps, sequences, true, param.paf_after_contained_filename);

    // resolve chimeric sequences
    if (graph_.stage() == Graph_Stage::Construct_Graph) {
      timer.Start();

      while (true) {
        auto components = connected_components(sequences, overlaps, graph_);
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
            if (overlap_update(overlaps[i][j], graph_)) {
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
              std::uint32_t type = overlap_type(jt, graph_);
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
        graph_.Store(param.cereal_filename);
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
                    param.hpc,
                    &filtered);
                graph_.piles_[sequences[i]->id]->AddKmers(filtered, param.kmer_len, sequences[i]); // NOLINT

                std::uint32_t k = 0;
                for (std::uint32_t j = 0; j < dst.size(); ++j) {
                  if (!overlap_update(dst[j], graph_) || overlap_length(dst[j]) < 1000) {
                    continue;
                  }

                  const auto &it = dst[j];

                  auto lhs_anno = annotation_extract(
                      it.lhs_id,
                      it.lhs_begin,
                      it.lhs_end,
                      sequences[sequences_map[it.lhs_id]]->inflated_len,
                      true,
                      graph_);

                  auto rhs_anno = annotation_extract(
                      it.rhs_id,
                      it.rhs_begin,
                      it.rhs_end,
                      sequences[sequences_map[it.rhs_id]]->inflated_len,
                      it.strand,
                      graph_);

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
            if (!overlap_update(jt, graph_)) {
              continue;
            }
            std::uint32_t type = overlap_type(jt, graph_);
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
          if (overlap_update(overlaps.back()[i], graph_)) {
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
        auto components = connected_components(sequences, overlaps, graph_);
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

        if (it->id() < param.split) {
          node->color = 1;
          node->pair->color = 1;
        }
      }

      std::cerr << "[raven::Graph::Construct] stored " << graph_.nodes_.size() << " nodes "  // NOLINT
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();

      for (auto &it: overlaps.back()) {  // create edges
        if (!overlap_finalize(it, graph_)) {
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

      graph_.PrintGfa(param.gfa_after_construction_filename, false);
    }

    // checkpoint
    if (graph_.stage() == Graph_Stage::Construct_Graph_2) {
      graph_.advance_stage();
      if (graph_.use_checkpoints()) {
        timer.Start();
        graph_.Store(param.cereal_filename);
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
        biosoup::Overlap overlap{lhs_seq_id, (std::uint32_t)std::stoi(items[2]), (std::uint32_t)std::stoi(items[3]),
                                  rhs_seq_id, (std::uint32_t)std::stoi(items[7]), (std::uint32_t)std::stoi(items[8]),
                                  255, items[4] == "+"};
        edlib_align tmp = {};
        extended_overlap total_ovlp{overlap, tmp};
        extended_overlaps[lhs_seq_id].emplace_back(total_ovlp);
        //extended_overlaps[rhs_seq_id].emplace_back(reverse_hifiasm_overlap(total_ovlp));
      }
    }
    std::cerr << "[raven::Graph::LoadHerroSNPs] loaded overlaps from: " << overlaps_path << std::endl;
  }

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
    }


  }


  void Graph_Constructor::LoadFromGfa(const std::string &gfa_path){
    try {
//        std::string gfa_path_without_leading_whitespace;
//        if (!gfa_path.empty()) {
//            gfa_path_without_leading_whitespace = gfa_path.substr(1);
//        }
        std::ifstream file(gfa_path);

        if (!file.is_open()) {
            throw std::runtime_error("Error opening file: " + gfa_path);
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
          std::cerr << "[raven::Graph::LoadFromGfa] loaded sequences from: " << gfa_path << std::endl;
        }
        //file.seekg(0, file.beg);

        std::ifstream file2(gfa_path);
        
        while (std::getline(file2, line)) {
          // Process each line here
          std::string single_line = line;
          std::istringstream iss(single_line);

          std::string item;
          std::string first_item;
          std::uint8_t counter = 0;
          std::uint32_t sequence_counter = 0;
          std::string seq_name;
          std::string nuc_sequence;
          std::string tail_node_name;
          std::string head_node_name;

          bool tail_node_strand;
          bool head_node_strand;

          //std::string edge_length;
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
              }else if(counter == 5) {
                std::stringstream ss(item);
                while (std::getline(ss, item2, 'M')) {
                  ol_length = item2;
                }
              }
//              }else if(counter == 6){
//                std::uint8_t mini_counter = 0;
//                std::stringstream ss(item);
//                while(std::getline(ss, item2 , ':')){
//                  if(mini_counter == 2) edge_length = item2;
//                  mini_counter++;
//                }
//              }
            }
            counter++;
          }
          if(first_item == "L"){
            
            auto tail_node = tail_node_strand ? sequence_to_node[tail_node_name].get() : sequence_to_node[tail_node_name]->pair;
            auto head_node = head_node_strand ? sequence_to_node[head_node_name].get() : sequence_to_node[head_node_name]->pair;

            //auto length = std::stoi(edge_length);
            auto length = head_node->sequence.inflated_len - std::stoi(ol_length);

            auto edge = std::make_shared<Edge>(tail_node, head_node, length);
            graph_.edges_.emplace_back(edge);
            graph_.edges_.emplace_back(std::make_shared<Edge>(head_node->pair, tail_node->pair, length));  // NOLINT
            edge->pair = graph_.edges_.back().get();
            edge->pair->pair = edge.get();
            
          }
        }
          //std::cout << line << std::endl;
        file2.close();
      if (file.eof()) {
          // File has been read successfully
          file.close();
          std::cerr << "[raven::Graph::LoadFromGfa] successfully loaded graph from: " << gfa_path << std::endl;
      } else {
          throw std::runtime_error("Error reading file: " + gfa_path);
      }
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    while(graph_.stage() != Graph_Stage::Assemble_Transitive_Edges){
      graph_.advance_stage();
    }
  }

  void Graph_Constructor::LoadFromPaf(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences, const std::string &paf_path){
    try {
      gzFile file = gzopen(paf_path.c_str(), "r");
      OverlapParser parser {file};
      std::vector<std::unique_ptr<Overlap>> overlaps = parser.ParsePaf((std::uint64_t )-1);

      Node::num_objects = 0;
      Edge::num_objects = 0;
      std::map<std::string, std::shared_ptr<Node>> sequence_to_node;
      std::uint32_t sequence_counter = 0;

      for (std::unique_ptr<biosoup::NucleicAcid> &seq_ptr : sequences) {
        biosoup::NucleicAcid seq = *seq_ptr.get();
        seq.id = ++sequence_counter;
        std::shared_ptr<Node> node = std::make_shared<Node>(seq);
        graph_.nodes_.emplace_back(node);
        seq.ReverseAndComplement();
        graph_.nodes_.emplace_back(std::make_shared<Node>(seq));
        node->pair = graph_.nodes_.back().get();
        node->pair->pair = node.get();
        sequence_to_node.emplace(seq.name, node);
      }

      bool tail_node_strand;
      bool head_node_strand;

      for (std::unique_ptr<Overlap> &overlap_ptr : overlaps) {
        Overlap overlap = *overlap_ptr.get();
        tail_node_strand = true;
        head_node_strand = overlap.strand;

        Node *tail_node = tail_node_strand ? sequence_to_node[overlap.q_name].get() : sequence_to_node[overlap.q_name]->pair;
        Node * head_node = head_node_strand ? sequence_to_node[overlap.t_name].get() : sequence_to_node[overlap.t_name]->pair;

        uint32_t length = overlap.q_len - overlap.overlap_len;
        uint32_t length_pair = overlap.t_len - overlap.overlap_len;

        std::shared_ptr<Edge> edge = std::make_shared<Edge>(tail_node, head_node, length);
        graph_.edges_.emplace_back(edge);
        graph_.edges_.emplace_back(std::make_shared<Edge>(head_node->pair, tail_node->pair, length_pair));  // NOLINT
        edge->pair = graph_.edges_.back().get();
        edge->pair->pair = edge.get();
      }

    } catch (const std::exception& e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }

    while(graph_.stage() != Graph_Stage::Assemble_Transitive_Edges){
      graph_.advance_stage();
    }
  }
  // NOLINT
} // raven