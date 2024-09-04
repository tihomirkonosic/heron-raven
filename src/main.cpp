// Copyright (c) 2020 Robert Vaser

#include <iostream>
#include <fstream>
#include <filesystem>

#include "biosoup/timer.hpp"

#include "option_manager.h"
#include "graph.hpp"
#include "graph_constructor.h"
#include "graph_assembler.h"
#include "parser.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};


int main(int argc, char **argv) {

  Program_Parameters param{};

  if (!ProcessParameters(argc, argv, param))
    return 0;



  raven::min_unitig_size = param.min_unitig_size;

  biosoup::Timer timer{};
  timer.Start();

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(param.num_threads);

  raven::Graph graph{param.checkpoints};

  if (param.resume) {
    try {
      graph.Load(param.cereal_filename);
    } catch (std::exception &exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    std::cerr << "[raven::] loaded previous run "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();
  }

  raven::Parser parser;
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;

  if ((graph.stage() == raven::Graph_Stage::Construct_Graph
      || graph.stage() == raven::Graph_Stage::Construct_Graph_2
      || param.num_polishing_rounds > 2) && !param.skip_contruction) {
    try {
      auto sparser = parser.CreateParser(param.sequence_path);
      if (sparser == nullptr) {
        return 1;
      }
      sequences = sparser->Parse(-1);
    } catch (const std::invalid_argument &exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    if (sequences.empty()) {
      std::cerr << "[raven::] error: empty sequences set" << std::endl;
      return 1;
    }

    std::cerr << "[raven::] loaded " << sequences.size() << " sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;
    } else {
      std::cerr << "[raven::] skipped sequence loading" << std::endl;
    }
    timer.Start();

    std::vector<biosoup::NucleicAcid> sequences2{sequences.size()};
    auto compress_homopolymers = [](biosoup::NucleicAcid &read){
      std::string sequence = read.InflateData();
      std::string new_sequence = "";
      std::size_t seq_ptr = 0;
      seq_ptr += 1;
      std::size_t diamer_ptr = 0;
      diamer_ptr += 3;
      std::uint32_t repeat_count = 1;
      std::uint32_t direpeat_count = 1;
      std::vector<std::pair<std::size_t, std::size_t>> homopolymer_positions;
      std::vector<std::pair<std::size_t, std::size_t>> diamer_positions;
      std::size_t current_index = 0;
      std::string tmp = "";

      for(diamer_ptr; diamer_ptr<sequence.size();){
          tmp = sequence[diamer_ptr];
          if(sequence[diamer_ptr] == sequence[diamer_ptr-2] && sequence[diamer_ptr-1] == sequence[diamer_ptr-3] && !(sequence[diamer_ptr] == sequence[diamer_ptr-1])){
              repeat_count++;
              diamer_ptr+=2;
          }else{
              if(repeat_count > 1){
                  diamer_positions.emplace_back(diamer_ptr-(2*repeat_count) - 1, diamer_ptr-2);
              }
              repeat_count = 1;
              diamer_ptr+=1;
          }
      }

      if(repeat_count > 1){
          diamer_positions.emplace_back(diamer_ptr-(2*repeat_count) - 1, diamer_ptr-2);
      }

      std::vector<std::pair<std::size_t, std::size_t>> new_diamer_pos;
      for(auto &positions : diamer_positions){
          new_sequence.append(sequence.substr(current_index, positions.first - current_index + 1));
          new_diamer_pos.emplace_back(new_sequence.size() - 1, new_sequence.size());
          current_index = positions.second;
      }
      if(current_index < sequence.size()){
        new_sequence.append(sequence.substr(current_index));
      }

      current_index = 0;
      sequence = "";
      std::size_t diamer_index = 0;
      repeat_count = 1;
      std::string tmp_2 = "";
      for(seq_ptr; seq_ptr < new_sequence.size();){
        tmp_2 = new_sequence[seq_ptr];
          if(new_sequence[seq_ptr] == new_sequence[seq_ptr-1] && (seq_ptr-1 != new_diamer_pos[diamer_index].second) && seq_ptr != new_diamer_pos[diamer_index].first){
              repeat_count++;
              seq_ptr++;
          }else if (seq_ptr == new_diamer_pos[diamer_index].first){
                if(diamer_index < new_diamer_pos.size()-1) diamer_index++;
                if(repeat_count > 1) homopolymer_positions.emplace_back(seq_ptr-repeat_count, seq_ptr-1);
                repeat_count = 1;
                seq_ptr += 3;
               
          }else if(seq_ptr == new_diamer_pos[diamer_index].second){
                if(diamer_index < new_diamer_pos.size()-1) diamer_index++;
                seq_ptr+=2;
                repeat_count = 1;
          }else {
            if(repeat_count > 1){
                  homopolymer_positions.emplace_back(seq_ptr-repeat_count, seq_ptr-1);
            }
            repeat_count = 1;
            seq_ptr++;
          }
      }
      
      if(repeat_count > 1){
          homopolymer_positions.emplace_back(seq_ptr-repeat_count+1, seq_ptr-1);
      }
      
      std::vector<std::pair<std::size_t, std::size_t>> new_homopolymer_pos;
      std::vector<std::size_t> diamer_pos_final;
      current_index = 0;
      for(auto &positions : homopolymer_positions){
          sequence.append(new_sequence.substr(current_index, positions.first - current_index + 1));
          new_homopolymer_pos.emplace_back(sequence.size() - 1, sequence.size());
          current_index = positions.second + 1;
      }
      if(current_index < new_sequence.size()){
        sequence.append(new_sequence.substr(current_index));
      }
      return biosoup::NucleicAcid(read.name, sequence);
    };

    for(std::uint32_t i = 0; i < sequences.size(); i++){
      sequences2[i] = compress_homopolymers(*sequences[i]);
    }


  // raven::Extended_reads compressed_reads{sequences, thread_pool};

  std::ofstream os("compressed_reads.fasta");
  for(auto &it: sequences2){
    os << ">" << it.name << std::endl;
    os << it.InflateData() << std::endl;
  }
  exit(0);

  raven::Graph_Constructor graph_constructor{graph, thread_pool};
  if (!param.skip_contruction){
    std::cout << "Constructing graph with params: kmer_size:" << param.kmer_len  << " winodw_size:" << param.window_len << " " << std::endl;
    graph_constructor.Construct(sequences, param);
  } else if (param.load_input_gfa) {
    graph_constructor.LoadFromGfa(param.input_gfa_path);
  } else if (param.load_input_paf) {
    graph_constructor.LoadFromPaf(sequences, param.input_paf_path);
  } else {
    std::cerr << "[raven::] error: unknown option" << std::endl;
    return 1;
  }


  graph.PrintGfa(param.gfa_post_construction_filename, param.print_gfa_seq);
  raven::Graph_Assembler graph_assembler{graph, param, thread_pool};
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> ul_sequences;
  if (!param.ul_read_path.empty()){
    try {
      auto ul_sequence_parser = parser.CreateParser(param.ul_read_path);
      ul_sequences = ul_sequence_parser->Parse(-1);
    } catch (const std::invalid_argument &exception) {
      std::cerr << exception.what() << std::endl;
    }

    if (ul_sequences.empty()) {
      std::cerr << "[raven::] error: ul read path set but the file appears empty\n" 
                << "[raven::] reverting to no ul assembly"<< std::endl;
    } else{
      std::cerr << "[raven::] loaded " << ul_sequences.size() << " ul sequences "
                << std::fixed << timer.Stop() << "s"
                << std::endl;
    }
  }

  if (ul_sequences.empty()) {
    graph_assembler.Assemble();
    graph.PrintGfa(param.gfa_post_cleaning_filename, param.print_gfa_seq);
    if(param.ploidy >= 2){
      graph_assembler.AssembleDiploids();
    } else {
      graph_assembler.AssembleHaploids();
    }
  } else {
    timer.Start();
    graph_assembler.UlAssemble(ul_sequences);
  }

  graph.PrintGfa(param.gfa_path, param.print_gfa_seq);

  if(param.ploidy >= 2) {
    // output to file
    std::filesystem::path path1 = param.root_path;
    std::filesystem::path path2 = param.root_path;
    path1 += "-1.fasta";
    path2 += "-2.fasta";
    std::ofstream outfile1, outfile2;
    outfile1.open(path1);
    if (!outfile1.is_open()) {
      std::cerr << "[raven::] error: cannot open file" << path1 << std::endl;
      return 1;
    }
    outfile2.open(path2);
    if (!outfile2.is_open()) {
      std::cerr << "[raven::] error: cannot open file" << path2 << std::endl;
      outfile1.close();
      return 1;
    }

    for (const auto &it: graph.GetAssembledData(true)) {
      outfile1 << ">" << it->name << std::endl;
      outfile1 << it->InflateData() << std::endl;
    }

    for (const auto &it: graph.GetAssembledData(false)) {
      outfile2 << ">" << it->name << std::endl;
      outfile2 << it->InflateData() << std::endl;
    }

    outfile1.close();
    outfile2.close();
  }
  else if(param.ploidy == 1){
    std::filesystem::path path1 = param.root_path;
    path1 += ".fasta";
    std::ofstream outfile1;

    outfile1.open(path1);
    if (!outfile1.is_open()) {
      std::cerr << "[raven::] error: cannot open file" << path1 << std::endl;
      return 1;
    }

    for (const auto &it: graph.GetAssembledData(true)) {
      outfile1 << ">" << it->name << std::endl;
      outfile1 << it->InflateData() << std::endl;
    }

    outfile1.close();
  }

  timer.Stop();
  std::cerr << "[raven::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return 0;
}
