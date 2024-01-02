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
      graph.Load();
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

  raven::Graph_Constructor graph_constructor{graph, thread_pool};
  if (!param.skip_contruction){
    std::cout << "Constructing graph with params: kmer_size:" << param.kmer_len  << " winodw_size:" << param.window_len << " " << std::endl;
    graph_constructor.Construct(sequences, param.disagreement, param.split, param.max_overlaps, param.ploidy, param.kmer_len, param.window_len, param.bandwidth, param.chain_n, param.match_n, param.gap_size, param.freq, param.hpc, param.paf, param.valid_region_size);
  } else {
    graph_constructor.LoadFromGfa(param.input_gfa_path);
  }

  graph.PrintGfa("post_construction.gfa");
  raven::Graph_Assembler graph_assembler{graph, thread_pool};
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
    graph.PrintGfa("post_cleaning.gfa");
    if(param.ploidy >= 2){
      graph_assembler.AssembleDiploids();
    } else {
      graph_assembler.AssembleHaploids();
    }
  } else {
    timer.Start();
    graph_assembler.UlAssemble(ul_sequences);
  }

  graph.PrintGfa(param.gfa_path);

  if (param.stdoutput) {
    // output to stdout
    for (const auto &it: graph.GetUnitigs(param.num_polishing_rounds > 0)) {
      std::cout << ">" << it->name << std::endl;
      std::cout << it->InflateData() << std::endl;
    }
  } else if(param.ploidy >= 2){
    // output to file
    std::filesystem::path root_path(param.output_path);
    std::filesystem::path noext("");
    root_path.replace_extension(noext);
    std::filesystem::path path1 = root_path;
    std::filesystem::path path2 = root_path;
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
    std::filesystem::path root_path(param.output_path);
    std::filesystem::path noext("");
    root_path.replace_extension(noext);
    std::filesystem::path path1 = root_path;
    std::ofstream outfile1;

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
