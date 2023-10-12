// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <iostream>
#include <fstream>
#include <filesystem>

#include "biosoup/timer.hpp"

#include "graph.hpp"
#include "graph_constructor.h"
#include "graph_assembler.h"
#include "parser.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

  static struct option options[] = {
      {"kmer-len", required_argument, nullptr, 'K'},
      {"window-len", required_argument, nullptr, 'W'},
      {"frequency", required_argument, nullptr, 'F'},
			{"hpc", no_argument, nullptr, 'H'},
			{"bandwidth", required_argument, nullptr, 'B'},
			{"chain_n", required_argument, nullptr, 'C'},
			{"match_n", required_argument, nullptr, 'M'},
			{"gap_size", required_argument, nullptr, 'G'},
			{"error-corrected-reads", required_argument, nullptr, 'E'},
      {"polishing-rounds", required_argument, nullptr, 'p'},
      {"split", required_argument, nullptr, 's'},
      {"disagreement", required_argument, nullptr, 'D'},
      {"graphical-fragment-assembly", required_argument, nullptr, 'f'},
      {"resume", no_argument, nullptr, 'r'},
      {"disable-checkpoints", no_argument, nullptr, 'd'},
      {"threads", required_argument, nullptr, 't'},
      {"version", no_argument, nullptr, 'v'},
      {"help", no_argument, nullptr, 'h'},
      {"output", required_argument, nullptr, 'o'},
      {"ultralong-phasing", required_argument, nullptr, 'u'},
      {"kMaxNumOverlaps", required_argument, nullptr, 'x'},
      {"ploidy", required_argument, nullptr, 'y'},
      {"min-unitig-size", required_argument, nullptr, 'U'},
      {"paf", no_argument, nullptr, 'P'},
      {"valid-region", required_argument, nullptr, 'R'},
      {"input", required_argument, nullptr, 'g'},
      {nullptr, 0, nullptr, 0}
  };

  std::string optstr = "K:W:F:Hp:s:D:f:rdt:vho:u:x:y:U:PR:g:";

  void Help() {
    std::cout <<
              "usage: raven [options ...] <sequences>\n"
              "\n"
              "  # default output is to stdout in FASTA format\n"
              "  <sequences>\n"
              "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
              "\n"
              "  options:\n"
              "    -K, --kmer-len <int>\n"
              "      default: 15\n"
              "      length of minimizers used to find overlaps\n"
              "    -W, --window-len <int>\n"
              "      default: 5\n"
              "      length of sliding window from which minimizers are sampled\n"
              "    -F, --frequency <double>\n"
              "      default: 0.001\n"
              "      threshold for ignoring most frequent minimizers\n"
              "    -p, --polishing-rounds <int>\n"
              "      default: 0\n"
              "      number of times racon is invoked\n"
							"    -H, --hpc\n"
							"     default: false\n"
							"     use HPC overlaps\n"
							"    -B, --bandwidth <int>\n"
							"      default: 500\n"
							"      bandwidth\n"
							"    -C, --chain_n <int>\n"
							"      default: 4\n"
							"      chain_n\n"
							"    -M, --match_n <int>\n"
							"      default: 100\n"
							"      match_n\n"
							"    -G, --gap_size <int>\n"
							"      default: 10000\n"
							"      gap_size\n"
							"    -E, --error-corrected-reads <string>\n"
							"      default: \"\"\n"
							"      path to error corrected reads\n"
              "    -s, --split <int>\n"
              "      default: 0\n"
              "      graph coloring\n"
              "    -D, --disagreement <double>\n"
              "      default: 0.1\n"
              "      maximal percentage of different anntoated bases in overlaps\n"
              "    -f, --graphical-fragment-assembly <string>\n"
              "      prints the assembly graph in GFA format\n"
              "    -r, --resume\n"
              "      resume previous run from last checkpoint\n"
              "    -d, --disable-checkpoints\n"
              "      disable checkpoint file creation\n"
              "    -t, --threads <int>\n"
              "      default: 1\n"
              "      number of threads\n"
              "    -v, --version\n"
              "      prints the version number\n"
              "    -o, --output <string>\n"
              "      output file name, if it is not set output is written to stdout\n"
              "      for diploid assembly, outputs will be written in 2 files with suffixes -1, -2\n"
              "    -u, --ultralong-phasing <string>\n"
              "       path to ul reads used for phasing\n"
              "    -x, --kMaxNumOverlaps <long unsigned int>\n"
              "      default: 32\n"
              "      maximum number of overlaps that will be taken during find overlaps and create piles stage\n"              "    -h, --help\n"
              "    -U, --min-unitig-size <int>\n"
              "      minimal uniting size (default 9999)\n"
              "    -P, --paf\n"
              "      overlaps are stored to files in PAF format\n"
              "    -R --valid-region <int>\n"
              "      default: 4\n"
              "      overlap valid region size\n"
              "      prints the usage\n"
              "    -g, --graph-input <string>\n"
              "      input graph path\n";
  }

}  // namespace

int main(int argc, char **argv) {
  unsigned split = 0;

	std::uint8_t ploidy = 2;

  std::uint8_t kmer_len = 15;
  std::uint8_t window_len = 5;
  std::uint16_t bandwidth = 500;
  std::uint16_t chain_n = 4;
  std::uint16_t match_n = 100;
  std::uint16_t gap_size = 10000;
  double freq = 0.001;
	bool hpc = false;

	std::string error_corrected_reads = "";

  std::int32_t num_polishing_rounds = 0;

  std::string ul_read_path;

  double disagreement = 0.1;
  std::string gfa_path = "";
  std::string input_gfa_path = "";
  bool skip_contruction = false;
  bool resume = false;
  bool checkpoints = true;

  bool stdoutput = true;
  std::string output_path = "";

  std::uint32_t num_threads = 1;

  std::size_t kMaxNumOverlaps = 16;
  std::uint32_t min_unitig_size = 9999;
  std::uint16_t valid_region_size = 4;
  bool paf = false;

  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 'K':
        kmer_len = std::atoi(optarg);
        break;
      case 'W':
        window_len = std::atoi(optarg);
        break;
      case 'F':
        freq = std::atof(optarg);
        break;
			case 'H': hpc = true; break;
			case 'B': bandwidth = std::atoi(optarg); break;
			case 'C': chain_n = std::atoi(optarg); break;
			case 'M': match_n = std::atoi(optarg); break;
			case 'G': gap_size = std::atoi(optarg); break;
			case 'E': error_corrected_reads = optarg; break;
      case 's':
        split = std::atoi(optarg);
        break;
      case 'p':
        num_polishing_rounds = atoi(optarg);
        break;
			case 'y': ploidy = std::atoi(optarg); break;
      case 'D':
        disagreement = std::atof(optarg);
        break;
      case 'f':
        gfa_path = optarg;
        break;
      case 'r':
        resume = true;
        break;
      case 'd':
        checkpoints = false;
        break;
      case 't':
        num_threads = atoi(optarg);
        break;
      case 'v':
        std::cout << VERSION << std::endl;
        return 0;
      case 'h':
        Help();
        return 0;
      case 'o':
        output_path = optarg;
        stdoutput = false;
        break;
      case 'u':
        ul_read_path = optarg;
        break;
      case 'x':
        kMaxNumOverlaps = std::atof(optarg);
        break;
      case 'U':
        min_unitig_size = std::atoi(optarg);
        break;
      case 'P':
        paf = true;
        break;
      case 'R':
        valid_region_size = std::atoi(optarg);
        break;
      case 'g':
        input_gfa_path = optarg;
        skip_contruction = true;
        break;
      default:
        return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind >= argc && !skip_contruction) {
    std::cerr << "[raven::] error: missing input file!" << std::endl;
    return 1;
  }

  raven::min_unitig_size = min_unitig_size;

  biosoup::Timer timer{};
  timer.Start();

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

  raven::Graph graph{checkpoints};




  if (resume) {
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
      || num_polishing_rounds > 2) && !skip_contruction) {
    try {
      auto sparser = parser.CreateParser(argv[optind]);
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
  if (!skip_contruction){
      graph_constructor.Construct(sequences, disagreement, split, kMaxNumOverlaps, ploidy, kmer_len, window_len, bandwidth, chain_n, match_n, gap_size, freq, hpc, paf, valid_region_size);
  } else {
      graph_constructor.LoadFromGfa(input_gfa_path);
  }

  graph.PrintGfa("post_construction.gfa");
  raven::Graph_Assembler graph_assembler{graph, thread_pool};
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> ul_sequences;
  if (!ul_read_path.empty()){
    try {
      auto ul_sequence_parser = parser.CreateParser(ul_read_path);
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
    if(ploidy >= 2){
      graph_assembler.AssembleDiploids();
    }
  } else {
    timer.Start();
    graph_assembler.UlAssemble(ul_sequences);
  }

  graph.PrintGfa(gfa_path);

  if (stdoutput) {
    // output to stdout
    for (const auto &it: graph.GetUnitigs(num_polishing_rounds > 0)) {
      std::cout << ">" << it->name << std::endl;
      std::cout << it->InflateData() << std::endl;
    }
  } else if(ploidy >= 2){
    // output to file
    std::filesystem::path root_path(output_path);
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
  else if(ploidy == 1){
    std::filesystem::path root_path(output_path);
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
