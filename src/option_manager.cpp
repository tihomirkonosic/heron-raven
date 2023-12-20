
#include <getopt.h>
#include <iostream>
#include "option_manager.h"

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

int ProcessParameters(int argc, char **argv, Program_Parameters& param) {
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 'K':
        param.kmer_len = std::atoi(optarg);
        break;
      case 'W':
        param.window_len = std::atoi(optarg);
        break;
      case 'F':
        param.freq = std::atof(optarg);
        break;
      case 'H': param.hpc = true; break;
      case 'B': param.bandwidth = std::atoi(optarg); break;
      case 'C': param.chain_n = std::atoi(optarg); break;
      case 'M': param.match_n = std::atoi(optarg); break;
      case 'G': param.gap_size = std::atoi(optarg); break;
      case 'E': param.error_corrected_reads = optarg; break;
      case 's':
        param.split = std::atoi(optarg);
        break;
      case 'p':
        param.num_polishing_rounds = atoi(optarg);
        break;
      case 'y': param.ploidy = std::atoi(optarg); break;
      case 'D':
        param.disagreement = std::atof(optarg);
        break;
      case 'f':
        param.gfa_path = optarg;
        break;
      case 'r':
        param.resume = true;
        break;
      case 'd':
        param.checkpoints = false;
        break;
      case 't':
        param.num_threads = atoi(optarg);
        break;
      case 'v':
        std::cout << RAVEN_VERSION << std::endl;
        return 0;
      case 'h':
        Help();
        return 0;
      case 'o':
        param.output_path = optarg;
        param.stdoutput = false;
        break;
      case 'u':
        param.ul_read_path = optarg;
        break;
      case 'x':
        param.kMaxNumOverlaps = std::atof(optarg);
        break;
      case 'U':
        param.min_unitig_size = std::atoi(optarg);
        break;
      case 'P':
        param.paf = true;
        break;
      case 'R':
        param.valid_region_size = std::atoi(optarg);
        break;
      case 'g':
        param.input_gfa_path = optarg;
        param.skip_contruction = true;
        break;
      default:
        return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  return 1;
}
