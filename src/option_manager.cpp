
#include <getopt.h>
#include <iostream>
#include "option_manager.h"

enum program_opt {
  opt_error_corrected_reads,
  opt_split,
  opt_graphical_fragment_assembly,
  opt_resume,
  opt_disable_checkpoints,
  opt_version,
  opt_help,
  opt_output,
  opt_input_gfa,
  opt_paf,

  opt_kmer_len,
  opt_window_len,
  opt_frequency,
  opt_hpc,
  opt_bandwidth,
  opt_chain_n,
  opt_match_n,
  opt_gap_size,
  opt_threads,

  opt_max_overlaps,
  opt_disagreement,
  opt_valid_region,
  opt_cleaning_rounds,
  opt_min_unitig_size,
  opt_ploidy,

  opt_ultralong_phasing,
  opt_polishing_rounds
};

static struct option options[] = {
    {"error-corrected-reads", required_argument, nullptr, opt_error_corrected_reads},
    {"split", required_argument, nullptr, opt_split},
    {"graphical-fragment-assembly", required_argument, nullptr, opt_graphical_fragment_assembly},
    {"resume", no_argument, nullptr, opt_resume},
    {"disable-checkpoints", no_argument, nullptr, opt_disable_checkpoints},
    {"version", no_argument, nullptr, opt_version},
    {"help", no_argument, nullptr, opt_help},
    {"output", required_argument, nullptr, opt_output},
    {"input-gfa", required_argument, nullptr, opt_input_gfa},
    {"paf", no_argument, nullptr, opt_paf},

    {"kmer-len", required_argument, nullptr, opt_kmer_len},
    {"window-len", required_argument, nullptr, opt_window_len},
    {"frequency", required_argument, nullptr, opt_frequency},
    {"hpc", no_argument, nullptr, opt_hpc},
    {"bandwidth", required_argument, nullptr, opt_bandwidth},
    {"chain_n", required_argument, nullptr, opt_chain_n},
    {"match_n", required_argument, nullptr, opt_match_n},
    {"gap_size", required_argument, nullptr, opt_gap_size},
    {"threads", required_argument, nullptr, opt_threads},

    {"max-overlaps", required_argument, nullptr, opt_max_overlaps},
    {"disagreement", required_argument, nullptr, opt_disagreement},
    {"valid-region", required_argument, nullptr, opt_valid_region},
    {"cleaning-rounds", required_argument, nullptr, opt_cleaning_rounds},
    {"min-unitig-size", required_argument, nullptr, opt_min_unitig_size},
    {"ploidy", required_argument, nullptr, opt_ploidy},

    {"ultralong-phasing", required_argument, nullptr, opt_ultralong_phasing},
    {"polishing-rounds", required_argument, nullptr, opt_polishing_rounds},

    {nullptr, 0, nullptr, 0}
};

std::string optstr = "K:W:F:HB:C:M:G:t:X:D:R:";

void Help() {
  std::cout <<
            "Usage: raven [options ...] <sequences>\n"
            "\n"
            "  # default output is to stdout in FASTA format\n"
            "  <sequences>\n"
            "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
            "\n"
            "Options:\n"
            "  Input/Output:\n"
            "    --error-corrected-reads <string>\n"
            "      default: \"\"\n"
            "      path to error corrected reads\n"
            "    --split <int>\n"
            "      default: 0\n"
            "      graph coloring\n"
            "    --graphical-fragment-assembly <string>\n"
            "      prints the assembly graph in GFA format\n"
            "    --resume\n"
            "      resume previous run from last checkpoint\n"
            "    --disable-checkpoints\n"
            "      disable checkpoint file creation\n"
            "    --version\n"
            "      prints the version number\n"
            "    --help\n"
            "      prints this help text\n"
            "    --output <string>\n"
            "      output file name, if it is not set output is written to stdout\n"
            "      for diploid assembly, outputs will be written in 2 files with suffixes -1, -2\n"
            "    --input-gfa <string>\n"
            "      input GFA file name, if it is set raven will skip construction phase\n"
            "    --paf\n"
            "      overlaps are stored to files in PAF format\n"
            "\n"
            "  Overlap:\n"
            "    -K, --kmer-len <int>\n"
            "      default: 15\n"
            "      length of minimizers used to find overlaps\n"
            "    -W, --window-len <int>\n"
            "      default: 5\n"
            "      length of sliding window from which minimizers are sampled\n"
            "    -F, --frequency <double>\n"
            "      default: 0.001\n"
            "      threshold for ignoring most frequent minimizers\n"
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
            "    -t, --threads <int>\n"
            "      default: 1\n"
            "      number of threads\n"
            "\n"
            "  Layout:\n"
            "    -X, --max-overlaps <long unsigned int>\n"
            "      default: 32\n"
            "      maximum number of overlaps that will be taken during find overlaps and create piles stage\n"
            "    -D, --disagreement <double>\n"
            "      default: 0.1\n"
            "      maximal percentage of different anntoated bases in overlaps\n"
            "    -R --valid-region <int>\n"
            "      default: 4\n"
            "      overlap valid region size\n"
            "    --cleaning-rounds <int>\n"
            "      default: 0\n"
            "      number of cleaning rounds\n"
            "    --min-unitig-size <int>\n"
            "      default: 9999\n"
            "      minimal unitig length\n"
            "    --ploidy <int>\n"
            "      default: 2\n"
            "      number of haplotypes\n"
            "\n"
            "  Phasing:\n"
            "    --ultralong-phasing <string>\n"
            "       path to ul reads used for phasing\n"
            "\n"
            "  Consensus:\n"
            "    --polishing-rounds <int>\n"
            "      default: 0\n"
            "      number of times racon is invoked\n";

}

int ProcessParameters(int argc, char **argv, Program_Parameters& param) {
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case opt_error_corrected_reads:
        param.error_corrected_reads = optarg;
        break;
      case opt_split:
        param.split = std::atoi(optarg);
        break;
      case opt_graphical_fragment_assembly:
        param.gfa_path = optarg;
        break;
      case opt_resume:
        param.resume = true;
        break;
      case opt_disable_checkpoints:
        param.checkpoints = false;
        break;
      case opt_version:
        std::cout << RAVEN_VERSION << std::endl;
        return 0;
      case opt_help:
        Help();
        return 0;
      case opt_output:
        param.output_path = optarg;
        param.stdoutput = false;
        break;
      case opt_input_gfa:
        param.input_gfa_path = optarg;
        param.skip_contruction = true;
        break;
      case opt_paf:
        param.paf = true;
        break;
      case opt_kmer_len:
      case 'K':
        param.kmer_len = std::atoi(optarg);
        break;
      case opt_window_len:
      case 'W':
        param.window_len = std::atoi(optarg);
        break;
      case opt_frequency:
      case 'F':
        param.freq = std::atof(optarg);
        break;
      case opt_hpc:
      case 'H':
        param.hpc = true;
        break;
      case opt_bandwidth:
      case 'B':
        param.bandwidth = std::atoi(optarg);
        break;
      case opt_chain_n:
      case 'C':
        param.chain_n = std::atoi(optarg);
        break;
      case opt_match_n:
      case 'M':
        param.match_n = std::atoi(optarg);
        break;
      case opt_gap_size:
      case 'G':
        param.gap_size = std::atoi(optarg);
        break;
      case opt_threads:
      case 't':
        param.num_threads = atoi(optarg);
        break;
      case opt_max_overlaps:
      case 'X':
        param.max_overlaps = std::atof(optarg);
        break;
      case opt_disagreement:
      case 'D':
        param.disagreement = std::atof(optarg);
        break;
      case opt_valid_region:
      case 'R':
        param.valid_region_size = std::atoi(optarg);
        break;
      case opt_cleaning_rounds:
        break;
      case opt_min_unitig_size:
        param.min_unitig_size = std::atoi(optarg);
        break;
      case opt_ploidy:
        param.ploidy = std::atoi(optarg);
        break;
      case opt_ultralong_phasing:
        param.ul_read_path = optarg;
        break;
      case opt_polishing_rounds:
        param.num_polishing_rounds = atoi(optarg);
        break;
      default:
        return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind >= argc && !param.skip_contruction) {
    std::cerr << "[raven::] error: missing input file!" << std::endl;
    return 0;
  }

  param.sequence_path = argv[optind];

  return 1;
}
