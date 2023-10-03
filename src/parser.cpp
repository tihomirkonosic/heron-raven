
#include <iostream>
#include "parser.h"

namespace raven {

  std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> Parser::CreateParser(const std::string &path) {
    auto is_suffix = [](const std::string &s, const std::string &suff) {
      return s.size() < suff.size() ? false :
             s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
    };

    if (is_suffix(path, ".fasta") || is_suffix(path, ".fa") ||
        is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
      try {
        return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);  // NOLINT
      } catch (const std::invalid_argument &exception) {
        std::cerr << exception.what() << std::endl;
        return nullptr;
      }
    }
    if (is_suffix(path, ".fastq") || is_suffix(path, ".fq") ||
        is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
      try {
        return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path);  // NOLINT
      } catch (const std::invalid_argument &exception) {
        std::cerr << exception.what() << std::endl;
        return nullptr;
      }
    }

    std::cerr << "[raven::CreateParser] error: file " << path
              << " has unsupported format extension (valid extensions: .fasta, "
              << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
              << std::endl;
    return nullptr;
  }
} // raven