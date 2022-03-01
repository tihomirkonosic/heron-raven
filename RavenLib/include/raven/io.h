#ifndef RAVEN_IO_H_
#define RAVEN_IO_H_

#include <string_view>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"

namespace raven {

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path);

}

#endif  // RAVEN_IO_H_
