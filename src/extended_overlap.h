#ifndef RAVEN_EXTENDED_OVERLAP_H
#define RAVEN_EXTENDED_OVERLAP_H

#include "biosoup/overlap.hpp"

struct edlib_align {
  std::uint32_t matches;
  std::uint32_t block_length;
  std::string cigar;
  std::int32_t edit_distance;
};

struct extended_overlap {
  biosoup::Overlap overlap;
  edlib_align edlib_alignment;
};


#endif // RAVEN_EXTENDED_OVERLAP_H


