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
  std::uint32_t total_overlap_snps;
  std::uint32_t total_overlap_snp_mismatches;
};

enum class OverlapType {
  perfect_heterozygous, // perfect(identity 99+) overlap with heterozygous SNPs
  perfect_homozygous, // perfect(identity 99+) overlap without SNPs
  high_heterozygous, // lower identity(90-99%) overlap with higher number of heterozygous SNPs(10+)
  low_heterozygous, // lower identity(90-99%) overlap with lower number of heterozygous SNPs(10+)
  high_homozygous, // lower identity overlap
  low_homozygous, // lower identity overlap
  repettive // high order repeat overlap
};

#endif // RAVEN_EXTENDED_OVERLAP_H


