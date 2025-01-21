#ifndef RAVEN_EXTENDED_OVERLAP_H
#define RAVEN_EXTENDED_OVERLAP_H

#include "biosoup/overlap.hpp"


enum class OverlapType {
  perfect_heterozygous_high_match, // perfect(identity 99.5+) overlap with a high number of heterozygous SNPs
  perfect_heterozygous_high_mismatch, // perfect(identity 99.5+) overlap with a high number of heterozygous SNPs with mismatches
  perfect_heterozygous_low_match, // perfect(identity  99.5+) overlap with a low number of heterozygous SNPs
  perfect_heterozygous_low_mismatch, // perfect(identity  99.5+) overlap with a high number of heterozygous SNPs with mismatches
  perfect_homozygous, // perfect(identity  99.5+) overlap without SNPs

  high_heterozygous_high_match, // high(identity [95, 99.5]) overlap with a high number of heterozygous SNPs
  high_heterozygous_high_mismatch, // high(identity [95, 99.5]) overlap with a high number of heterozygous SNPs with mismatches
  high_heterozygous_low_match, // high(identity [95, 99.5]) overlap with a low number of heterozygous SNPs
  high_heterozygous_low_mismatch, // high(identity [95, 99.5]) overlap with a high number of heterozygous SNPs with mismatches
  high_homozygous, // high(identity [95, 99.5]) overlap without SNPs

  mid_heterozygous_high_match, // mid(identity [85, 95]) overlap with a high number of heterozygous SNPs
  mid_heterozygous_high_mismatch, // mid(identity [85, 95]) overlap with a high number of heterozygous SNPs with mismatches
  mid_heterozygous_low_match, // mid(identity [85, 95]) overlap with a low number of heterozygous SNPs
  mid_heterozygous_low_mismatch, // mid(identity [85, 95]) overlap with a high number of heterozygous SNPs with mismatches
  mid_homozygous, // mid(identity [85, 95]) overlap without SNPs

  low_heterozygous_high_match, // perfect(identity 99+) overlap with a high number of heterozygous SNPs
  low_heterozygous_high_mismatch,
  low_heterozygous_low_match, // perfect(identity 99+) overlap with a low number of heterozygous SNPs
  low_heterozygous_low_mismatch, // perfect(identity 99+) overlap with a high number of heterozygous SNPs with mismatches
  low_homozygous, // perfect(identity 99+) overlap without SNPs

  other,
  repetitive, // high order repeat overlap
  undefined
};

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
  OverlapType ol_type;
};




#endif // RAVEN_EXTENDED_OVERLAP_H


