    #include <cstdint>
    #include "extended_overlap.h"
    #include <vector>
    
    
    struct edlib_align {
      std::uint32_t matches;
      std::uint32_t block_legnth;
      std::string cigar;
      std::int32_t edit_distance;
    };

    struct extended_overlap {
      biosoup::Overlap overlap;
      edlib_align edlib_alignment;
      std::uint32_t total_overlap_snps;
      std::uint32_t total_overlap_snp_mismatches;
      float identity;
      float heterozygosity_rate;
      std::uint32_t graph_overlap_type;
      std::uint8_t ol_class;
      bool ground_truth = NULL;
      std::vector<std::pair<std::uint32_t, std::uint32_t>> gap_positions; // odd entries are start positions, even entries are end positions
      std::uint32_t found_length = 0;
      std::uint32_t found_matches = 0;
    };