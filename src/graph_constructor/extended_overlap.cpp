    #include <cstdint>
    #include "extended_overlap.h"
    
    
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
      OverlapType ol_type = OverlapType::undefined;
      std::uint8_t ol_class;
      bool ground_truth = NULL;
    };