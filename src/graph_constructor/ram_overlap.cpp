#include "ram_overlap.h"
#include <algorithm>



// ram_overlap::ram_overlap() {
//     n_fragments = 1;
// };

extended_overlap ram_overlap::resolve() {
    // is a contained in b
    auto is_contained = [](const biosoup::Overlap& a, const biosoup::Overlap& b) -> bool {
        return (a.lhs_begin >= b.lhs_begin && a.lhs_end <= b.lhs_end); //&&
               //(a.rhs_begin >= b.rhs_begin && a.rhs_end <= b.rhs_end);
    };

    std::sort(overlap_fragments.begin(), overlap_fragments.end(),
        [](const extended_overlap& lhs, const extended_overlap& rhs) -> bool {
            return lhs.overlap.lhs_begin < rhs.overlap.lhs_begin;
        }
    );

    // check if two overlaps are overlapping, also check for overlaps on the target side
    auto is_overlapping = [](const biosoup::Overlap& a, const biosoup::Overlap& b) -> bool {
        if((std::max(a.lhs_begin, b.lhs_begin) <= std::min(a.lhs_end, b.lhs_end)) &&
               (std::max(a.rhs_begin, b.rhs_begin) <= std::min(a.rhs_end, b.rhs_end))){
                return 1;
            } else if((std::max(a.lhs_begin, b.lhs_begin) <= std::min(a.lhs_end, b.lhs_end)) &&
               !(std::max(a.rhs_begin, b.rhs_begin) <= std::min(a.rhs_end, b.rhs_end))){
                return 2;
            } else {
                return 0;
               }
    };

    auto normalize = [](biosoup::Overlap& o) -> void {
        if (o.lhs_end < o.lhs_begin) std::swap(o.lhs_begin, o.lhs_end);
        if (o.rhs_end < o.rhs_begin) std::swap(o.rhs_begin, o.rhs_end);
    };

    auto axis_overlap = [](std::int64_t a0, std::int64_t a1,
                        std::int64_t b0, std::int64_t b1,
                        std::int64_t tol = 0) -> bool{
    // overlap if max(start) <= min(end) + tol
        return std::max(a0, b0) <= std::min(a1, b1) + tol;
    };

    auto merge_into = [](biosoup::Overlap& anchor,
                       const biosoup::Overlap& other) -> void {
    // Prefix/suffix (or general) merge: expand to cover both on both axes
        anchor.lhs_begin = std::min(anchor.lhs_begin, other.lhs_begin);
        anchor.lhs_end   = std::max(anchor.lhs_end,   other.lhs_end);
        anchor.rhs_begin = std::min(anchor.rhs_begin, other.rhs_begin);
        anchor.rhs_end   = std::max(anchor.rhs_end,   other.rhs_end);
        anchor.score += other.score;
    // If you keep additional fields (score, matches, length, etc.), update them here.
    // e.g., anchor.score += other.score; or anchor.len = (anchor.lhs_end - anchor.lhs_begin) ...
    };

    enum class OverlapType : std::uint8_t {
        None = 0,
        QueryOnly = 2,
        Both = 1
    };

    auto classify_overlap = [&axis_overlap](const biosoup::Overlap& a,
                                    const biosoup::Overlap& b,
                                    std::int64_t tol = 0) -> OverlapType {
        bool q = axis_overlap(a.lhs_begin, a.lhs_end, b.lhs_begin, b.lhs_end, tol);
        bool t = axis_overlap(a.rhs_begin, a.rhs_end, b.rhs_begin, b.rhs_end, tol);
        if (q && t) return OverlapType::Both;
        if (q && !t) return OverlapType::QueryOnly;
        return OverlapType::None;
    };

    auto merge_prefix_suffix_overlaps = [&normalize, &merge_into, &classify_overlap](std::vector<extended_overlap>& overlap_fragments,
                                  std::int64_t tol = 0) -> void {
    // Optional: normalize all intervals once
    for (auto& e : overlap_fragments) normalize(e.overlap);

    // We’ll treat each item as a potential anchor and greedily merge any
    // follower that overlaps on both axes.
    for (std::size_t i = 0; i < overlap_fragments.size(); /* no ++ here */) {
        auto& anchor = overlap_fragments[i].overlap;
        bool merged_any = false;

        for (std::size_t j = i + 1; j < overlap_fragments.size(); /* no ++ here */) {
            auto type = classify_overlap(anchor, overlap_fragments[j].overlap, tol);

            if (type == OverlapType::Both) {
                // Merge j into i, erase j, and keep scanning j at same index
                merge_into(overlap_fragments[i].overlap, overlap_fragments[j].overlap);
                overlap_fragments.erase(overlap_fragments.begin() + j);
                merged_any = true;
                // Do not increment j; the next element just slid into position j
            } else {
                // Skip non-mergeable (QueryOnly or None) for now
                ++j;
            }
        }

        // if (merged_any) {
        //     // We expanded anchor; optionally, you could re-run normalization (not needed here).
        //     // Keep i at same position only if you need to re-scan prior elements;
        //     // since we only ever look forward, just advance i.
        // }
        ++i;
    }
    };


    // resolve contained overlaps
    std::vector<extended_overlap> resolved_fragments;
    resolved_fragments.reserve(overlap_fragments.size());

    resolved_fragments.emplace_back(overlap_fragments.front());
    if (n_fragments == 1){
        extended_overlap eo = resolved_fragments.front();
        return eo;
    } else {
        for(auto& ov: overlap_fragments){
            bool is_contained_flag = false;
            for(auto& res_ov: resolved_fragments){
                if(is_contained(ov.overlap, res_ov.overlap)){
                    is_contained_flag = true;
                    break;
                }
            }
            if(!is_contained_flag){
                resolved_fragments.emplace_back(ov);
            }
        }
    }

    overlap_fragments.swap(resolved_fragments);
    n_fragments = overlap_fragments.size();

    
    extended_overlap eo = overlap_fragments.front();
    return eo;

    if(n_fragments == 1){
        return overlap_fragments.front();
    }
    // // merge fragments
    // // Resolve overlapping
    merge_prefix_suffix_overlaps(overlap_fragments, 100); // tolerance of 100 bp


};


// class ram_overlap {
// public:
//     std::vector<biosoup::Overlap> overlaps;
//     std::uint32_t lhs_id;
//     std::uint32_t rhs_id;

//     std::uint32_t lhs_begin;
//     std::uint32_t lhs_end;
//     std::uint32_t rhs_begin;
//     std::uint32_t rhs_end;

//     std::uint16_t n_fragments = 0; // number of fragments in the overlap
//     std::vector<std::uint16_t> fragment_gaps; //gaps between fragments

// };