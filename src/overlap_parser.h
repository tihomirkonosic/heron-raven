
#ifndef RAVEN_OVERLAP_PARSER_H
#define RAVEN_OVERLAP_PARSER_H

#include <memory>
#include <vector>
#include <cstring>
#include "zlib.h"

namespace raven {

  struct Overlap {
  public:
    Overlap() = default;

    Overlap(
      const char *q_name, std::uint32_t q_name_len,
      std::uint32_t q_len,
      std::uint32_t q_begin,
      std::uint32_t q_end,
      char orientation,
      const char *t_name, std::uint32_t t_name_len,
      std::uint32_t t_len,
      std::uint32_t t_begin,
      std::uint32_t t_end,
      std::uint32_t score,
      std::uint32_t overlap_len,
      std::uint32_t quality)
      :
      lhs_id(0),
      lhs_begin(q_begin),
      lhs_end(q_end),
      rhs_id(0),
      rhs_begin(t_begin),
      rhs_end(t_end),
      score(score),
      strand(orientation == '+'),
      alignment(),
      q_name(q_name, q_name_len),
      q_len(q_len),
      t_name(t_name, t_name_len),
      t_len(t_len),
      overlap_len(overlap_len),
      paf_quality(quality) {}

    Overlap(
      const char *q_name, std::uint32_t q_name_len,
      std::uint32_t flag,
      const char *t_name, std::uint32_t t_name_len,
      std::uint32_t t_begin,
      std::uint32_t map_quality,
      const char *cigar, std::uint32_t cigar_len,
      const char *t_next_name, std::uint32_t t_next_name_len,
      std::uint32_t t_next_begin,
      std::uint32_t template_len,
      const char *data, std::uint32_t data_len,
      const char *quality, std::uint32_t quality_len)
      :
      lhs_id(0),
      lhs_begin(0),
      lhs_end(0),
      rhs_id(0),
      rhs_begin(t_begin),
      rhs_end(0),
      score(0),
      strand(strand),
      alignment(cigar, cigar_len),
      q_name(q_name, q_name_len),
      sam_flag(flag),
      t_name(t_name, t_name_len),
      map_quality(map_quality),
      t_next_name(t_next_name, t_next_name_len),
      t_next_begin(t_next_begin),
      template_len(template_len),
      data(data, data_len),
      sam_quality(quality, quality_len) {}

    Overlap(const Overlap &) = default;

    Overlap &operator=(const Overlap &) = default;

    Overlap(Overlap &&) = default;

    Overlap &operator=(Overlap &&) = default;

    ~Overlap() = default;

    std::uint32_t lhs_id;
    std::uint32_t lhs_begin;
    std::uint32_t lhs_end;
    std::uint32_t rhs_id;
    std::uint32_t rhs_begin;
    std::uint32_t rhs_end;
    std::uint32_t score;  // based on k-mer matches or alignment score
    bool strand;  // (optional) Watson-Crick strand
    std::string alignment;  // (optional) cigar string

    //PAF
    std::string q_name;
    std::uint32_t q_len;
    std::string t_name;
    std::uint32_t t_len;
    std::uint32_t overlap_len;
    std::uint32_t paf_quality;

    //SAM
    std::uint32_t sam_flag;
    std::uint32_t map_quality;
    std::string t_next_name;
    std::uint32_t t_next_begin;
    std::uint32_t template_len;
    std::string data;
    std::string sam_quality;
  };

  class OverlapParser {
  public:
    OverlapParser(const OverlapParser &) = delete;

    OverlapParser &operator=(const OverlapParser &) = delete;

    OverlapParser(OverlapParser &&) = delete;

    OverlapParser &operator=(OverlapParser &&) = delete;

    OverlapParser(gzFile file)
      : file_(file, gzclose),
        buffer_(storage_size_, 0),  // 64 kB
        buffer_ptr_(0),
        buffer_bytes_(0),
        line_(storage_size_, 0),
        line_ptr_(0) {}

    virtual ~OverlapParser() {}

    static std::unique_ptr<OverlapParser> Create(const std::string &path);

    void ResetFilePointer();

    std::vector<std::unique_ptr<Overlap>> ParsePaf(std::uint64_t bytes, bool shorten_names = true);
    std::vector<std::unique_ptr<Overlap>> ParseSam(std::uint64_t bytes, bool shorten_names = true);

  private:
    std::vector<std::unique_ptr<Overlap>> Parse(std::uint64_t bytes, bool shorten_names,
                                          void (OverlapParser::*add_overlap_to_list)(std::vector<std::unique_ptr<Overlap>> &dst, bool shorten_names));

    void AddPafOverlapToList(std::vector<std::unique_ptr<Overlap>> &dst, bool shorten_names);
    void AddSamOverlapToList(std::vector<std::unique_ptr<Overlap>> &dst, bool shorten_names);

    bool ReadBuffer();

    void ReadLineFromBuffer(std::uint32_t count, bool strip = false);

    void TerminateLine(std::uint32_t i);

    void ClearLine();

    static std::uint32_t RightStrip(const char *str, std::uint32_t str_len);

    static std::uint32_t Shorten(const char *str, std::uint32_t str_len);

    std::unique_ptr<gzFile_s, int (*)(gzFile)> file_;
    std::vector<char> buffer_;
    std::uint32_t buffer_ptr_;
    std::uint32_t buffer_bytes_;
    std::vector<char> line_;
    std::uint32_t line_ptr_;

    const std::uint32_t storage_size_ = 65536; // 64 kB
  };

} // raven

#endif //RAVEN_OVERLAP_PARSER_H
