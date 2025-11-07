#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>        // std::system
#include <algorithm>

#include "kmer_api.h"     // KMC: CKmerAPI
#include "kmc_file.h"     // KMC: CKMCFile

#ifndef KMC_BIN_DIR
#  define KMC_BIN_DIR "kmc_bin" // fallback; replaced by CMake define
#endif

namespace kmcwrap {

// Build a KMC DB from FASTA/FASTQ files.
// - inputs: list of paths (can be gzipped)
// - db_prefix: directory/prefix for output (creates files <prefix>.kmc_pre/.kmc_suf)
// - tmp_dir: temporary directory for KMC (created if needed)
// - k: k-mer length (supports >31)
// - mem_gb: RAM limit for KMC in GB
// - threads: number of threads
inline void build_database(const std::vector<std::string>& inputs,
                           const std::string& db_prefix,
                           const std::string& tmp_dir,
                           int k,
                           int mem_gb = 8,
                           int threads = 8,
                           bool canonical = true,
                           uint32_t min_count = 1,
                           uint32_t max_count = 0) // 0 = no upper cap
{
    if (inputs.empty()) throw std::runtime_error("KMC: no input files");

    // KMC expects either a single file or a file-list prefixed by '@'
    std::string list_path = db_prefix + ".kmc_input.lst";
    {
        std::ofstream out(list_path);
        if (!out) throw std::runtime_error("KMC: cannot write list file " + list_path);
        for (auto& p : inputs) out << p << "\n";
    }

    // Build command; key options:
    //  -k<k>         k-mer length
    //  -t<threads>   threads
    //  -m<mem_gb>    memory (GB)
    //  -ci<min>      min count
    //  -cx<max>      max count (0 = unlimited; we omit if 0)
    //  -fa/-fq       auto-detected normally; not required
    //  -b            both strands counted canonically (optional, depends on KMC version)
    std::ostringstream cmd;
    cmd << KMC_BIN_DIR << "/kmc"
        << " -k" << k
        << " -t" << threads
        << " -m" << mem_gb
        << " -ci" << min_count;
    if (max_count > 0) cmd << " -cx" << max_count;
    if (canonical) cmd << " -b"; // count canonical (F/R) into one (if supported)
    cmd << " @" << list_path
        << " " << db_prefix
        << " " << tmp_dir;

    int ret = std::system(cmd.str().c_str());
    if (ret != 0) throw std::runtime_error("KMC build failed: " + cmd.str());
}

// Produce a histogram file using kmc_tools, then parse it.
// Returns a vector where each element is (count, how_many_kmers_with_that_count).
inline std::vector<std::pair<uint64_t,uint64_t>>
histogram(const std::string& db_prefix, const std::string& out_hist_path = "")
{
    std::string hist_path = out_hist_path.empty()
                          ? (db_prefix + ".hist.txt")
                          : out_hist_path;

    std::ostringstream cmd;
    cmd << KMC_BIN_DIR << "/kmc_tools"
        << " transform " << db_prefix
        << " histogram " << hist_path;
    int ret = std::system(cmd.str().c_str());
    if (ret != 0) throw std::runtime_error("kmc_tools histogram failed: " + cmd.str());

    std::ifstream in(hist_path);
    if (!in) throw std::runtime_error("cannot read histogram " + hist_path);

    std::vector<std::pair<uint64_t,uint64_t>> hist;
    hist.reserve(1<<20);
    uint64_t c=0, n=0;
    while (in >> c >> n) hist.emplace_back(c, n);
    return hist;
}

// Open a KMC DB for random access queries (fast).
class Database {
public:
    explicit Database(const std::string& db_prefix) {
        if (!db_.OpenForRA(db_prefix)) {
            throw std::runtime_error("KMC: failed to open DB " + db_prefix);
        }
        k_ = db_.KmerLength();
    }
    ~Database() { db_.Close(); }

    // Return count of a k-mer string (A/C/G/T). If canonical=true, we’ll
    // query both orientations by taking the canonical representation.
    uint32_t count(const std::string& kmer, bool canonical=true) {
        if ((int)kmer.size() != k_) throw std::runtime_error("KMC: wrong k length for query");

        CKmerAPI q(k_);
        // KMC API expects uppercase DNA letters
        std::string s = kmer;
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return (char)std::toupper(c); });
        q.from_string(s.c_str());

        if (canonical) q.canonicalize();

        uint32_t c = 0;
        db_.CheckKmer(q, c);
        return c;
    }

private:
    CKMCFile db_;
    int k_ = 0;
};

} // namespace kmcwrap
