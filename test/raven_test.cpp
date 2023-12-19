// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"
#include "graph_constructor.h"
#include "graph_assembler.h"

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "edlib.h"  // NOLINT
#include "gtest/gtest.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace raven::test {

class RavenTest: public ::testing::Test {
 public:
  void SetUp() {
//    biosoup::NucleicAcid::num_objects = 0;
//    auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(  // NOLINT
//        TEST_DATA + std::string("ERA476754.fastq.gz"));
//    s = p->Parse(-1);
//
//    biosoup::NucleicAcid::num_objects = 0;
//    p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(  // NOLINT
//        TEST_DATA + std::string("NC_001416.fasta.gz"));
//    r = std::move(p->Parse(-1).front());
  }

//  static std::uint32_t EditDistance(
//      const std::string& lhs,
//      const std::string &rhs) {
//    EdlibAlignResult result = edlibAlign(
//        lhs.c_str(), lhs.size(),
//        rhs.c_str(), rhs.size(),
//        edlibDefaultAlignConfig());
//    std::uint32_t ed = result.editDistance;
//    edlibFreeAlignResult(result);
//    return ed;
//  }
//
//  std::vector<std::unique_ptr<biosoup::NucleicAcid>> s;
//  std::unique_ptr<biosoup::NucleicAcid> r;
};

//TEST_F(RavenTest, Assemble) {
//  Graph g{false, nullptr};
//  Graph_Constructor graph_constructor{g, nullptr};
//  graph_constructor.Construct(s);
//  Graph_Assembler graph_assembler{g, nullptr};
//  graph_assembler.Assemble();
//  //g.Polish(s, 3, -5, -4, 0, false, 0, 2);
//  auto u = std::move(g.GetUnitigs(true).front());
//  u->ReverseAndComplement();
//  EXPECT_EQ(1137, EditDistance(u->InflateData(), r->InflateData()));
//}
//
//TEST_F(RavenTest, Checkpoints) {
//  Graph g{false, nullptr};
//  Graph_Constructor graph_constructor{g, nullptr};
//  graph_constructor.Construct(s);
//  Graph_Assembler graph_assembler{g, nullptr};
//  graph_assembler.Assemble();
//  //g.Polish(s, 2, -5, -2, 0, false, 0, 2);
//  auto u = std::move(g.GetUnitigs(true).front());
//
//  SetUp();
//
//  g = {true, nullptr};
//  graph_constructor.Construct(s);
//
//  g = {true, nullptr};
//  g.Load();
//  graph_assembler.Assemble();
//
//  // g = {true, nullptr};
//  // g.Load();
//  // g.Polish(s, 2, -5, -2, 0, false, 0, 1);
//
//  // g = {true, nullptr};
//  // g.Load();
//  // g.Polish(s, 2, -5, -2, 0, false, 0, 2);
//
//  EXPECT_EQ(u->deflated_data, g.GetUnitigs(true).front()->deflated_data);
//}

    TEST(RavenTest, RemoveTips) {
      raven::Graph graph;

      raven::Graph_Assembler assembler(graph);

      std::uint32_t res = assembler.RemoveTips();
    }

}  // namespace raven::test
