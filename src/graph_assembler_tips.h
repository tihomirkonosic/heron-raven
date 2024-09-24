
#ifndef RAVEN_SRC_GRAPH_ASSEMBLER_TIPS_H_
#define RAVEN_SRC_GRAPH_ASSEMBLER_TIPS_H_

#include <cstdint>
#include "graph.hpp"

namespace raven {

class GraphAssemblerTips {
 public:
  explicit GraphAssemblerTips(Graph &graph);

  std::uint32_t RemoveTips();

 private:
  Graph &graph_;
};

}

#endif //RAVEN_SRC_GRAPH_ASSEMBLER_TIPS_H_
