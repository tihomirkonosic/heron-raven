
#ifndef RAVEN_SRC_GRAPH_ASSEMBLER_TRANSITIVE_H_
#define RAVEN_SRC_GRAPH_ASSEMBLER_TRANSITIVE_H_

#include <cstdint>
#include "graph.hpp"

namespace raven {

class GraphAssemblerTransitive {
 public:
  explicit GraphAssemblerTransitive(Graph &graph);

  std::uint32_t RemoveTransitiveEdges();

 private:
  Graph &graph_;
};

}

#endif //RAVEN_SRC_GRAPH_ASSEMBLER_TRANSITIVE_H_
