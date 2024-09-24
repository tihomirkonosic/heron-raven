
#ifndef RAVEN_SRC_GRAPH_ASSEMBLER_BUBBLES_H_
#define RAVEN_SRC_GRAPH_ASSEMBLER_BUBBLES_H_

#include <cstdint>
#include "graph.hpp"

namespace raven {

class GraphAssemblerBubbles {
 public:
  explicit GraphAssemblerBubbles(Graph &graph);

  std::uint32_t RemoveBubbles();
  std::uint32_t RemoveSnpBubbles();

 private:
  Graph &graph_;
};

}

#endif //RAVEN_SRC_GRAPH_ASSEMBLER_BUBBLES_H_
