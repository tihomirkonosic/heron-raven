
#include <fstream>
#include <random>
#include "graph_assembler.h"
#include "biosoup/timer.hpp"
#include "graph_assembler_transitive.h"
#include "graph_assembler_tips.h"
#include "graph_assembler_bubbles.h"

namespace raven {
  Graph_Assembler::Graph_Assembler(Graph &graph, Program_Parameters &param, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
      : graph_(graph), param_(param), thread_pool_(thread_pool ?
                                    thread_pool :
                                    std::make_shared<thread_pool::ThreadPool>(1)) {
  }

  void Graph_Assembler::Assemble() {
    if (graph_.stage() != Graph_Stage::Assemble_Transitive_Edges
        && graph_.stage() != Graph_Stage::Assemble_Tips_Bubbles
        && graph_.stage() != Graph_Stage::Assemble_Force_Directed) {
      return;
    }

    biosoup::Timer timer{};

    // remove transitive edges
    if (graph_.stage() == Graph_Stage::Assemble_Transitive_Edges) {
      timer.Start();

      RemoveTransitiveEdges();

      std::cerr << "[raven::Graph::Assemble] removed transitive edges "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa(param_.gfa_after_transitive_filename, true);
    }

    if (graph_.stage() == Graph_Stage::Assemble_Transitive_Edges) {
      timer.Start();

      //CreateUnitigGraph();

      std::cerr << "[raven::Graph::Assemble] created bubble chain "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      // PrintGfa("after_bubble_chain.gfa");
    }

    // checkpoint
    if (graph_.stage() == Graph_Stage::Assemble_Transitive_Edges) {
      graph_.advance_stage();
      if (graph_.use_checkpoints()) {
        timer.Start();
        //Store();
        std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }

    // remove tips and bubbles
    if (graph_.stage() == Graph_Stage::Assemble_Tips_Bubbles) {
      timer.Start();

      while (true) {
        std::uint32_t num_changes = RemoveTips();
        num_changes += RemoveBubbles();
        if (num_changes == 0) {
          break;
        }
      }

      std::cerr << "[raven::Graph::Assemble] removed tips and bubbles "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa(param_.gfa_after_bubble_filename, false);
    }

    // checkpoint
    if (graph_.stage() == Graph_Stage::Assemble_Tips_Bubbles) {
      graph_.advance_stage();
      if (graph_.use_checkpoints()) {
        timer.Start();
        //Store();
        std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }

    // remove long edges
    if (graph_.stage() == Graph_Stage::Assemble_Force_Directed) {
      timer.Start();

      if (graph_.annotations_.empty()) {
        graph_.CreateUnitigs(42);  // speed up force directed layout
      }
      RemoveLongEdges(16);

      std::cerr << "[raven::Graph::Assemble] removed long edges "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa(param_.gfa_after_force_filename, false);
    }

    // checkpoint
    if (graph_.stage() == Graph_Stage::Assemble_Force_Directed) {
      graph_.advance_stage();
      if (graph_.use_checkpoints()) {
        timer.Start();
        //Store();
        std::cerr << "[raven::Graph::Assemble] reached checkpoint "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;
      }
    }

    std::cerr << "[raven::Graph::Assemble] "
              << std::fixed << timer.elapsed_time() << "s"
              << std::endl;
  }

  void Graph_Assembler::AssembleDiploids() {
    if (graph_.stage() != Graph_Stage::Assemble_Diploids) {
      return;
    }
    biosoup::Timer timer{};

    timer.Start();

    graph_.DuplicateGraph();

    std::cerr << "[raven::Graph::AssembleDiploids] Duplicated graph "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    while (true) {
      std::uint32_t num_changes = RemoveSnpBubbles();

      if (num_changes == 0) {
        break;
      }
    }

    std::cerr << "[raven::Graph::AssembleDiploids] removed SNP bubbles "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    timer.Start();

    //graph_.SalvageHaplotypesPrimary();
    //graph_.SalvageHaplotypesAlternative();
    graph_.SalvageHaplotypes();

    std::cerr << "[raven::Graph::AssembleDiploids] Created haplotypes "
              << std::fixed << timer.Stop() << "s"
              << std::endl;

    std::cerr << "[raven::Graph::AssembleDiploids] "
              << std::fixed << timer.elapsed_time() << "s"
              << std::endl;
  }

  void Graph_Assembler::AssembleHaploids(){
      if (graph_.stage() != Graph_Stage::Assemble_Diploids) {
      return;
    }
    biosoup::Timer timer{};

    timer.Start();

    graph_.SalvageHaplotypes();

    std::cerr << "[raven::Graph::AssembleHaploids] Created haplotypes "
              << std::fixed << timer.Stop() << "s"
              << std::endl;
  }

  void Graph_Assembler::UlAssemble(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_sequences) {
    if (graph_.stage() != Graph_Stage::Assemble_Transitive_Edges
        && graph_.stage() != Graph_Stage::Assemble_Tips_Bubbles
        && graph_.stage() != Graph_Stage::Assemble_Force_Directed) {
      return;
    }

    biosoup::Timer timer{};
    graph_.PrintGfa(param_.gfa_after_construction_filename, false);

    // remove transitive edges
    if (graph_.stage() == Graph_Stage::Assemble_Transitive_Edges) {
      timer.Start();

      RemoveTransitiveEdges();

      std::cerr << "[raven::Graph::Assemble] removed transitive edges "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      graph_.PrintGfa(param_.gfa_after_transitive_filename, false);
    }

    if (graph_.stage() == Graph_Stage::Assemble_Transitive_Edges) {
      timer.Start();

      graph_.CreateUnitigGraph(param_.gfa_unitig_graph_filename);

      std::cerr << "[raven::Graph::Assemble] created bubble chain "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

      // PrintGfa("after_bubble_chain.gfa");
    }

    if (graph_.stage() == Graph_Stage::Assemble_Transitive_Edges) {
      timer.Start();

      ResolveGraphWithUl(ul_sequences);

      std::cerr << "[raven::Graph::Assemble] created bubble chain "
                << std::fixed << timer.Stop() << "s"
                << std::endl;

    }

  }

  std::uint32_t Graph_Assembler::RemoveTransitiveEdges() {
    GraphAssemblerTransitive asm_transitive{graph_};
    return asm_transitive.RemoveTransitiveEdges();
  }

  std::uint32_t Graph_Assembler::RemoveTips() {
    GraphAssemblerTips asm_tips{graph_};
    return asm_tips.RemoveTips();
  }

  std::uint32_t Graph_Assembler::RemoveBubbles() {
    GraphAssemblerBubbles asm_bubbles{graph_};
    return asm_bubbles.RemoveBubbles();
  }

  std::uint32_t Graph_Assembler::RemoveSnpBubbles() {
    GraphAssemblerBubbles asm_bubbles{graph_};
    return asm_bubbles.RemoveSnpBubbles();
  }

  std::uint32_t Graph_Assembler::RemoveLongEdges(std::uint32_t num_rounds) {
    std::uint32_t num_long_edges = 0;

    for (std::uint32_t i = 0; i < num_rounds; ++i) {
      CreateForceDirectedLayout();

      std::unordered_set<std::uint32_t> marked_edges;
      for (const auto &it: graph_.nodes_) {
        if (it == nullptr || it->outdegree() < 2) {
          continue;
        }
        for (auto jt: it->outedges) {
          for (auto kt: it->outedges) {
            if (jt != kt && jt->weight * 2.0 < kt->weight) {  // TODO(rvaser)
              marked_edges.emplace(kt->id);
              marked_edges.emplace(kt->pair->id);
            }
          }
        }
      }
      graph_.RemoveEdges(marked_edges);
      num_long_edges += marked_edges.size() / 2;

      RemoveTips();
    }

    return num_long_edges;
  }

  void Graph_Assembler::CreateForceDirectedLayout(const std::string &path) {
    std::ofstream os;
    bool is_first = true;
    if (!path.empty()) {
      os.open(path);
      os << "{" << std::endl;
    }

    std::vector<std::unordered_set<std::uint32_t>> components;
    std::vector<char> is_visited(graph_.nodes_.size(), 0);
    for (std::uint32_t i = 0; i < graph_.nodes_.size(); ++i) {
      if (graph_.nodes_[i] == nullptr || is_visited[i]) {
        continue;
      }

      components.resize(components.size() + 1);

      std::deque<std::uint32_t> que = {i};
      while (!que.empty()) {
        std::uint32_t j = que.front();
        que.pop_front();

        if (is_visited[j]) {
          continue;
        }
        const auto &node = graph_.nodes_[j];
        is_visited[node->id] = 1;
        is_visited[node->pair->id] = 1;
        components.back().emplace((node->id >> 1) << 1);

        for (auto it: node->inedges) {
          que.emplace_back(it->tail->id);
        }
        for (auto it: node->outedges) {
          que.emplace_back(it->head->id);
        }
      }
    }
    std::vector<char>().swap(is_visited);

    std::sort(components.begin(), components.end(),
              [](const std::unordered_set<std::uint32_t> &lhs,
                 const std::unordered_set<std::uint32_t> &rhs) {
                return lhs.size() > rhs.size();
              });

    static std::uint64_t seed = 21;
    seed <<= 1;

    std::mt19937 generator(seed);
    std::uniform_real_distribution<> distribution(0., 1.);

    struct Point {
      Point() = default;

      Point(double x, double y)
          : x(x),
            y(y) {}

      bool operator==(const Point &other) const {
        return x == other.x && y == other.y;
      }

      Point operator+(const Point &other) const {
        return Point(x + other.x, y + other.y);
      }

      Point &operator+=(const Point &other) {
        x += other.x;
        y += other.y;
        return *this;
      }

      Point operator-(const Point &other) const {
        return Point(x - other.x, y - other.y);
      }

      Point operator*(double c) const {
        return Point(x * c, y * c);
      }

      Point &operator/=(double c) {
        x /= c;
        y /= c;
        return *this;
      }

      double Norm() const {
        return sqrt(x * x + y * y);
      }

      double x;
      double y;
    };

    struct Quadtree {
      Quadtree(Point nucleus, double width)
          : nucleus(nucleus),
            width(width),
            center(0, 0),
            mass(0),
            subtrees() {
      }

      bool Add(const Point &p) {
        if (nucleus.x - width > p.x || p.x > nucleus.x + width ||
            nucleus.y - width > p.y || p.y > nucleus.y + width) {
          return false;
        }
        ++mass;
        if (mass == 1) {
          center = p;
        } else if (subtrees.empty()) {
          if (center == p) {
            return true;
          }
          double w = width / 2;
          subtrees.emplace_back(Point(nucleus.x + w, nucleus.y + w), w);
          subtrees.emplace_back(Point(nucleus.x - w, nucleus.y + w), w);
          subtrees.emplace_back(Point(nucleus.x - w, nucleus.y - w), w);
          subtrees.emplace_back(Point(nucleus.x + w, nucleus.y - w), w);
          for (auto &it: subtrees) {
            if (it.Add(center)) {
              break;
            }
          }
        }
        for (auto &it: subtrees) {
          if (it.Add(p)) {
            break;
          }
        }
        return true;
      }

      void Centre() {
        if (subtrees.empty()) {
          return;
        }
        center = Point(0, 0);
        for (auto &it: subtrees) {
          it.Centre();
          center += it.center * it.mass;
        }
        center /= mass;
      }

      Point Force(const Point &p, double k) const {
        auto delta = p - center;
        auto distance = delta.Norm();
        if (width * 2 / distance < 1) {
          return delta * (mass * (k * k) / (distance * distance));
        }
        delta = Point(0, 0);
        for (const auto &it: subtrees) {
          delta += it.Force(p, k);
        }
        return delta;
      }

      Point nucleus;
      double width;
      Point center;
      std::uint32_t mass;
      std::vector<Quadtree> subtrees;
    };

    std::uint32_t c = 0;
    for (const auto &component: components) {
      if (component.size() < 6) {
        continue;
      }

      bool has_junctions = false;
      for (const auto &it: component) {
        if (graph_.nodes_[it]->is_junction()) {
          has_junctions = true;
          break;
        }
      }
      if (has_junctions == false) {
        continue;
      }

      // update transitive edges
      for (const auto &n: component) {
        std::unordered_set<std::uint32_t> valid;
        for (const auto &m: graph_.nodes_[n]->transitive) {
          if (component.find(m) != component.end()) {
            valid.emplace(m);
          }
        }
        graph_.nodes_[n]->transitive.swap(valid);
      }

      std::uint32_t num_iterations = 100;
      double k = sqrt(1. / static_cast<double>(component.size()));
      double t = 0.1;
      double dt = t / static_cast<double>(num_iterations + 1);

      std::vector<Point> points(graph_.nodes_.size());
      for (const auto &it: component) {
        points[it].x = distribution(generator);
        points[it].y = distribution(generator);
      }

      for (std::uint32_t i = 0; i < num_iterations; ++i) {
        Point x = {0, 0}, y = {0, 0};
        for (const auto &n: component) {
          x.x = std::min(x.x, points[n].x);
          x.y = std::max(x.y, points[n].x);
          y.x = std::min(y.x, points[n].y);
          y.y = std::max(y.y, points[n].y);
        }
        double w = (x.y - x.x) / 2, h = (y.y - y.x) / 2;

        Quadtree tree(Point(x.x + w, y.x + h), std::max(w, h) + 0.01);
        for (const auto &n: component) {
          tree.Add(points[n]);
        }
        tree.Centre();

        std::vector<std::future<void>> thread_futures;
        std::vector<Point> displacements(graph_.nodes_.size(), Point(0, 0));

        auto thread_task = [&](std::uint32_t n) -> void {
          auto displacement = tree.Force(points[n], k);
          for (auto e: graph_.nodes_[n]->inedges) {
            auto m = (e->tail->id >> 1) << 1;
            auto delta = points[n] - points[m];
            auto distance = delta.Norm();
            if (distance < 0.01) {
              distance = 0.01;
            }
            displacement += delta * (-1. * distance / k);
          }
          for (auto e: graph_.nodes_[n]->outedges) {
            auto m = (e->head->id >> 1) << 1;
            auto delta = points[n] - points[m];
            auto distance = delta.Norm();
            if (distance < 0.01) {
              distance = 0.01;
            }
            displacement += delta * (-1. * distance / k);
          }
          for (const auto &m: graph_.nodes_[n]->transitive) {
            auto delta = points[n] - points[m];
            auto distance = delta.Norm();
            if (distance < 0.01) {
              distance = 0.01;
            }
            displacement += delta * (-1. * distance / k);
          }
          auto length = displacement.Norm();
          if (length < 0.01) {
            length = 0.1;
          }
          displacements[n] = displacement * (t / length);
          return;
        };

        for (const auto &n: component) {
          thread_futures.emplace_back(thread_pool_->Submit(thread_task, n));
        }
        for (const auto &it: thread_futures) {
          it.wait();
        }
        for (const auto &n: component) {
          points[n] += displacements[n];
        }

        t -= dt;
      }

      for (const auto &it: graph_.edges_) {
        if (it == nullptr || it->id & 1) {
          continue;
        }
        auto n = (it->tail->id >> 1) << 1;
        auto m = (it->head->id >> 1) << 1;

        if (component.find(n) != component.end() &&
            component.find(m) != component.end()) {
          it->weight = (points[n] - points[m]).Norm();
          it->pair->weight = it->weight;
        }
      }

      if (!path.empty()) {
        if (!is_first) {
          os << "," << std::endl;
        }
        is_first = false;

        os << "    \"component_" << c++ << "\": {" << std::endl;

        bool is_first_node = true;
        os << "      \"nodes\": {" << std::endl;
        for (const auto &it: component) {
          if (!is_first_node) {
            os << "," << std::endl;
          }
          is_first_node = false;
          os << "        \"" << it << "\": [";
          os << points[it].x << ", ";
          os << points[it].y << ", ";
          os << (graph_.nodes_[it]->is_junction() ? 1 : 0) << ", ";
          os << graph_.nodes_[it]->count << "]";
        }
        os << std::endl << "      }," << std::endl;

        bool is_first_edge = true;
        os << "      \"edges\": [" << std::endl;
        for (const auto &it: component) {
          for (auto e: graph_.nodes_[it]->inedges) {
            auto o = (e->tail->id >> 1) << 1;
            if (it < o) {
              continue;
            }
            if (!is_first_edge) {
              os << "," << std::endl;
            }
            is_first_edge = false;
            os << "        [\"" << it << "\", \"" << o << "\", 0]";
          }
          for (auto e: graph_.nodes_[it]->outedges) {
            auto o = (e->head->id >> 1) << 1;
            if (it < o) {
              continue;
            }
            if (!is_first_edge) {
              os << "," << std::endl;
            }
            is_first_edge = false;
            os << "        [\"" << it << "\", \"" << o << "\", 0]";
          }
          for (const auto &o: graph_.nodes_[it]->transitive) {
            if (it < o) {
              continue;
            }
            if (!is_first_edge) {
              os << "," << std::endl;
            }
            is_first_edge = false;
            os << "        [\"" << it << "\", \"" << o << "\", 1]";
          }
        }
        os << std::endl << "      ]" << std::endl;
        os << "    }";
      }
    }

    if (!path.empty()) {
      os << std::endl << "}";
      os << std::endl;
      os.close();
    }
  }

  void Graph_Assembler::ResolveGraphWithUl(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_reads) {
    //get unitig sequences
    //
  }
} // raven