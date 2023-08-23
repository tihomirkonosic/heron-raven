//
// Created by Tiho on 9.8.2023.
//

#include "unitig.h"

namespace raven {

    std::uint32_t Unitig::CreateUnitigs(std::uint32_t epsilon) {
        std::unordered_set<std::uint32_t> marked_edges;
        std::vector<std::shared_ptr<Node>> unitigs;
        std::vector<std::shared_ptr<Edge>> unitig_edges;
        std::vector<std::uint32_t> node_updates(nodes_.size(), 0);
        std::vector<char> is_visited(nodes_.size(), 0);

        for (const auto& it : nodes_) {
            if (it == nullptr || is_visited[it->id] || it->is_junction()) {
                continue;
            }

            std::uint32_t extension = 1;

            bool is_circular = false;
            auto begin = it.get();
            while (!begin->is_junction()) {  // extend left
                is_visited[begin->id] = 1;
                is_visited[begin->pair->id] = 1;
                if (begin->indegree() == 0 ||
                    begin->inedges.front()->tail->is_junction()) {
                    break;
                }
                begin = begin->inedges.front()->tail;
                ++extension;
                if (begin == it.get()) {
                    is_circular = true;
                    break;
                }
            }

            auto end = it.get();
            while (!end->is_junction()) {  // extend right
                is_visited[end->id] = 1;
                is_visited[end->pair->id] = 1;
                if (end->outdegree() == 0 ||
                    end->outedges.front()->head->is_junction()) {
                    break;
                }
                end = end->outedges.front()->head;
                ++extension;
                if (end == it.get()) {
                    is_circular = true;
                    break;
                }
            }

            if (!is_circular && begin == end) {
                continue;
            }
            if (!is_circular && extension < 2 * epsilon + 2) {
                continue;
            }

            if (begin != end) {  // remove nodes near junctions
                for (std::uint32_t i = 0; i < epsilon; ++i) {
                    begin = begin->outedges.front()->head;
                }
                for (std::uint32_t i = 0; i < epsilon; ++i) {
                    end = end->inedges.front()->tail;
                }
            }

            auto unitig = std::make_shared<Node>(begin, end);
            unitigs.emplace_back(unitig);
            unitigs.emplace_back(std::make_shared<Node>(end->pair, begin->pair));
            unitig->pair = unitigs.back().get();
            unitig->pair->pair = unitig.get();

            if (begin != end) {  // connect unitig to graph
                if (begin->indegree()) {
                    marked_edges.emplace(begin->inedges.front()->id);
                    marked_edges.emplace(begin->inedges.front()->pair->id);

                    auto edge = std::make_shared<Edge>(
                            begin->inedges.front()->tail,
                            unitig.get(),
                            begin->inedges.front()->length);
                    unitig_edges.emplace_back(edge);
                    unitig_edges.emplace_back(std::make_shared<Edge>(
                            unitig->pair,
                            begin->inedges.front()->pair->head,
                            begin->inedges.front()->pair->length + unitig->pair->sequence.inflated_len - begin->pair->sequence.inflated_len));  // NOLINT
                    edge->pair = unitig_edges.back().get();
                    edge->pair->pair = edge.get();
                }
                if (end->outdegree()) {
                    marked_edges.emplace(end->outedges.front()->id);
                    marked_edges.emplace(end->outedges.front()->pair->id);

                    auto edge = std::make_shared<Edge>(
                            unitig.get(),
                            end->outedges.front()->head,
                            end->outedges.front()->length + unitig->sequence.inflated_len - end->sequence.inflated_len);  // NOLINT
                    unitig_edges.emplace_back(edge);
                    unitig_edges.emplace_back(std::make_shared<Edge>(
                            end->outedges.front()->pair->tail,
                            unitig->pair,
                            end->outedges.front()->pair->length));
                    edge->pair = unitig_edges.back().get();
                    edge->pair->pair = edge.get();
                }
            }

            auto jt = begin;
            while (true) {
                marked_edges.emplace(jt->outedges.front()->id);
                marked_edges.emplace(jt->outedges.front()->pair->id);

                // update transitive edges
                node_updates[jt->id & ~1UL] = unitig->id;
                unitig->transitive.insert(
                        nodes_[jt->id & ~1UL]->transitive.begin(),
                        nodes_[jt->id & ~1UL]->transitive.end());

                if ((jt = jt->outedges.front()->head) == end) {
                    break;
                }
            }
        }

        nodes_.insert(nodes_.end(), unitigs.begin(), unitigs.end());
        edges_.insert(edges_.end(), unitig_edges.begin(), unitig_edges.end());
        RemoveEdges(marked_edges, true);

        for (const auto& it : nodes_) {  // update transitive edges
            if (it) {
                std::unordered_set<std::uint32_t> valid;
                for (auto jt : it->transitive) {
                    valid.emplace(node_updates[jt] == 0 ? jt : node_updates[jt]);
                }
                it->transitive.swap(valid);
            }
        }

        return unitigs.size() / 2;
    }

    std::uint32_t Unitig::CreateUnitigsAlternate(std::uint32_t epsilon) {
        //std::unordered_set<std::uint32_t> marked_edges;
        std::vector<MarkedEdge> marked_edges;
        std::vector<std::shared_ptr<Node>> unitigs;
        std::vector<std::shared_ptr<Edge>> unitig_edges;
        std::vector<std::uint32_t> node_updates(nodes_alternate_.size(), 0);
        std::vector<char> is_visited(nodes_alternate_.size(), 0);

        for (const auto& it : nodes_alternate_) {
            if (it == nullptr || is_visited[it->id] || it->is_junction()) {
                continue;
            }

            std::uint32_t extension = 1;

            bool is_circular = false;
            auto begin = it.get();
            while (!begin->is_junction()) {  // extend left
                is_visited[begin->id] = 1;
                is_visited[begin->pair->id] = 1;
                if (begin->indegree() == 0 ||
                    begin->inedges.front()->tail->is_junction()) {
                    break;
                }
                begin = begin->inedges.front()->tail;
                ++extension;
                if (begin == it.get()) {
                    is_circular = true;
                    break;
                }
            }

            auto end = it.get();
            while (!end->is_junction()) {  // extend right
                is_visited[end->id] = 1;
                is_visited[end->pair->id] = 1;
                if (end->outdegree() == 0 ||
                    end->outedges.front()->head->is_junction()) {
                    break;
                }
                end = end->outedges.front()->head;
                ++extension;
                if (end == it.get()) {
                    is_circular = true;
                    break;
                }
            }

            if (!is_circular && begin == end) {
                continue;
            }
            if (!is_circular && extension < 2 * epsilon + 2) {
                continue;
            }

            if (begin != end) {  // remove nodes near junctions
                for (std::uint32_t i = 0; i < epsilon; ++i) {
                    begin = begin->outedges.front()->head;
                }
                for (std::uint32_t i = 0; i < epsilon; ++i) {
                    end = end->inedges.front()->tail;
                }
            }

            auto unitig = std::make_shared<Node>(begin, end, Node::num_objects_alternate++);
            unitigs.emplace_back(unitig);
            unitigs.emplace_back(std::make_shared<Node>(end->pair, begin->pair, Node::num_objects_alternate++));
            unitig->pair = unitigs.back().get();
            unitig->pair->pair = unitig.get();

            if (begin != end) {  // connect unitig to graph
                if (begin->indegree()) {
                    marked_edges.emplace_back(begin->inedges.front()->id, 2);
                    marked_edges.emplace_back(begin->inedges.front()->pair->id, 2);

                    auto edge = std::make_shared<Edge>(
                            begin->inedges.front()->tail,
                            unitig.get(),
                            begin->inedges.front()->length, Edge::num_objects_alternate++);
                    unitig_edges.emplace_back(edge);
                    unitig_edges.emplace_back(std::make_shared<Edge>(
                            unitig->pair,
                            begin->inedges.front()->pair->head,
                            begin->inedges.front()->pair->length + unitig->pair->sequence.inflated_len - begin->pair->sequence.inflated_len,
                            Edge::num_objects_alternate++));  // NOLINT
                    edge->pair = unitig_edges.back().get();
                    edge->pair->pair = edge.get();
                }
                if (end->outdegree()) {
                    marked_edges.emplace_back(end->outedges.front()->id, 2);
                    marked_edges.emplace_back(end->outedges.front()->pair->id, 2);

                    auto edge = std::make_shared<Edge>(
                            unitig.get(),
                            end->outedges.front()->head,
                            end->outedges.front()->length + unitig->sequence.inflated_len - end->sequence.inflated_len,
                            Edge::num_objects_alternate++);  // NOLINT
                    unitig_edges.emplace_back(edge);
                    unitig_edges.emplace_back(std::make_shared<Edge>(
                            end->outedges.front()->pair->tail,
                            unitig->pair,
                            end->outedges.front()->pair->length,
                            Edge::num_objects_alternate++));
                    edge->pair = unitig_edges.back().get();
                    edge->pair->pair = edge.get();
                }
            }

            auto jt = begin;
            while (true) {
                marked_edges.emplace_back(jt->outedges.front()->id, 2);
                marked_edges.emplace_back(jt->outedges.front()->pair->id, 2);

                // update transitive edges
                node_updates[jt->id & ~1UL] = unitig->id;
                unitig->transitive.insert(
                        nodes_alternate_[jt->id & ~1UL]->transitive.begin(),
                        nodes_alternate_[jt->id & ~1UL]->transitive.end());

                if ((jt = jt->outedges.front()->head) == end) {
                    break;
                }
            }
        }

        nodes_alternate_.insert(nodes_alternate_.end(), unitigs.begin(), unitigs.end());
        edges_alternate_.insert(edges_alternate_.end(), unitig_edges.begin(), unitig_edges.end());
        RemoveAlternateEdges(marked_edges, true);

        for (const auto& it : nodes_alternate_) {  // update transitive edges
            if (it) {
                std::unordered_set<std::uint32_t> valid;
                for (auto jt : it->transitive) {
                    valid.emplace(node_updates[jt] == 0 ? jt : node_updates[jt]);
                }
                it->transitive.swap(valid);
            }
        }

        return unitigs.size() / 2;
    }



    void Unitig::CreateUnitigGraph(){
        std::unordered_set<std::uint32_t> marked_edges;
        std::vector<std::shared_ptr<Node>> unitigs;
        std::vector<std::shared_ptr<Edge>> unitig_edges;
        std::vector<std::uint32_t> node_updates(nodes_.size(), 0);
        std::vector<char> is_visited(nodes_.size(), 0);


        for (const auto& it : nodes_) {

            if (it == nullptr || is_visited[it->id]) {
                continue;
            }

            std::uint32_t extension = 1;

            bool is_circular = false;
            auto begin = it.get();

            while (true) {  // extend left
                is_visited[begin->id] = 1;
                is_visited[begin->pair->id] = 1;

                if(begin->indegree() > 1 || begin->indegree() == 0){
                    break;
                }

                auto next = begin->inedges.front()->tail;
                if(next->outdegree() != 1){
                    break;
                }

                begin = next;
                ++extension;
                if (begin == it.get()) {
                    is_circular = true;
                    break;
                }
            }


            auto end = it.get();
            while (true) {  // extend right
                is_visited[end->id] = 1;
                is_visited[end->pair->id] = 1;

                if(end->outdegree() > 1 || end->outdegree() == 0){
                    break;
                }

                auto next = end->outedges.front()->head;

                if(next->indegree() != 1){
                    break;
                }
                end = next;
                ++extension;
                if (end == it.get()) {
                    is_circular = true;
                    break;
                }
            }

            auto unitig = std::make_shared<Node>(begin, end, true);
            unitigs.emplace_back(unitig);
            unitigs.emplace_back(std::make_shared<Node>(end->pair, begin->pair, true));
            unitig->pair = unitigs.back().get();
            unitig->pair->pair = unitig.get();
        }

        auto GetUnitigFromNode = [&] (Node* query_node) -> Node* {
            for(auto& unitig : unitigs){
                for(auto unitig_node : unitig->unitig_nodes){
                    if(unitig_node->id == query_node->id){
                        return unitig.get();
                    };
                };
            };
            return nullptr;
        };

        std::uint32_t counter = 0;
        for(auto& unitig : unitigs){
            if (!unitig->is_circular) {  // connect unitig to graph
                counter++;
                for(auto& inedge : unitig->front_inedges){
                    marked_edges.emplace(inedge->id);
                    marked_edges.emplace(inedge->pair->id);

                    auto tail_unitig = GetUnitigFromNode(inedge->tail);

                    auto edge = std::make_shared<Edge>(
                            tail_unitig,
                            unitig.get(),
                            tail_unitig->sequence.inflated_len - inedge->tail->sequence.inflated_len + inedge->length
                    );
                    unitig_edges.emplace_back(edge);

                    auto rc_head_unitig = GetUnitigFromNode(inedge->pair->head);
                    auto rc_edge = std::make_shared<Edge>(
                            unitig->pair,
                            rc_head_unitig,
                            unitig->pair->sequence.inflated_len - inedge->pair->head->sequence.inflated_len + inedge->length
                    );
                    unitig_edges.emplace_back(rc_edge);
                    edge->pair = unitig_edges.back().get();
                    edge->pair->pair = edge.get();
                };

                for(auto& outedge : unitig->back_outedges){
                    marked_edges.emplace(outedge->id);
                    marked_edges.emplace(outedge->pair->id);

                    auto head_unitig = GetUnitigFromNode(outedge->head);
                    auto edge = std::make_shared<Edge>(
                            unitig.get(),
                            head_unitig,
                            unitig->sequence.inflated_len - outedge->tail->sequence.inflated_len + outedge->length
                    );
                    unitig_edges.emplace_back(edge);

                    auto rc_tail_unitig = GetUnitigFromNode(outedge->pair->head);
                    auto rc_edge = std::make_shared<Edge>(
                            rc_tail_unitig,
                            unitig->pair,
                            rc_tail_unitig->sequence.inflated_len - outedge->pair->head->sequence.inflated_len + outedge->length
                    );
                    unitig_edges.emplace_back(rc_edge);
                    edge->pair = unitig_edges.back().get();
                    edge->pair->pair = edge.get();
                };
            };

        }
        unitig_nodes_.insert(unitig_nodes_.end(), unitigs.begin(), unitigs.end());
        unitig_edges_.insert(unitig_edges_.end(), unitig_edges.begin(), unitig_edges.end());

        std::string output_path = "unitig_graph.gfa";
        PrintUnitigGfa(output_path);
    }

    void Unitig::ResolveGraphWithUl(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &ul_reads){
        //get unitig sequences
        //
    }

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> Unitig::GetUnitigPairs(
            bool drop_unpolished) {

        //CreateUnitigs();

        biosoup::NucleicAcid::num_objects = 0;

        std::vector<std::unique_ptr<biosoup::NucleicAcid>> dst;
        for (const auto& it : nodes_alternate_) {
            if (it == nullptr || it->is_rc() || !it->is_unitig) {
                continue;
            }
            if (drop_unpolished && !it->is_polished) {
                continue;
            }

            std::string name = it->sequence.name +
                               " LN:i:" + std::to_string(it->sequence.inflated_len) +
                               " RC:i:" + std::to_string(it->count) +
                               " XO:i:" + std::to_string(it->is_circular);

            dst.emplace_back(new biosoup::NucleicAcid(
                    name,
                    it->sequence.InflateData()));
        }

        return dst;
    }

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> Unitig::GetAssembledData(bool primary) {

        biosoup::NucleicAcid::num_objects = 0;
        auto node_list = nodes_;
        if (!primary)
            node_list = nodes_alternate_;

        if (node_list.size() == 0)
            return std::vector<std::unique_ptr<biosoup::NucleicAcid>> {};

        std::vector<std::unique_ptr<biosoup::NucleicAcid>> dst;
        for (const auto& it : node_list) {
            if (it == nullptr || it->is_rc()) {
                continue;
            }

//      if (!it->is_unitig) {
//        continue;
//      }

            std::string name = it->sequence.name +
                               " LN:i:" + std::to_string(it->sequence.inflated_len) +
                               " RC:i:" + std::to_string(it->count) +
                               " XO:i:" + std::to_string(it->is_circular);

            dst.emplace_back(new biosoup::NucleicAcid(
                    name,
                    it->sequence.InflateData()));
        }

        return dst;
    }

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> Unitig::GetUnitigs(
            bool drop_unpolished) {

        CreateUnitigs();

        biosoup::NucleicAcid::num_objects = 0;

        std::vector<std::unique_ptr<biosoup::NucleicAcid>> dst;
        for (const auto &it: nodes_) {
            if (it == nullptr || it->is_rc() || !it->is_unitig) {
                continue;
            }
            if (drop_unpolished && !it->is_polished) {
                continue;
            }

            std::string name = it->sequence.name +
                               " LN:i:" + std::to_string(it->sequence.inflated_len) +
                               " RC:i:" + std::to_string(it->count) +
                               " XO:i:" + std::to_string(it->is_circular);

            dst.emplace_back(new biosoup::NucleicAcid(
                    name,
                    it->sequence.InflateData()));
        }

        return dst;
    }

}