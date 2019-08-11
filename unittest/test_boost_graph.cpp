//
// Created by xiamr on 8/11/19.
//

#include "std.hpp"
#include <gmock/gmock.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/labeled_graph.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/reverse_graph.hpp>

using namespace std;
using namespace testing;

struct VertexNode {
    typedef boost::vertex_property_tag kind;
};


class MyVisitor : public boost::default_dfs_visitor {

public:
    MyVisitor(std::unordered_set<int> &all_connected_atoms) : all_connected_atoms(all_connected_atoms) {};

    template<typename Vertex, typename Graph>
    void discover_vertex(Vertex v, const Graph &g) const {
        all_connected_atoms.emplace(boost::get(VertexNode(), g, v));
    }

private:

    std::unordered_set<int> &all_connected_atoms;
};


TEST(BoostGraph, Start) {

    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::property<VertexNode, int>> graph_t;

    graph_t g;

    auto v1 = boost::add_vertex(10, g);
    auto v2 = boost::add_vertex(13, g);
    auto v3 = boost::add_vertex(14, g);
    auto v4 = boost::add_vertex(15, g);


    boost::add_edge(v1, v2, g);
    boost::add_edge(v2, v1, g);
    boost::add_edge(v2, v4, g);

    ASSERT_THAT(boost::num_vertices(g), Eq(4));
    ASSERT_THAT(boost::num_edges(g), Eq(2));


    std::unordered_set<int> all_connected_atoms;
    MyVisitor vis(all_connected_atoms);

    boost::depth_first_visit(g, v1, vis, boost::make_vector_property_map<boost::default_color_type>(
            boost::get(boost::vertex_index, g)));

    ASSERT_THAT(all_connected_atoms, UnorderedElementsAre(10, 13, 15));


}
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/depth_first_search.hpp>
//#include <iostream>
//using namespace std;
//
////======================================================================================================
//// information representing each vertex
//struct Vertex {
//
//    int id = 0;
//    int parent_id = -1;
//    int l_child = -1, r_child = -1;
//
//    Vertex(int id = -1) : id(id) {}
//};
//
////======================================================================================================
//// information representing each weight
//// it carries the boundary length and also the distance
//struct Edge {
//
//    // distance
//    float boundary_length = 0;
//    float weight = 1;
//    // float L, a, b = 1;
//
//    Edge(float boundary_length = 1) : boundary_length(boundary_length) {}
//};
//
//typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, Vertex, Edge> Graph;
//
//class MyVisitor : public boost::default_dfs_visitor {
//  public:
//    MyVisitor() : vv(new std::vector<int>()) {}
//
//    void discover_vertex(int v, const Graph &g) const {
//        vv->push_back(g[v].id);
//        return;
//    }
//
//    std::vector<int> &GetVector() const { return *vv; }
//
//  private:
//    boost::shared_ptr<std::vector<int> > vv;
//};
//
//int main() {
//
//    Graph g;
//    boost::add_edge(0, 1, g);
//    boost::add_edge(0, 2, g);
//    boost::add_edge(1, 2, g);
//    boost::add_edge(1, 3, g);
//
//    boost::add_edge(5, 6, g);
//    boost::add_edge(5, 8, g);
//
//    for (auto v : boost::make_iterator_range(boost::vertices(g)))
//        g[v] = Vertex(v);
//
//    auto indexmap = boost::get(boost::vertex_index, g);
//    auto colormap = boost::make_vector_property_map<boost::default_color_type>(indexmap);
//
//    MyVisitor vis;
//    boost::depth_first_search(g, vis, colormap, 1);
//
//    std::vector<int> vctr = vis.GetVector();
//
//    for(auto id : vctr)
//        std::cout << id << " ";
//}




