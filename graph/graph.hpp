#ifndef _GRAPH_H_
#define _GRAPH_H_
#include <algorithm>
#include <vector>
#include "edge.hpp"
#include "graph_algorithm.hpp"

template<typename edge>
struct general_graph{
  using weight = typename edge::weight;
  template<typename T>
  using vec = std::vector<T>;
  using graph = vec<vec<edge>>;
  using simple_graph = general_graph<simple_edge<int>>;
  int n;
  graph g;
  general_graph(int n): n(n), g(n){}
  general_graph(const vec<vec<edge>> &G): n(G.size()), g(G){}

  void add_edge(int a, edge e){
    g[a].push_back(e);
  }
  graph_algorithm::bfs_shortest_path<edge> bfs_shortest_path(){
    return graph_algorithm::bfs_shortest_path<edge>(g);
  }
  graph_algorithm::zero_one_bfs_shortest_path<edge> zero_one_bfs_shortest_path(){
    return graph_algorithm::zero_one_bfs_shortest_path<edge>(g);
  }
  graph_algorithm::dijkstra<edge> dijkstra(){
    return graph_algorithm::dijkstra<edge>(g);
  }
  graph_algorithm::bellman_ford<edge> bellman_ford(){
    return graph_algorithm::bellman_ford<edge>(g);
  }
  graph_algorithm::warshall_floyd<edge> warshall_floyd(){
    return graph_algorithm::warshall_floyd<edge>(g);
  }
  static simple_graph make_simple_graph(const vec<vec<int>> &G){
    int n = G.size();
    vec<vec<simple_edge<int>>> G2(n);
    for(int i = 0; i < n; i++) for(int j : G[i]) G2[i].push_back(simple_edge<int>(i, j));
    return simple_graph(G2);
  }
  // {連結成分, グラフ} 有向サイクル
  std::pair<vec<int>, simple_graph> scc(){
    vec<int> cmp = graph_algorithm::scc(g).first;
    int m = *std::max_element(cmp.begin(), cmp.end()) + 1;
    vec<vec<int>> G(m);
    for(int i = 0; i < n; i++){
      for(auto &e : g[i]){
        if(cmp[i] != cmp[e.t]){
          G[cmp[i]].push_back(cmp[e.t]);
        }
      }
    }
    for(int i = 0; i < m; i++){
      std::sort(G[i].begin(), G[i].end());
      G[i].erase(std::unique(G[i].begin(), G[i].end()), G[i].end());
    }
    return {cmp, make_simple_graph(G)};
  }
  // {連結成分, グラフ} 1つの辺が消えても連結, 森か木になる
  std::pair<vec<int>, vec<vec<int>>> two_edge_connected(){
    vec<int> cmp = graph_algorithm::two_edge_connected(g).first;
    int m = *std::max_element(cmp.begin(), cmp.end()) + 1;
    vec<vec<int>> G(m);
    for(int i = 0; i < n; i++){
      for(auto &e : g[i]){
        if(cmp[i] != cmp[e.t]){
          G[cmp[i]].push_back(cmp[e.t]);
        }
      }
    }
    return {cmp, G};
  }
  // {関節点フラグ, 各連結成分が含む頂点(関節点は複数の連結成分に含まれる))} 1つの頂点が消えても連結
  std::pair<vec<bool>, vec<vec<int>>> bcc(){
    return graph_algorithm::bcc<edge>(g);
  }
  // {関節点フラグ, cmp, 森} 1つの頂点が消えても連結, 森か木になる
  std::tuple<vec<bool>, vec<int>, vec<vec<edge>>> block_cut_tree(){
    return graph_algorithm::block_cut_tree<edge>(g);
  }
  // 閉路が存在するなら空のvector
  vec<int> topological_sort(){
    return graph_algorithm::topological_sort(g);
  }
  // 最小全域木, グラフが連結ならsの値は関係ない
  vec<edge> undirected_mst(int s = 0){
    return graph_algorithm::undirected_mst<edge>(g, s);
  }
  // rを根とするbfs木O(V + E)
   vec<vec<edge>> bfs_tree(int r){
    return graph_algorithm::bfs_tree(g, r);
  }
  // rを根とするbfs木, rからの最短経路 O((V + E)logV)
  std::pair<vec<vec<edge>>, vec<weight>> bfs_tree_shortest(int r){
    return graph_algorithm::bfs_tree_shortest(g, r);
  }
  // g[i]の辺を{同じcmpへの辺, 異なるcmpへの辺}に並び替える, O(V + E)
  void cmp_edge_arrange(const vec<int> &cmp){
    graph_algorithm::cmp_edge_arrange(cmp, g);
  }
  vec<edge> &operator [](int i){return g[i];}
};

using simple_graph = general_graph<simple_edge<int>>;
template<typename T>
using weighted_graph = general_graph<weighted_edge<T>>;
template<typename T>
using labeled_graph =  general_graph<labeled_edge<T>>;
template<typename T>
using weighted_labeled_graph =  general_graph<weighted_labeled_edge<T>>;

#endif
