#ifndef _FAST_LCA_H_
#define _FAST_LCA_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include "../data_structure/range_query/sparse_table.hpp"

struct fast_lca{
private:
  using graph = std::vector<std::vector<int>>;
  int N;
  std::vector<int> begin, end, begin_compressed, tour_compressed, dep, id2order, order2id;
  rmq<int> table;
  void init_dfs(int v, int p, const graph &g, int &a, int &b, int &c){
    id2order[v] = b, order2id[b++] = v, begin[v] = end[v] = c++;
    for(int t : g[v]){
      if(t == p) continue;
      dep[t] = dep[v] + 1;
      init_dfs(t, v, g, a, b, c);
      end[v] = c++;
      tour_compressed.push_back(id2order[v]);
      if(begin_compressed[v] == -1) begin_compressed[v] = a;
      a++;
    }
    if(begin_compressed[v] == -1) begin_compressed[v] = a;
  }
public:
  fast_lca(){}
  // init O(NlogN)
  fast_lca(const graph &g, int root): N(g.size()), begin(N), end(N),
  begin_compressed(N, -1), dep(N), id2order(N), order2id(N){
    dep[root] = 0;
    int a = 0, b = 0, c = 0;
    init_dfs(root, -1, g, a, b, c);
    table = rmq<int>(tour_compressed);
  }
  template<typename edge>
  fast_lca(const std::vector<std::vector<edge>> &g, int root): N(g.size()), begin(N), end(N), 
  begin_compressed(N, -1), dep(N), id2order(N), order2id(N){
    std::vector<std::vector<int>> g2(N);
    for(int i = 0; i < N; i++){
      for(auto &e : g[i]){
        g2[i].push_back(e.t);
      }
    }
    dep[root] = 0;
    int a = 0, b = 0, c = 0;
    init_dfs(root, -1, g2, a, b, c);
    table = rmq<int>(tour_compressed);
  }
  // preorderでvを訪れる順番
  inline int pre_order(int v)const{
    return id2order[v];
  }
  // preorderでk番目の頂点
  inline int pre_kth(int k)const{
    return order2id[k];
  }
  // euler_tourでvが最初に現れる順番
  inline int et_begin(int v)const{
    return begin[v];
  }
  // eouler_tourでvが最後に現れる位置
  inline int et_end(int v)const{
    return end[v];
  }
  // pがcの祖先か(p == cの場合も含む)
  inline bool is_ancesstor(int p, int c)const{
    return begin[p] <= begin[c] && end[c] <= end[p];
  }
  inline bool is_leaf(int v)const{
    return begin[v] == end[v];
  }
  // O(1) per query
  inline int lca(int u, int v){
    if(begin_compressed[u] == begin_compressed[v]) return (dep[u] > dep[v] ? v : u);
    if(begin_compressed[u] > begin_compressed[v]) std::swap(u, v);
    return order2id[table.query(begin_compressed[u], begin_compressed[v] + 1)];
  }
  inline int depth(int v)const{
    return dep[v];
  }
  inline int dist(int u, int v){
    return dep[u] + dep[v] - 2 * dep[lca(u, v)];
  }
};
#endif