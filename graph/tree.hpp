#ifndef _TREE_H_
#define _TREE_H_
#include <vector>
#include <algorithm>
#include "edge.hpp"
#include "tree_algorithm.hpp"

// 双方向に辺を張らない場合や多重辺がある場合壊れる
template<typename edge>
struct tree{
  using __edge = edge;
  using weight = typename edge::weight;
  template<typename T>
  using vec = std::vector<T>;
  using graph = vec<vec<edge>>;
  using simple_tree = tree<simple_edge<int>>;
public:
  graph g;
  int n, r;
  vec<int> subsize, depth, parent;
  tree_algorithm::hld<edge> *hld_p;
  tree_algorithm::dfs_order<edge> *dfs_p;
  tree_algorithm::bfs_order<edge> *bfs_p;
  tree_algorithm::weighted_dist<tree<edge>> *dist_p;
  tree(int n, int r = 0): g(n), n(n), r(r), hld_p(nullptr), dfs_p(nullptr), bfs_p(nullptr), dist_p(nullptr){}
  tree(graph &g, int r = 0): g(g), n(g.size()), r(r), hld_p(nullptr), dfs_p(nullptr), bfs_p(nullptr), dist_p(nullptr){}
  void add_edge(int a, edge e){
    assert(0 <= a && a < n);
    g[a].push_back(e);
  }
  void simple_dfs(){
    auto [s, d, p] = tree_algorithm::simple_dfs(g, r);
    subsize = s, depth = d, parent = p;
  }
  void hld_build(){
    hld_p = new tree_algorithm::hld<edge>(g, r);
    subsize = hld_p->subsize, depth = hld_p->depth, parent = hld_p->parent;
  }
  void dfs_build(){
    dfs_p = new tree_algorithm::dfs_order(g, r);
    subsize = dfs_p->subsize, depth = dfs_p->depth, parent = dfs_p->parent;
  }
  void bfs_build(){
    bfs_p = new tree_algorithm::bfs_order(g, r);
    parent = bfs_p->parent;
  }
  void dist_build(){
    dist_p = new tree_algorithm::weighted_dist(*this);
  }
  // O(logN)
  int lca(int u, int v){
    if(!hld_p) hld_build();
    return hld_p->lca(u, v);
  }
  // k個上の祖先, k = 0なら自身, k > depth[u]なら-1
  // O(logN)
  int la(int u, int k){
    if(!hld_p) hld_build();
    return hld_p->la(u, k);
  }
  // O(1)
  int dep(int v){
    if(depth.empty()) dfs_build();
    return depth[v];
  }
  // O(1)
  int par(int v){
    if(parent.empty()) dfs_build();
    return parent[v];
  }
  // O(1)
  int size(int v){
    if(subsize.empty()) dfs_build();
    return subsize[v];
  }
  // 重みなしの距離
  // O(logN)
  int dist_unweighted(int u, int v){
    if(!hld_p) hld_build();
    return hld_p->dist(u, v);
  }
  //　重み付き距離
  // O(logN)
  weight dist_weighted(int u, int v){
    if(!hld_p) hld_build();
    if(!dist_p) dist_build();
    return dist_p->dist(u, v);
  }
  // 辺(u, v)(双方向)にw足す(ない場合は何もしない)
  void add_weight(int u, int v, weight w){
    if(!hld_p) hld_build();
    if(!dist_p) dist_build();
    if(parent[u] != v && parent[v] != u) return;
    if(parent[u] == v) std::swap(u, v); // u->vにする
    int uv = find_edge(u, v), vu = find_edge(v, u);
    g[u][uv].w += w, g[v][vu].w += w;
    dist_p->add_weight(u, v, w);
  }
  // u -> vのパスのk番目, k = 0ならu, k > パス長なら-1
  // O(logN)
  int kth_vertex_on_path(int u, int v, int k){
    if(!hld_p) hld_build();
    return hld_p->kth_vertex_on_path(u, v, k);
  }
  // vがuの部分木に含まれるか(u自身も部分木)
  // O(1)(dfs)またはO(logN)(hld)
  bool is_contained_subtree(int u, int v){
    if(!dfs_p && !hld_p) dfs_build();
    if(dfs_p) return dfs_p->is_contained_subtree(u, v);
    else return hld_p->is_contained_subtree(u, v);
  }
  // wがlca(u, v) or wがuかvのどちらか片方のみ部分木として含む
  bool is_contained_path(int u, int v, int w){
    if(!hld_p) hld_build();
    return hld_p->is_contained_path(u, v, w);
  }
  // 辺{u, v}がある場合, g[u]のインデックスを返す, 無い場合は-1
  int find_index(int u, int v){
    if(!bfs_p) bfs_build();
    int k = bfs_p->find_index(u, v);
    if(k >= g[u].size()) return -1;
    return k;
  }
  // 辺{u, v}がある場合辺を返す, 無い場合エラー
  edge find_edge(int u, int v){
    int k = find_index(u, v);
    assert(k != -1);
    return g[u][k];
  }
  // s->tパスの辺
  vec<edge> get_path(int s, int t){
    int l = lca(s, t);
    vec<edge> L, R;
    while(s != l){
      int p = parent[s];
      L.push_back(find_edge(s, p));
      s = p;
    }
    while(t != l){
      int p = parent[t];
      R.push_back(find_edge(p, t));
      t = parent[t];
    }
    std::reverse(R.begin(), R.end());
    L.insert(L.end(), R.begin(), R.end());
    return L;
  }
   // hldでの頂点vの探索順
  int index_hld(int v){
    if(!hld_p) hld_build();
    return hld_p->index_vertex(v);
  }
  // hldの探索順でu-vパスに相当するO(logN)個の区間
  std::vector<std::pair<int, int>> unordered_path_to_segment(int u, int v, bool is_edge = false){
    if(!hld_p) hld_build();
    return hld_p->unordered_path(u, v, is_edge);
  }
  /*
  // hldの探索順でu-vパスに相当するO(logN)個の区間
  std::vector<std::pair<int, int>> ordered_path_to_segment(int u, int v, bool is_edge = false){
    if(!hld_p) hld_build();
    return hld_p->unordered_path(u, v, is_edge);
  }
  */
  // f : 区間[l, r)に対する更新
  template<typename F>
  void update_path(int u, int v, F f, bool is_edge = false){
    if(!hld_p) hld_build();
    hld_p->update_path(u, v, f, is_edge);
  }
  // f : 区間[l, r)に対する求値, g : 値のマージ, z : 単位元
  template<typename F, typename G, typename Val>
  Val query_path(int u, int v, F f, G g, Val z, bool is_edge = false){
    if(!hld_p) hld_build();
    return hld_p->query_path(u, v, f, g, z, is_edge);
  }
  // f : 区間[l, r)に対する求値, g : 値のマージ, h : パスの反転, z : 単位元
  template<typename F, typename G, typename H, typename Val>
  Val query_path_flip(int u, int v, F f, G g, H h, Val z, bool is_edge = false){
    if(!hld_p) hld_build();
    return hld_p->query_path_flip(u, v, f, g, h, z, is_edge);
  }
  // dfsのpreorderで頂点vの探索順
  int index_dfs_preorder(int v){
    if(!dfs_p) dfs_build();
    return dfs_p->in_pre[v];
  }
  // dfsのpreorderでvの部分木に相当する区間[l, r)
  std::pair<int, int> subtree_to_segment(int v){
    if(!dfs_p) dfs_build();
    return dfs_p->subtree_to_segment(v);
  }
  // f : 区間[l, r)に対する更新
  template<typename F>
  void update_subtree(int u, F f, bool is_edge = false){
    if(!dfs_p) dfs_build();
    dfs_p->update_subtree(u, f, is_edge);
  }
  // f : 区間[l, r)に対する求値
  template<typename F>
  auto query_subtree(int u, F f, bool is_edge = false){
    if(!dfs_p) dfs_build();
    return dfs_p->query_subtree(u, f, is_edge);
  }
  // {重み付き距離, s, t}
  std::tuple<weight, int, int> diameter(){
    return tree_algorithm::tree_diameter(g);
  }
  static simple_tree make_simple_tree(const vec<vec<int>> &G, int root){
    int n = G.size();
    vec<vec<simple_edge<int>>> G2(n);
    for(int i = 0; i < n; i++) for(int j : G[i]) G2[i].push_back(simple_edge<int>(i, j));
    return simple_tree(G2, root);
  }
  // {頂点の集合, 木}
  std::pair<vec<int>, simple_tree> lca_tree(vec<int> v){
    if(!dfs_p) dfs_build();
    auto __lca = [&](int u, int v){return lca(u, v);};
    auto __dfs = [&](int v){return dfs_p->in_pre[v];};
    auto __dep = [&](int v){return dep(v);};
    auto [V, G, root] = tree_algorithm::lca_tree(v, __lca, __dfs, __dep);
    return {V, make_simple_tree(G, root)};
  }
  simple_tree centroid_decomposition(){
    auto [G, root, size_i, dep_i, par_i] = tree_algorithm::centroid_decomposition<edge>(g);
    simple_tree ret = make_simple_tree(G, root);
    ret.subsize = size_i;
    ret.depth = dep_i;
    ret.parent = par_i;
    return ret;
  }
  vec<edge> &operator [](int i){return g[i];}
};
using simple_tree = tree<simple_edge<int>>;

#endif