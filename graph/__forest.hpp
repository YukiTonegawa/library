#ifndef _FOREST_H_
#define _FOREST_H_
#include <vector>
#include "edge.hpp"
#include "forest_algorithm.hpp"

// 双方向に辺を張らない場合や多重辺がある場合壊れる
template<typename edge>
struct forest{
  using _edge = edge;
  using weight = typename edge::weight;
  template<typename T>
  using vec = std::vector<T>;
  using graph = vec<vec<edge>>;
  using simple_forest = forest<simple_edge<int>>;
public:
  graph g;
  int n;
  vec<int> subsize, depth, parent, root;
  forest_algorithm::hld<forest<edge>> *hld_p;
  forest_algorithm::dfs_order<forest<edge>> *dfs_p;
  forest_algorithm::bfs_order<forest<edge>> *bfs_p;

  forest(int n): g(n), n(n), hld_p(nullptr), dfs_p(nullptr), bfs_p(nullptr){}
  forest(graph &g): g(g), n(g.size()), hld_p(nullptr), dfs_p(nullptr), bfs_p(nullptr){}
  static simple_forest make_simple_forest(const vec<vec<int>> &G){
    int n = G.size();
    vec<vec<simple_edge<int>>> G2(n);
    for(int i = 0; i < n; i++) for(int j : G[i]) G2[i].push_back(simple_edge<int>(i, j));
    return simple_forest(G2);
  }
  void add_edge(int a, edge e){
    assert(0 <= a && a < n);
    g[a].push_back(e);
  }
  // rを含む連結成分をrを根としてdfs
  void simple_dfs(int r){
    forest_algorithm::simple_dfs(*this, r);
  }
  // 全ての連結成分にdfs
  void simple_dfs_all(){
    for(int i = 0; i < n; i++) forest_algorithm::simple_dfs(*this, i);
  }
  // rを含む連結成分をrを根として初期化　, すでに連結成分が初期化されている場合は何もしない
  void hld_build(int r){
    if(!hld_p) hld_p = new forest_algorithm::hld<forest<edge>>(*this);
    hld_p->build(r);
  }
  void dfs_build(int r){
    if(!dfs_p) dfs_p = new forest_algorithm::dfs_order<forest<edge>>(*this);
    dfs_p->build(r);
  }
  void bfs_build(int r){
    if(!bfs_p) bfs_p = new forest_algorithm::bfs_order<forest<edge>>(*this);
    bfs_p->build(r);
  }
  // O(logN), 異なる連結成分の場合は-1
  int lca(int u, int v){
    assert(hld_p);
    return hld_p->lca(u, v);
  }
  // k個上の祖先, k = 0なら自身, k > depth[u]なら-1
  // O(logN)
  int la(int u, int k){
    assert(hld_p);
    return hld_p->la(u, k);
  }
  // O(1)
  int dep(int v){
    assert(!depth.empty());
    return depth[v];
  }
  // O(1)
  int par(int v){
    assert(!parent.empty());
    return parent[v];
  }
  // O(1)
  int size(int v){
    assert(!subsize.empty());
    return subsize[v];
  }
  int find_root(int v){
    assert(!root.empty());
    return root[v];
  }
  // 重みなしの距離
  // O(logN)
  int dist_unweighted(int u, int v){
    assert(hld_p);
    return hld_p->dist(u, v);
  }
  
  //　重み付き距離
  //
  //  


  // u -> vのパスのk番目, k = 0ならu, k > パス長なら-1
  // O(logN)
  int kth_vertex_on_path(int u, int v, int k){
    assert(hld_p);
    return hld_p->kth_vertex_on_path(u, v, k);
  }
  // vがuの部分木に含まれるか(u自身も部分木)
  // O(1)(dfs)またはO(logN)(hld)
  bool is_contained_subtree(int u, int v){
    assert(dfs_p);
    return dfs_p->is_contained_subtree(u, v);
  }
  // u->vパス(最短経路)にwが含まれるか(端点も含む)
  // vがuの部分木 -> wがuの部分木 && vがwの部分木
  // それ以外 -> wがuかvのどちらかを部分木として含む
  // O(1)
  bool is_contained_path(int u, int v, int w){
    assert(dfs_p);
    return dfs_p->is_contained_path(u, v, w);
  }
  // 辺{u, v}がある場合, g[u]のインデックスを返す, 無い場合は-1
  int find_index(int u, int v){
    assert(bfs_p);
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
  // 同じ連結成分でないなら空のvector
  vec<edge> get_path(int s, int t){
    int l = lca(s, t);
    if(l == -1) return {};
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
    assert(hld_p);
    return hld_p->index_vertex(v);
  }
  // hldの探索順でu-vパスに相当するO(logN)個の区間
  std::vector<std::pair<int, int>> path_to_segment(int u, int v){
    assert(hld_p);
    return hld_p->ordered_path();
  }
  // f : 区間[l, r)に対する更新
  template<typename F>
  void update_path(int u, int v, F f, bool is_edge = false){
    assert(hld_p);
    hld_p->update_path(u, v, f, is_edge);
  }
  // f : 区間[l, r)に対する求値, g : 値のマージ, z : 単位元
  template<typename F, typename G, typename Val>
  Val query_path(int u, int v, F f, G g, Val z, bool is_edge = false){
    assert(hld_p);
    return hld_p->query_path(u, v, f, g, z, is_edge);
  }
  // f : 区間[l, r)に対する求値, g : 値のマージ, h : パスの反転, z : 単位元
  template<typename F, typename G, typename H, typename Val>
  Val query_path_flip(int u, int v, F f, G g, H h, Val z, bool is_edge = false){
    assert(hld_p);
    return hld_p->query_path_flip(u, v, f, g, h, z, is_edge);
  }
  // dfsのpreorderで頂点vの探索順
  int index_dfs_preorder(int v){
    assert(dfs_p);
    return dfs_p->in_pre[v];
  }
  // dfsのpreorderでvの部分木に相当する区間[l, r)
  std::pair<int, int> subtree_to_segment(int v){
    assert(dfs_p);
    return dfs_p->subtree_to_segment(v);
  }
  // f : 区間[l, r)に対する更新
  template<typename F>
  void update_subtree(int u, F f, bool is_edge = false){
    assert(dfs_p);
    dfs_p->update_subtree(u, f, is_edge);
  }
  // f : 区間[l, r)に対する求値
  template<typename F>
  auto query_subtree(int u, F f, bool is_edge = false){
    assert(dfs_p);
    return dfs_p->query_subtree(u, f, is_edge);
  }
  // rを含む連結成分の直径を求める{重み付き距離, s, t}
  std::tuple<weight, int, int> diameter(int r){
    return forest_algorithm::diameter(*this, r);
  }
  simple_forest make_simple_tree(const vec<vec<int>> &G){
    int n = G.size();
    vec<vec<simple_edge<int>>> G2(n);
    for(int i = 0; i < n; i++) for(int j : G[i]) G2[i].push_back(simple_edge<int>(i, j));
    return simple_forest(G2);
  }
  // {頂点の集合, 森}
  std::pair<vec<int>, simple_forest> lca_tree(vec<int> v){
    assert(dfs_p && hld_p);
    auto __lca = [&](int u, int v){
      return lca(u, v);
    };
    auto __dfs = [&](int v){
      return dfs_p->in_pre[v];
    };
    auto __dep = [&](int v){
      return dep(v);
    };
    
    //auto [V, G, root] = tree_algorithm::lca_tree(v, __lca, __dfs, __dep);
    //return {V, make_simple_forest(G)};
  }
  /*
  simple_tree centroid_decomposition(){
    auto [G, root, size_i, dep_i, par_i] = tree_algorithm::centroid_decomposition<edge>(g);
    simple_tree ret = make_simple_tree(G, root);
    ret.subsize = size_i;
    ret.depth = dep_i;
    ret.parent = par_i;
    return ret;
  }
  */
  vec<edge> &operator [](int i){return g[i];}
};
using simple_forest = forest<simple_edge<int>>;

#endif