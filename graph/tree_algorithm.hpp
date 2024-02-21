#ifndef _TREE_ALGORITHM_H_
#define _TREE_ALGORITHM_H_
#include <vector>
#include <queue>
#include <stack>
#include "../data_structure/segment_tree/binary_indexed_tree.hpp"
#include "../data_structure/segment_tree/segment_tree.hpp"

namespace tree_algorithm{
  template<typename T>
  using vec = std::vector<T>;
  // 部分木のサイズ, 深さ, 親
  template<typename edge>
  std::tuple<vec<int>, vec<int>, vec<int>> simple_dfs(const vec<vec<edge>> &g, int root){
    int n = g.size();
    vec<int> sz(n, 1), de(n), pa(n);
    auto simple_dfs_f = [&](auto simple_dfs_f, int cur, int par, int dep)->void{
      pa[cur] = par;
      de[cur] = dep;
      for(edge &e : g[cur]){
        if(e.t == par) continue;
        simple_dfs_f(simple_dfs_f, e.t, cur, dep + 1);
        sz[cur] += sz[e.t];
      }
    };
    simple_dfs_f(simple_dfs_f, root, -1, 0);
    return {sz, de, pa};
  }
  template<typename edge>
  struct hld{
    vec<int> subsize, depth, parent, in, out, head, rev, heavy;
    hld(vec<vec<edge>> &g, int root){
      build(g, root);
    }
    void dfs_sz(int cur, int par, int dep, vec<vec<edge>> &g){
      depth[cur] = dep;
      parent[cur] = par;
      subsize[cur] = 1;
      for(int i = 0; i < g[cur].size(); i++){
        edge &e = g[cur][i];
        if(e.t == par) continue;
        dfs_sz(e.t, cur, dep + 1, g);
        subsize[cur] += subsize[e.t];
        if(heavy[cur] == -1 || subsize[g[cur][heavy[cur]].t] < subsize[e.t]){
          heavy[cur] = i;
        }
      }
    }
    void dfs_hld(int cur, int par, int &times, vec<vec<edge>> &g){
      in[cur] = times++;
      rev[in[cur]] = cur;
      int h = heavy[cur];
      if(h != -1){
        head[g[cur][h].t] = head[cur];
        dfs_hld(g[cur][h].t, cur, times, g);
      }
      for(int i = 0; i < h; i++){
        if(g[cur][i].t == par) continue;
        head[g[cur][i].t] = g[cur][i].t;
        dfs_hld(g[cur][i].t, cur, times, g);
      }
      for(int i = h + 1; i < g[cur].size(); i++){
        if(g[cur][i].t == par) continue;
        head[g[cur][i].t] = g[cur][i].t;
        dfs_hld(g[cur][i].t, cur, times, g);
      }
      out[cur] = times;
    }
    void build(vec<vec<edge>> &g, int root){
      int n = g.size();
      subsize.resize(n), depth.resize(n), parent.resize(n);
      in.resize(n), out.resize(n), head.resize(n), rev.resize(n), heavy.resize(n, -1);
      dfs_sz(root, -1, 0, g);
      head[root] = root;
      int t = 0;
      dfs_hld(root, -1, t, g);
    }
    // k個上の祖先, k = 0なら自身, k > depth[v]なら-1
    int la(int v, int k){
      if(depth[v] < k) return -1;
      while(true){
        int u = head[v];
        if(in[v] - k >= in[u]) return rev[in[v] - k];
        k -= in[v] - in[u] + 1;
        v = parent[u];
      }
    }
    int lca(int u, int v){
      for(;; v = parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v);
        if(head[u] == head[v]) return u;
      }
    }
    // 重みなしの距離
    int dist(int u, int v){
      int l = lca(u, v);
      return depth[u] + depth[v] - 2 * depth[l];
    }
    // vがuの部分木に含まれるか(u自身も部分木)
    bool is_contained_subtree(int u, int v){
      if(depth[u] > depth[v]) return false;
      return u == la(v, depth[v] - depth[u]);
    }
    // u->vパス(最短経路)にwが含まれるか(端点も含む)
    // vがuの部分木 -> wがuの部分木 && vがwの部分木
    // それ以外 -> wがlca(u, v) or wがuかvのどちらか片方のみ部分木として含む
    bool is_contained_path(int u, int v, int w){
      if(lca(u, v) == w) return true;
      return is_contained_subtree(w, u) ^ is_contained_subtree(w, v);
    }
    // u -> vのパスのk番目, k = 0ならu, k > パス長なら-1
    int kth_vertex_on_path(int u, int v, int k){
      int l = lca(u, v), dlu = depth[u] - depth[l];
      if(dlu > k) return la(u, k);
      k = depth[v] - depth[l] - k + dlu;
      if(k < 0) return -1;
      return la(v, k);
    }
    // 任意のパスはO(log(n))個の区間になる
    // 頂点[0, n)の中でuの位置
    int index_vertex(int u){
      return in[u];
    }
    // 辺[1, n)の中でe{s, t}の位置, 辺が存在しない場合は-1
    int index_edge(int s, int t){
      if(in[s] > in[t]) std::swap(s, t);
      if(parent[t] != s) return -1;
      return in[s] + 1;
    }
    using path = vec<std::pair<int, int>>;
    // 順序を気にせずO(log(n))個の区間を列挙
    path unordered_path(int u, int v, bool is_edge = false){
      path ret;
      for(;; v = parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v);
        if(head[u] == head[v]) break;
        ret.push_back({in[head[v]], in[v] + 1});
      }
      ret.push_back({in[u] + is_edge, in[v] + 1});
      return ret;
    }
    // {u->lcaのパス, lca->vのパス}
    std::pair<path, path> ordered_path(int u, int v, bool is_edge = false){
      bool is_swaped = false;
      std::pair<path, path> ret;
      path &a = ret.first, &b = ret.second;
      for(;; v = parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v), std::swap(a, b), is_swaped ^= 1;
        if(head[u] == head[v]) break;
        b.push_back({in[head[v]], in[v] + 1});
      }
      b.push_back({in[u] + is_edge, in[v] + 1});
      if(is_swaped) std::swap(a, b);
      std::reverse(b.begin(), b.end());
      return {a, b};
    }
    template<typename F>
    void update_path(int u, int v, F &f, bool is_edge = false){
      for(;; v = parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v);
        if(head[u] == head[v]) break;
        f(in[head[v]], in[v] + 1);
      }
      f(in[u] + is_edge, in[v] + 1);
    }
    template<typename F, typename G, typename Val>
    Val query_path(int u, int v, F &f, G &g, Val z, bool is_edge = false){
      Val l = z, r = z;
      for(;; v = parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v), std::swap(l, r);
        if(head[u] == head[v]) break;
        r = g(f(in[head[v]], in[v] + 1), r);
      }
      r = g(f(in[u] + is_edge, in[v] + 1), r);
      return g(l, r);
    }
    template<typename F, typename G, typename H, typename Val>
    Val query_path_flip(int u, int v, F &f, G &g, H &h, Val z, bool is_edge = false){
      bool is_swaped = false;
      Val l = z, r = z;
      for(;; v = parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v), std::swap(l, r), is_swaped ^= 1;
        if(head[u] == head[v]) break;
        r = g(f(in[head[v]], in[v] + 1), r);
      }
      r = g(f(in[u] + is_edge, in[v] + 1), r);
      if(is_swaped) std::swap(l, r);
      return g(h(l), r);
    }
  };

  template<typename tree>
  struct weighted_dist{
  private:
    using edge = typename tree::__edge;
    using weight = typename edge::weight;
    segment_tree<range_sum<weight>> seg;
    tree &t;
  public:
    weighted_dist(tree &_t): t(_t){
      int n = t.n;
      vec<weight> tmp(n, 0);
      for(int i = 0; i < n; i++){
        for(auto &e : t[i]){
          if(t.dep(i) > t.dep(e.t)) continue;
          tmp[t.index_hld(e.t)] = e.wei();
        }
      }
      seg = segment_tree<range_sum<weight>>(tmp);
    }
    // 辺{u, v}(両方向)の重みをw足す　
    void add_weight(int u, int v, weight w){
      auto tmp = seg.get(t.index_hld(v));
      seg.set(t.index_hld(v), tmp + w);
    }
    // u-vパスの重み付きの距離
    weight dist(int u, int v){
      return t.query_path(
        u, v, [&](int  l, int r){return seg.query(l, r);},
        [&](weight a, weight b){return a + b;}, 0, true
      );
    }
    /*
    // u-vパスで重みの総和がw以上になる初めての頂点, ない場合は-1, w <= 0ならu
    int lower_bound(int u, int v, weight w){
      if(w <= 0) return u;
      int l = t.lca(u, v);
      for(;; v = t.hld_p->parent[t.hld_p->head[u]]){
        int a = (t.hld_p->head[u] == t.hld_p->head[l] ? l : t.hld_p->in[t.hld_p->head[u]]);
        int b = t.hld_p->in[u];
        // a -> b
        weight z = seg.query(t.hld_p->in[a], t.hld_p->in[b] + 1);
        if(z >= w){
          
        }
        w -= z;
        if(t.hld_p->head[u] == t.hld_p->head[l]) break;
      }
      return -1;
    }
    */
  };

  // pre_order 頂点を初めて訪れた時刻を記録
  // in_pre: 初めて訪れた時刻
  // out_pre: in_pre以降に初めてvより上のノードが現れる時刻, 区間[in_pre, out_pre)は部分木
  // rev_pre: in_preの順番に頂点を並び替えた状態
  template<typename edge>
  struct dfs_order{
    vec<int> subsize, depth, parent;
    vec<int> in_pre, out_pre, rev_pre;// 訪れた順番(サイズN)
    vec<int> in_path, out_path, rev_path;// 戻る辺も考慮する(サイズ2N-1)
    dfs_order(vec<vec<edge>> &g, int root){
      build(g, root);
    }
    void dfs_build_inner(int cur, int par, int dep, int &tpath, int &tpre, vec<vec<edge>> &g){
      depth[cur] = dep;
      parent[cur] = par;
      in_path[cur] = out_path[cur] = tpath;
      rev_path[tpath++] = cur;
      in_pre[cur] = out_pre[cur] = tpre;
      rev_pre[tpre++] = cur;
      for(int i = 0; i < g[cur].size(); i++){
        int to = g[cur][i].to();
        if(to == par) continue;
        dfs_build_inner(to, cur, dep + 1, tpath, tpre, g);
        subsize[cur] += subsize[to];
        out_path[cur] = tpath;
        rev_path[tpath++] = cur;
      }
      out_pre[cur] = tpre;
    }
    void build(vec<vec<edge>> &g, int root){
      int n = g.size();
      depth.resize(n), parent.resize(n), subsize.resize(n, 1);
      in_pre.resize(n), out_pre.resize(n), rev_pre.resize(n);
      in_path.resize(n), out_path.resize(n), rev_path.resize(2 * n - 1);
      int a = 0, b = 0;
      dfs_build_inner(root, -1, 0, a, b, g);
    }
    // vがuの部分木に含まれるか(u自身も部分木)
    bool is_contained_subtree(int u, int v){
      return in_path[u] <= in_path[v] && out_path[v] <= out_path[u];
    }
    // [in_pre, out_pre)がuの部分木中に存在する頂点番号
    std::pair<int, int> subtree_to_segment(int u){
      return {in_pre[u], out_pre[u]};
    }
    template<typename F>
    void update_subtree(int u, F f, bool is_edge = false){
      f(in_pre[u] + is_edge, out_pre[u]);
    }
    template<typename F>
    auto query_subtree(int u, F f, bool is_edge = false){
      return f(in_pre[u] + is_edge, out_pre[u]);
    }
  };

  template<typename edge>
  struct bfs_order{
    vec<int> parent;
    vec<int> in_bfs, rev_bfs, child_in, parent_index;
    bfs_order(vec<vec<edge>> &g, int root){
      build(root, g);
    }
    void build(int root, vec<vec<edge>> &g){
      int n = g.size();
      in_bfs.resize(n);
      rev_bfs.resize(n);
      child_in.resize(n, -1);
      parent_index.resize(n, -1);
      parent.resize(n);
      std::queue<std::pair<int, int>> q;
      q.push({root, -1});
      int t = 0;
      while(!q.empty()){
        auto [v, p] = q.front();
        q.pop();
        parent[v] = p;
        if(p != -1 && child_in[p] == -1) child_in[p] = t;
        rev_bfs[t] = v;
        in_bfs[v] = t++;
        for(int i = 0; i < g[v].size(); i++){
          if(g[v][i].to() == p){
            parent_index[v] = i;
          }else{
            q.push({g[v][i].to(), v});
          }
        }
      }
    }
    // 辺{u, v}がある場合, g[u]のインデックスを返す, 無い場合は-1
    int find_index(int u, int v){
      if(parent[u] == v) return parent_index[u];
      if(parent[v] != u) return -1;
      int k = in_bfs[v] - child_in[u];
      // k >= g[u]の場合辺がない
      if(parent_index[u] != -1 && parent_index[u] <= k) return k + 1;
      return k;
    }
  };
  template<typename edge>
  void tree_diameter_dfs(int cur, int par, typename edge::weight d, typename edge::weight &dmax, int &vmax, vec<vec<edge>> &g){
    if(d > dmax) dmax = d, vmax = cur;
    for(edge &e : g[cur]){
      if(e.to() == par) continue;
      tree_diameter_dfs(e.to(), cur, d + e.wei(), dmax, vmax, g);
    }
  }
  // {直径, s, t}
  template<typename edge>
  std::tuple<typename edge::weight, int, int> tree_diameter(vec<vec<edge>> &g){
    int s = 0, t = 0;
    typename edge::weight d = edge::z();
    tree_diameter_dfs(s, -1, 0, d, t, g);
    s = t, t = 0, d = edge::z();
    tree_diameter_dfs(s, -1, 0, d, t, g);
    return {d, s, t};
  }

  template<typename LCA, typename DFS, typename DEP>
  std::tuple<vec<int>, vec<vec<int>>, int> lca_tree(vec<int> v, LCA lca, DFS dfs_in, DEP dep){
    if(v.empty()) return {{}, {}, -1};
    std::sort(v.begin(), v.end(), [&](int a, int b){return dfs_in(a) < dfs_in(b);});
    v.erase(std::unique(v.begin(), v.end()), v.end());
    std::stack<int> st;
    vec<std::pair<int, int>> E;
    vec<int> V;
    st.push(v[0]);
    for(int i = 1; i < v.size(); i++){
      if(v[i] == v[i - 1]) continue;
      int l = lca(v[i], v[i - 1]);
      while(true){
        int c = st.top();
        st.pop();
        if(st.empty() || dep(st.top()) <= dep(l)){
          st.push(l);
          st.push(v[i]);
          if(dep(c) > dep(l)){
            E.push_back({l, c});
            V.push_back(c);
            V.push_back(l);
          }
          break;
        }
        int p = st.top();
        if(dep(c) > dep(p)){
          E.push_back({p, c});
          V.push_back(c), V.push_back(p);
        }
      }
    }
    while(st.size() >= 2){
      int c = st.top();
      st.pop();
      int p = st.top();
      if(c != p) E.push_back({p, c}), V.push_back(c), V.push_back(p);
    }
    if(!st.empty()) V.push_back(st.top());
    std::sort(V.begin(), V.end());
    V.erase(std::unique(V.begin(), V.end()), V.end());
    int root = 0;
    int n = V.size();
    for(int i = 1; i < n; i++){
      if(dep(V[root]) > dep(V[i])){
        root = i;
      }
    }
    vec<vec<int>> G(n);
    for(auto [u, v] : E){
      int u2 = std::lower_bound(V.begin(), V.end(), u) - V.begin();
      int v2 = std::lower_bound(V.begin(), V.end(), v) - V.begin();
      G[u2].push_back(v2);
      G[v2].push_back(u2);
    }
    return {V, G, root};
  }
  template<typename edge>
  std::tuple<vec<vec<int>>, int, vec<int>, vec<int>, vec<int>> centroid_decomposition(vec<vec<edge>> &g){
    int n = g.size();
    assert(n);
    vec<vec<int>> G(n);
    std::vector<int> size_i(n, 0), dep_i(n, std::numeric_limits<int>::max()), par_i(n, -1);

    auto add_edge = [&](int p, int c)->void{
      G[p].push_back(c);
      G[c].push_back(p);
      par_i[c] = p;
    };
    auto find_centroid = [&](auto &&find_centroid, int v, int p, const int N, const int8_t rank)->std::pair<int, int>{
      int sz = 1;
      for(edge &e: g[v]){
        if(e.t == p || dep_i[e.t] < rank) continue;
        auto [sz_c, cent_c] = find_centroid(find_centroid, e.t, v, N, rank);
        if(sz_c == -1) return {-1, cent_c};
        size_i[e.t] = sz_c, sz += sz_c;
      }
      //サイズが半分以上になったとき
      if(sz * 2 >= N){
        size_i[v] = N;
        dep_i[v] = rank;
        for(edge &e: g[v]){
          if(e.t == p || dep_i[e.t] < rank) continue;
          auto [sz_c, cent_c] = find_centroid(find_centroid, e.t, -1, size_i[e.t], rank + 1);
          assert(sz_c == -1);
          add_edge(v, cent_c);
        }
        if(p != -1){
          auto [sz_c, cent_c] = find_centroid(find_centroid, p, -1, N - sz, rank + 1);
          assert(sz_c == -1);
          add_edge(v, cent_c);
        }
        return {-1, v};// 重心を発見
      }
      return {sz, -1};
    };
    int root = find_centroid(find_centroid, 0, -1, n, 0).second;
    return {G, root, size_i, dep_i, par_i};
  }
};
#endif