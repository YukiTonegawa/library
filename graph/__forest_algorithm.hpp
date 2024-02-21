#ifndef _FOREST_ALGORITHM_H_
#define _FOREST_ALGORITHM_H_
#include <vector>
#include <algorithm>
#include <queue>
#include <numeric>
#include "edge.hpp"

namespace forest_algorithm{
  template<typename T>
  using vec = std::vector<T>;
  // 部分木のサイズ, 深さ, 親, 根
  template<typename forest>
  void simple_dfs(forest &f, int r){
    using edge = typename forest::edge;
    int n = f.n;
    auto simple_dfs_f = [&](auto simple_dfs_f, int cur, int par, int dep)->void{
      f.parent[cur] = par;
      f.depth[cur] = dep;
      f.root[cur] = r;
      for(edge &e : f[cur]){
        if(e.to() == par) continue;
        simple_dfs_f(simple_dfs_f, e.to(), cur, dep + 1);
        f.size[cur] += f.size[e.to()];
      }
    };
    if(f.depth.empty()){
      f.depth.resize(n, -1);
      f.parent.resize(n, -1);
      f.size.resize(n, -1);
      f.root.resize(n, -1);
    }
    if(f.depth[r] == -1) simple_dfs_f(simple_dfs_f, r, -1, 0);
  }

  template<typename forest>
  struct hld{
    using edge = typename forest::_edge;
    forest &f;
    vec<int> in, out, head, rev, heavy;
    hld(){}
    hld(forest &_f): f(_f){
      int n = f.n;
      f.subsize.resize(n, -1);
      f.depth.resize(n, -1);
      f.parent.resize(n, -1);
      f.root.resize(n, -1);
      in.resize(n, -1);
      out.resize(n);
      head.resize(n);
      rev.resize(n);
      heavy.resize(n);
    }
    void dfs_sz(int cur, int par, int dep, int r){
      f.root[cur] = r;
      f.depth[cur] = dep;
      f.parent[cur] = par;
      f.subsize[cur] = 1;
      for(int i = 0; i < f[cur].size(); i++){
        edge &e = f[cur][i];
        if(e.t == par) continue;
        dfs_sz(e.t, cur, dep + 1, r);
        f.subsize[cur] += f.subsize[e.t];
        if(heavy[cur] == -1 || f.subsize[f[cur][heavy[cur]].t] < f.subsize[e.t]){
          heavy[cur] = i;
        }
      }
    }
    void dfs_hld(int cur, int par, int &times){
      in[cur] = times++;
      rev[in[cur]] = cur;
      int h = heavy[cur];
      if(h != -1){
        head[f[cur][h].t] = head[cur];
        dfs_hld(f[cur][h].t, cur, times);
      }
      for(int i = 0; i < h; i++){
        if(f[cur][i].t == par) continue;
        head[f[cur][i].t] = f[cur][i].t;
        dfs_hld(f[cur][i].t, cur, times);
      }
      for(int i = h + 1; i < f[cur].size(); i++){
        if(f[cur][i].t == par) continue;
        head[f[cur][i].t] = f[cur][i].t;
        dfs_hld(f[cur][i].t, cur, times);
      }
      out[cur] = times;
    }
    // uを含む連結成分がすでに初期化されているか
    bool initialized(int u){
      return in.size() == f.n && in[u] != -1;
    }
    // rを含む連結成分をrを根として初期化　, すでに連結成分が初期化されている場合は何もしない
    void build(int r){
      if(initialized(r)) return;
      dfs_sz(r, -1, 0, r);
      static int t = 0;
      dfs_hld(r, -1, t);
    }
    bool is_root(int u){
      assert(initialized(u));
      return u == f.root[u];
    }
    int find_root(int u){
      assert(initialized(u));
      return f.root[u];
    }
    bool same_tree(int u, int v){
      assert(initialized(u) && initialized(v));
      return f.root[u] == f.root[v];
    }
    int la(int v, int k){
      assert(initialized(v));
      if(f.depth[v] < k) return -1;
      while(true){
        int u = head[v];
        if(in[v] - k >= in[u]) return rev[in[v] - k];
        k -= in[v] - in[u] + 1;
        v = f.parent[u];
      }
    }
    // 同じ木に無い場合は-1
    int lca(int u, int v){
      assert(initialized(u) && initialized(v));
      for(;v != -1; v = f.parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v);
        if(head[u] == head[v]) return u;
      }
      return -1;
    }
    bool is_contained_subtree(int u, int v){
      assert(initialized(u) && initialized(v));
      if(f.depth[u] > f.depth[v]) return false;
      return u == f.la(v, f.depth[v] - f.depth[u]);
    }
    int kth_vertex_on_path(int u, int v, int k){
      int l = lca(u, v);
      if(l == -1) return -1;
      int dlu = f.depth[u] - f.depth[l];
      if(dlu > k) return la(u, k);
      k = f.depth[v] - f.depth[l] - k + dlu;
      if(k < 0) return -1;
      return la(v, k);
    }
    // 重みなしの距離
    int dist(int u, int v){
      int l = lca(u, v);
      return f.depth[u] + f.depth[v] - 2 * f.depth[l];
    }
    // 任意のパスはO(log(n))個の区間になる

    // 任意のパスはO(log(n))個の区間になる
    // 頂点[0, n)の中でuの位置
    int index_vertex(int u){
      assert(initialized(u));
      return in[u];
    }
    // 辺[1, n)の中でe{s, t}の位置, 辺が存在しない場合は-1
    int index_edge(int s, int t){
      assert(initialized(s) && initialized(t));
      if(in[s] > in[t]) std::swap(s, t);
      if(f.parent[t] != s) return -1;
      return in[s] + 1;
    }
    using path = vec<std::pair<int, int>>;
    // 順序を気にせずO(log(n))個の区間を列挙
    // 連結成分が違う場合は空のvector
    path unordered_path(int u, int v, bool is_edge = false){
      if(!same_tree(u, v)) return {};
      path ret;
      for(;; v = f.parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v);
        if(head[u] == head[v]) break;
        ret.push_back({in[head[v]], in[v] + 1});
      }
      ret.push_back({in[u] + is_edge, in[v] + 1});
      return ret;
    }
    // {lca->uのパス, lca->vのパス}
    // 連結成分が違う場合は空のvector
    std::pair<path, path> ordered_path(int u, int v, bool is_edge = false){
      if(!same_tree(u, v)) return {};
      bool is_swaped = false;
      std::pair<path, path> ret;
      path &a = ret.first, &b = ret.second;
      for(;; v = f.parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v), std::swap(a, b), is_swaped ^= 1;
        if(head[u] == head[v]) break;
        b.push_back({in[head[v]], in[v] + 1});
      }
      b.push_back({in[u] + is_edge, in[v] + 1});
      if(is_swaped) std::swap(a, b);
      std::reverse(a.begin(), a.end());
      std::reverse(b.begin(), b.end());
      return {a, b};
    }
    // 連結成分が違う場合は何もしない
    template<typename F>
    void update_path(int u, int v, F &f, bool is_edge = false){
      if(!same_tree(u, v)) return;
      for(;; v = f.parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v);
        if(head[u] == head[v]) break;
        f(in[head[v]], in[v] + 1);
      }
      f(in[u] + is_edge, in[v] + 1);
    }
    // 連結成分が違う場合は単位元
    template<typename F, typename G, typename Val>
    Val query_path(int u, int v, F &f, G &g, Val z, bool is_edge = false){
      if(!same_tree(u, v)) return z;
      Val l = z, r = z;
      for(;; v = f.parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v), std::swap(l, r);
        if(head[u] == head[v]) break;
        r = g(f(in[head[v]], in[v] + 1), r);
      }
      r = g(f(in[u] + is_edge, in[v] + 1), r);
      return g(l, r);
    }
    // 連結成分が違う場合は単位元
    template<typename F, typename G, typename H, typename Val>
    Val query_path_flip(int u, int v, F &f, G &g, H &h, Val z, bool is_edge = false){
      if(!same_tree(u, v)) return z;
      bool is_swaped = false;
      Val l = z, r = z;
      for(;; v = f.parent[head[v]]){
        if(in[u] > in[v]) std::swap(u, v), std::swap(l, r), is_swaped ^= 1;
        if(head[u] == head[v]) break;
        r = g(f(in[head[v]], in[v] + 1), r);
      }
      r = g(f(in[u] + is_edge, in[v] + 1), r);
      if(is_swaped) std::swap(l, r);
      return g(h(l), r);
    }
  };

  // pre_order 頂点を初めて訪れた時刻を記録
  // in_pre: 初めて訪れた時刻
  // out_pre: in_pre以降に初めてvより上のノードが現れる時刻, 区間[in_pre, out_pre)は部分木
  // rev_pre: in_preの順番に頂点を並び替えた状態
  template<typename forest>
  struct dfs_order{
    vec<int> in_pre, out_pre, rev_pre;// 訪れた順番(サイズN)
    vec<int> in_path, out_path, rev_path;// 戻る辺も考慮する(サイズ2N-1)
    forest &f;
    dfs_order(forest &_f): f(_f){
      int n = f.n;
      f.subsize.resize(n, -1);
      f.depth.resize(n, -1);
      f.parent.resize(n, -1);
      f.root.resize(n, -1);
      in_pre.resize(n, -1);
      out_pre.resize(n);
      rev_pre.resize(n);
      in_path.resize(n);
      out_path.resize(n);
      rev_path.resize(2 * n - 1);
    }
    void dfs_build_inner(int cur, int par, int r, int dep, int &tpath, int &tpre){
      f.depth[cur] = dep;
      f.parent[cur] = par;
      f.root[cur] = r;
      in_path[cur] = out_path[cur] = tpath;
      rev_path[tpath++] = cur;
      in_pre[cur] = out_pre[cur] = tpre;
      rev_pre[tpre++] = cur;
      for(int i = 0; i < f[cur].size(); i++){
        int to = f[cur][i].to();
        if(to == par) continue;
        dfs_build_inner(to, cur, r, dep + 1, tpath, tpre);
        f.subsize[cur] += f.subsize[to];
        out_path[cur] = tpath;
        rev_path[tpath++] = cur;
      }
      out_pre[cur] = tpre;
    }
    bool initialized(int r){
      return in_pre.size() == f.n && in_pre[r] != -1;
    }
    void build(int r){
      if(initialized(r)) return;
      static int a = 0, b = 0;
      dfs_build_inner(r, -1, r, 0, a, b);
    }
    // vがuの部分木に含まれるか(u自身も部分木)
    bool is_contained_subtree(int u, int v){
      assert(initialized(u) && initialized(v));
      return in_path[u] <= in_path[v] && out_path[v] <= out_path[u];
    }
    // u->vパス(最短経路)にwが含まれるか(端点も含む)
    // vがuの部分木 -> wがuの部分木 && vがwの部分木
    // それ以外 -> wがuかvのどちらかを部分木として含む
    bool is_contained_path(int u, int v, int w){
      assert(initialized(u) && initialized(v) && initialized(w));
      if(in_path[u] > in_path[v]) std::swap(u, v);
      if(is_contained_subtree(u, v)) return in_path[u] <= in_path[w] && is_contained_subtree(w, v);
      return is_contained_subtree(w, u) || is_contained_subtree(w, v);
    }
    // [in_pre, out_pre)がuの部分木中に存在する頂点番号
    std::pair<int, int> subtree_to_segment(int u){
      assert(initialized(u));
      return {in_pre[u], out_pre[u]};
    }
    template<typename F>
    void update_subtree(int u, F f, bool is_edge = false){
      assert(initialized(u));
      f(in_pre[u] + is_edge, out_pre[u]);
    }
    template<typename F>
    auto query_subtree(int u, F f, bool is_edge = false){
      assert(initialized(u));
      return f(in_pre[u] + is_edge, out_pre[u]);
    }
  };

  template<typename forest>
  struct bfs_order{
    vec<int> in_bfs, rev_bfs, child_in, parent_index;
    forest &f;
    bfs_order(forest &_f): f(_f){
      int n = f.n;
      in_bfs.resize(n, -1);
      rev_bfs.resize(n);
      child_in.resize(n, -1);
      parent_index.resize(n, -1);
      f.parent.resize(n, -1);
      f.root.resize(n, -1);
    }
    void build(int r){
      if(initialized(r)) return;
      std::queue<std::pair<int, int>> q;
      q.push({r, -1});
      static int t = 0;
      while(!q.empty()){
        auto [v, p] = q.front();
        q.pop();
        f.parent[v] = p;
        if(p != -1 && child_in[p] == -1) child_in[p] = t;
        rev_bfs[t] = v;
        in_bfs[v] = t++;
        for(int i = 0; i < f[v].size(); i++){
          if(f[v][i].to() == p){
            parent_index[v] = i;
          }else{
            q.push({f[v][i].to(), v});
          }
        }
      }
    }
    bool initialized(int r){
      return in_bfs.size() == f.n && in_bfs[r] != -1;
    }
    // 辺{u, v}がある場合, g[u]のインデックスを返す, 無い場合は-1
    int find_index(int u, int v){
      assert(initialized(u) && initialized(v));
      if(f.parent[u] == v) return parent_index[u];
      if(f.parent[v] != u) return -1;
      int k = in_bfs[v] - child_in[u];
      // k >= g[u]の場合辺がない
      if(parent_index[u] != -1 && parent_index[u] <= k) return k + 1;
      return k;
    }
  };
  template<typename forest>
  void diameter_dfs(int cur, int par, typename forest::_edge::weight d, typename forest::_edge::weight &dmax, int &vmax, forest &f){
    if(d > dmax) dmax = d, vmax = cur;
    for(auto &e : f[cur]){
      if(e.to() == par) continue;
      diameter_dfs(e.to(), cur, d + e.wei(), dmax, vmax, f);
    }
  }
  // {直径, s, t}
  template<typename forest>
  std::tuple<typename forest::_edge::weight, int, int> diameter(forest &f, int r){
    int s = r, t = 0;
    typename forest::_edge::weight d = forest::_edge::z();
    diameter_dfs(s, -1, 0, d, t, f);
    s = t, t = 0, d = forest::_edge::z();
    diameter_dfs(s, -1, 0, d, t, f);
    return {d, s, t};
  }
};
#endif