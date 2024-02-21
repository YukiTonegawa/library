#ifndef _TREE_CONTOUR_H_
#define _TREE_CONTOUR_H_
#include "tree.hpp"
#include "fast_lca.hpp"
#include "../data_structure/range_query/accumulate1d.hpp"
#include "../data_structure/segment_tree/binary_indexed_tree_range_add.hpp"

namespace tree_algorithm{
  // 各重心が管理する頂点の番号が連続になるように頂点番号を振り直す
  struct centroid_decomposition_structure{
    int n, r;
    std::vector<int> cnv, rev; // cnv: 元の頂点番号 -> 振り直された番号, rev: 振り直された番号 -> 元の番号
    std::vector<int> begin; // st[g] := gが管理する頂点番号の最小値
    std::vector<int> dep, par; // dep[g] := 重心の深さ, par[g] := 1つ上の重心
    simple_tree traw;
    
    // 重心gからbfs(変更後の頂点番号)
    // {頂点, 親}
    std::vector<std::pair<int, int>> bfs_raw(int g){
      int j = 0;
      std::vector<std::pair<int, int>> ret(size_raw(g));
      std::queue<std::pair<int, int>> q;
      q.push({g, -1});
      while(!q.empty()){
        ret[j++] = q.front();
        auto [v, par] = q.front();
        q.pop();
        for(auto &e : traw[v]){
          if(e.t == par || dep[e.t] <= dep[g]) continue;
          q.push({e.t, v});
        }
      }
      return ret;
    }
    // bfs_rawと比べて出力が子部分木ごとに分かれている
    std::vector<std::vector<std::pair<int, int>>> bfs_raw2(int g){
      std::vector<std::vector<std::pair<int, int>>> ret;
      std::queue<std::pair<int, int>> q;
      for(auto &e : traw[g]){
        if(dep[e.t] <= dep[g]) continue;
        q.push({e.t, g});
        ret.push_back({});
        while(!q.empty()){
          ret.back().push_back(q.front());
          auto [v, par] = q.front();
          q.pop();
          for(auto &e : traw[v]){
            if(e.t == par || dep[e.t] <= dep[g]) continue;
            q.push({e.t, v});
          }
        }
      }
      return ret;
    }
    centroid_decomposition_structure(): traw(0){}

    template<typename edge>
    centroid_decomposition_structure(tree<edge> &_t): n(_t.n), cnv(n), rev(n), begin(n), dep(n), par(n), traw(n){
      simple_tree c = _t.centroid_decomposition();
      r = c.r;
      int j = 0;
      auto f = [&](auto &&f, int cur, int p, int d)->void{
        int mn = j;
        for(auto e : c[cur]){
          if(e.t == p) continue;
          f(f, e.t, cur, d + 1);
        }
        cnv[cur] = j;
        begin[j] = mn;
        par[j] = p;
        dep[j] = d;
        rev[j++] = cur;
      };
      f(f, c.r, -1, 0);
      for(int i = 0; i < n; i++) if(par[i] != -1) par[i] = cnv[par[i]];
      traw.r = encode(_t.r);
      for(int i = 0; i < n; i++){
        int ci = cnv[i];
        for(auto &e : _t[i]) traw.add_edge(ci, {ci, cnv[e.t]});
      }
    }
    // gの管理する頂点の番号の最小値　(変更後の頂点番号)
    int begin_raw(int g){
      return begin[g];
    }
    // 重心gの管理する部分グラフのサイズ(変更後の頂点番号)
    int size_raw(int g){
      return g - begin[g] + 1;
    }
    // 重心gの管理する部分グラフのサイズ(元の頂点番号)
    int size(int g){
      int g2 = cnv[g];
      return g2 - begin[g2] + 1;
    }
    // 変更後の頂点番号でgが　cを含むか
    bool contain_raw(int g, int c){
      return begin[g] <= c && c <= g;
    }
    // 元の頂点番号でgが　cを含むか
    bool contain(int g, int c){
      return contain_raw(cnv[g], cnv[c]);
    }
    // 元の頂点番号 -> 変更後の番号
    int encode(int v){
      return cnv[v];
    }
    // 変更後の頂点番号 -> 元の番号
    int decode(int v){
      return rev[v];
    }
  };
  // 各重心が管理する頂点の番号が連続になるように頂点番号を振り直す
  template<typename edge>
  struct centroid_decomposition_structure_weighted{
    using weight = typename edge::weight;
    int n, r;
    std::vector<int> cnv, rev; // cnv: 元の頂点番号 -> 振り直された番号, rev: 振り直された番号 -> 元の番号
    std::vector<int> begin; // st[g] := gが管理する頂点番号の最小値
    std::vector<int> dep, par; // dep[g] := 重心の深さ, par[g] := 1つ上の重心
    tree<edge> traw;
    
    // 重心gからbfs(変更後の頂点番号)
    // {頂点, 親}
    std::vector<std::pair<int, int>> bfs_raw(int g){
      int j = 0;
      std::vector<std::pair<int, int>> ret(size_raw(g));
      std::queue<std::pair<int, int>> q;
      q.push({g, -1});
      while(!q.empty()){
        ret[j++] = q.front();
        auto [v, par] = q.front();
        q.pop();
        for(auto &e : traw[v]){
          if(e.t == par || dep[e.t] <= dep[g]) continue;
          q.push({e.t, v});
        }
      }
      return ret;
    }
    // bfs_rawと比べて出力が子部分木ごとに分かれている
    std::vector<std::vector<std::pair<int, int>>> bfs_raw2(int g){
      std::vector<std::vector<std::pair<int, int>>> ret;
      std::queue<std::pair<int, int>> q;
      for(auto &e : traw[g]){
        if(dep[e.t] <= dep[g]) continue;
        q.push({e.t, g});
        ret.push_back({});
        while(!q.empty()){
          ret.back().push_back(q.front());
          auto [v, par] = q.front();
          q.pop();
          for(auto &e : traw[v]){
            if(e.t == par || dep[e.t] <= dep[g]) continue;
            q.push({e.t, v});
          }
        }
      }
      return ret;
    }
    // 重心から距離の短い順に探索 {重心からの距離, 頂点, 親}
    std::vector<std::tuple<weight, int, int>> shortest_raw(int g){
      int j = 0;
      using t = std::tuple<weight, int, int>; // {距離, 頂点, 親}
      std::vector<t> ret(size_raw(g));
      std::priority_queue<t, std::vector<t>, std::greater<t>> pq;
      pq.push({0, g, -1});
      while(!pq.empty()){
        ret[j++] = pq.top();
        auto [d, v, par] = pq.top();
        pq.pop();
        for(auto &e : traw[v]){
          if(e.t == par || dep[e.t] <= dep[g]) continue;
          pq.push({d + e.w, e.t, v});
        }
      }
      return ret;
    }
    // 重心から距離の短い順に探索 {重心からの距離, 頂点, 親}
    std::vector<std::vector<std::tuple<weight, int, int>>> shortest_raw2(int g){
      using t = std::tuple<weight, int, int>; // {距離, 頂点, 親}
      std::vector<std::vector<t>> ret;
      std::priority_queue<t, std::vector<t>, std::greater<t>> pq;
      for(auto &e : traw[g]){
        if(dep[e.t] <= dep[g]) continue;
        pq.push({e.w, e.t, g});
        ret.push_back({});
        while(!pq.empty()){
          ret.back().push_back(pq.top());
          auto [d, v, par] = pq.top();
          pq.pop();
          for(auto &e : traw[v]){
            if(e.t == par || dep[e.t] <= dep[g]) continue;
            pq.push({d + e.w, e.t, v});
          }
        }
      }
      return ret;
    }
    centroid_decomposition_structure_weighted(): traw(0){}

    centroid_decomposition_structure_weighted(tree<edge> &_t): n(_t.n), cnv(n), rev(n), begin(n), dep(n), par(n), traw(n){
      simple_tree c = _t.centroid_decomposition();
      r = c.r;
      int j = 0;
      auto f = [&](auto &&f, int cur, int p, int d)->void{
        int mn = j;
        for(auto e : c[cur]){
          if(e.t == p) continue;
          f(f, e.t, cur, d + 1);
        }
        cnv[cur] = j;
        begin[j] = mn;
        par[j] = p;
        dep[j] = d;
        rev[j++] = cur;
      };
      f(f, c.r, -1, 0);
      for(int i = 0; i < n; i++) if(par[i] != -1) par[i] = cnv[par[i]];
      traw.r = encode(_t.r);
      for(int i = 0; i < n; i++){
        int ci = cnv[i];
        for(auto e : _t[i]){
          e.s = ci, e.t = cnv[e.t];
          traw.add_edge(ci, e);
        }
      }
    }
    // gの管理する頂点の番号の最小値　(変更後の頂点番号)
    int begin_raw(int g){
      return begin[g];
    }
    // 重心gの管理する部分グラフのサイズ(変更後の頂点番号)
    int size_raw(int g){
      return g - begin[g] + 1;
    }
    // 重心gの管理する部分グラフのサイズ(元の頂点番号)
    int size(int g){
      int g2 = cnv[g];
      return g2 - begin[g2] + 1;
    }
    // 変更後の頂点番号でgが　cを含むか
    bool contain_raw(int g, int c){
      return begin[g] <= c && c <= g;
    }
    // 元の頂点番号でgが　cを含むか
    bool contain(int g, int c){
      return contain_raw(cnv[g], cnv[c]);
    }
    // 元の頂点番号 -> 変更後の番号
    int encode(int v){
      return cnv[v];
    }
    // 変更後の頂点番号 -> 元の番号
    int decode(int v){
      return rev[v];
    }
  };
  
  template<typename Val>
  struct contour_sum{
    int n;
    contour_sum(){}
    struct node{
      vec<accumulate1d<Val>> sum;
      vec<vec<int>> dep_start;
    };
    std::vector<std::vector<std::array<int, 4>>> in; // [頂点番号][深さ][全体でのin / 部分木でのin / 部分木番号 / 重心との距離]
    std::vector<node> nodes;
    fast_lca fl;
    centroid_decomposition_structure cd;

    template<typename edge>
    void build(tree<edge> t, const std::vector<Val> &val){
      n = t.n;
      cd = centroid_decomposition_structure(t);
      fl = fast_lca(cd.traw.g, cd.r);
      nodes.resize(n);
      in.resize(n);
      std::vector<Val> val2(n);
      for(int i = 0; i < n; i++){
        in[i].resize(cd.dep[i] + 1);
        val2[cd.encode(i)] = val[i];
      }
      std::vector<int> dtmp(n);
      for(int i = 0; i < n; i++){
        auto v1 = cd.bfs_raw(i);
        auto v2 = cd.bfs_raw2(i);
        nodes[i].sum.resize(v2.size() + 1);
        nodes[i].dep_start.resize(v2.size() + 1);
        int dep = cd.dep[i];
        for(int j = 0; j < v2.size(); j++){
          std::vector<Val> tmp(v2[j].size());
          for(int k = 0; k < v2[j].size(); k++){
            int u = v2[j][k].first;
            tmp[k] = val2[u];
            dtmp[u] = (k == 0 ? 0 : dtmp[v2[j][k].second] + 1);
            if(dtmp[u] == nodes[i].dep_start[j].size()) nodes[i].dep_start[j].push_back(k);
            in[u][dep][1] = k;
            in[u][dep][2] = j;
            in[u][dep][3] = fl.dist(u, i);
          }
          nodes[i].sum[j] = accumulate1d<Val>(tmp);
        }
        std::vector<Val> tmp(v1.size());
        for(int k = 0; k < v1.size(); k++){
          int u = v1[k].first;
          tmp[k] = val2[u];
          dtmp[u] = (k == 0 ? 0 : dtmp[v1[k].second] + 1);
          if(dtmp[u] == nodes[i].dep_start.back().size()) nodes[i].dep_start.back().push_back(k);
          in[u][dep][0] = k;
        }
        nodes[i].sum.back() = accumulate1d<Val>(tmp);
      }
    }
    // vからの距離が[ldist, rdist)の頂点の数
    int count(int v, int ldist, int rdist){
      int ret = 0;
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        int dist = (u == v ? 0 : in[v][d][3]);
        int M = nodes[u].dep_start.back().size(), L = std::max(0, ldist - dist), R = std::min(rdist - dist, M);
        if(L < R) ret += (R == M ? nodes[u].sum.back().sum.size() : nodes[u].dep_start.back()[R]) - nodes[u].dep_start.back()[L];
        if(u != v){
          int kc = in[v][d][2];
          M = nodes[u].dep_start[kc].size(), L = std::max(0, L - 1), R = std::min(rdist - dist - 1, M);
          if(L < R) ret -= (R == M ? nodes[u].sum[kc].sum.size(): nodes[u].dep_start[kc][R]) - nodes[u].dep_start[kc][L];
        }
      }
      return ret;
    }
    // vからの距離が[ldist, rdist)の頂点の和
    Val query(int v, int ldist, int rdist){
      Val ret = 0;
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        int dist = (u == v ? 0 : in[v][d][3]);
        int M = nodes[u].dep_start.back().size(), L = std::max(0, ldist - dist), R = std::min(rdist - dist, M);
        if(L < R) ret += nodes[u].sum.back().query(nodes[u].dep_start.back()[L], R == M ? nodes[u].sum.back().sum.size() : nodes[u].dep_start.back()[R]);
        if(u != v){
          int kc = in[v][d][2];
          M = nodes[u].dep_start[kc].size(), L = std::max(0, L - 1), R = std::min(rdist - dist - 1, M);
          if(L < R) ret -= nodes[u].sum[kc].query(nodes[u].dep_start[kc][L], R == M ? nodes[u].sum[kc].sum.size() : nodes[u].dep_start[kc][R]);
        }
      }
      return ret;
    }
  };
  template<typename Val>
  struct point_add_contour_sum{
    int n;
    point_add_contour_sum(){}
    struct node{
      vec<binary_indexed_tree<Val>> sum;
      vec<vec<int>> dep_start;
    };
    std::vector<std::vector<std::array<int, 4>>> in; // [頂点番号][深さ][全体でのin / 部分木でのin / 部分木番号 / 重心との距離]
    std::vector<node> nodes;
    fast_lca fl;
    centroid_decomposition_structure cd;

    template<typename edge>
    void build(tree<edge> t, const std::vector<Val> &val){
      n = t.n;
      cd = centroid_decomposition_structure(t);
      fl = fast_lca(cd.traw.g, cd.r);
      nodes.resize(n);
      in.resize(n);
      std::vector<Val> val2(n);
      for(int i = 0; i < n; i++){
        in[i].resize(cd.dep[i] + 1);
        val2[cd.encode(i)] = val[i];
      }
      std::vector<int> dtmp(n);
      for(int i = 0; i < n; i++){
        auto v1 = cd.bfs_raw(i);
        auto v2 = cd.bfs_raw2(i);
        nodes[i].sum.resize(v2.size() + 1);
        nodes[i].dep_start.resize(v2.size() + 1);
        int dep = cd.dep[i];
        for(int j = 0; j < v2.size(); j++){
          std::vector<Val> tmp(v2[j].size());
          for(int k = 0; k < v2[j].size(); k++){
            int u = v2[j][k].first;
            tmp[k] = val2[u];
            dtmp[u] = (k == 0 ? 0 : dtmp[v2[j][k].second] + 1);
            if(dtmp[u] == nodes[i].dep_start[j].size()) nodes[i].dep_start[j].push_back(k);
            in[u][dep][1] = k;
            in[u][dep][2] = j;
            in[u][dep][3] = fl.dist(u, i);
          }
          nodes[i].sum[j] = binary_indexed_tree<Val>(tmp);
        }
        std::vector<Val> tmp(v1.size());
        for(int k = 0; k < v1.size(); k++){
          int u = v1[k].first;
          tmp[k] = val2[u];
          dtmp[u] = (k == 0 ? 0 : dtmp[v1[k].second] + 1);
          if(dtmp[u] == nodes[i].dep_start.back().size()) nodes[i].dep_start.back().push_back(k);
          in[u][dep][0] = k;
        }
        nodes[i].sum.back() = binary_indexed_tree<Val>(tmp);
      }
    }
    // vからの距離が[ldist, rdist)の頂点の数
    int count(int v, int ldist, int rdist){
      int ret = 0;
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        int dist = (u == v ? 0 : in[v][d][3]);
        int M = nodes[u].dep_start.back().size(), L = std::max(0, ldist - dist), R = std::min(rdist - dist, M);
        if(L < R) ret += (R == M ? nodes[u].sum.back().M : nodes[u].dep_start.back()[R]) - nodes[u].dep_start.back()[L];
        if(u != v){
          int kc = in[v][d][2];
          M = nodes[u].dep_start[kc].size(), L = std::max(0, L - 1), R = std::min(rdist - dist - 1, M);
          if(L < R) ret -= (R == M ? nodes[u].sum[kc].M : nodes[u].dep_start[kc][R]) - nodes[u].dep_start[kc][L];
        }
      }
      return ret;
    }
    // v += x
    void update(int v, Val x){
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        nodes[u].sum.back().update(in[v][d][0], x);
        if(u != v) nodes[u].sum[in[v][d][2]].update(in[v][d][1], x);
      }
    }
    // vからの距離が[ldist, rdist)の頂点の和
    Val query(int v, int ldist, int rdist){
      Val ret = 0;
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        int dist = (u == v ? 0 : in[v][d][3]);
        int M = nodes[u].dep_start.back().size(), L = std::max(0, ldist - dist), R = std::min(rdist - dist, M);
        if(L < R) ret += nodes[u].sum.back().query(nodes[u].dep_start.back()[L], R == M ? nodes[u].sum.back().M : nodes[u].dep_start.back()[R]);
        if(u != v){
          int kc = in[v][d][2];
          M = nodes[u].dep_start[kc].size(), L = std::max(0, L - 1), R = std::min(rdist - dist - 1, M);
          if(L < R) ret -= nodes[u].sum[kc].query(nodes[u].dep_start[kc][L], R == M ? nodes[u].sum[kc].M : nodes[u].dep_start[kc][R]);
        }
      }
      return ret;
    }
  };
  template<typename Val>
  struct contour_add_point_get{
    int n;
    contour_add_point_get(){}
    struct node{
      vec<binary_indexed_tree_range_add<Val>> sum;
      vec<vec<int>> dep_start;
    };
    std::vector<std::vector<std::array<int, 4>>> in; // [頂点番号][深さ][全体でのin / 部分木でのin / 部分木番号 / 重心との距離]
    std::vector<node> nodes;
    fast_lca fl;
    centroid_decomposition_structure cd;

    template<typename edge>
    void build(tree<edge> t, const std::vector<Val> &val){
      n = t.n;
      cd = centroid_decomposition_structure(t);
      fl = fast_lca(cd.traw.g, cd.r);
      nodes.resize(n);
      in.resize(n);
      std::vector<Val> val2(n);
      for(int i = 0; i < n; i++){
        in[i].resize(cd.dep[i] + 1);
        val2[cd.encode(i)] = val[i];
      }
      std::vector<int> dtmp(n);
      for(int i = 0; i < n; i++){
        auto v1 = cd.bfs_raw(i);
        auto v2 = cd.bfs_raw2(i);
        nodes[i].sum.resize(v2.size() + 1);
        nodes[i].dep_start.resize(v2.size() + 1);
        int dep = cd.dep[i];
        for(int j = 0; j < v2.size(); j++){
          std::vector<Val> tmp(v2[j].size());
          for(int k = 0; k < v2[j].size(); k++){
            int u = v2[j][k].first;
            tmp[k] = val2[u];
            dtmp[u] = (k == 0 ? 0 : dtmp[v2[j][k].second] + 1);
            if(dtmp[u] == nodes[i].dep_start[j].size()) nodes[i].dep_start[j].push_back(k);
            in[u][dep][1] = k;
            in[u][dep][2] = j;
            in[u][dep][3] = fl.dist(u, i);
          }
          nodes[i].sum[j] = binary_indexed_tree_range_add_range_sum<Val>(tmp);
        }
        std::vector<Val> tmp(v1.size());
        for(int k = 0; k < v1.size(); k++){
          int u = v1[k].first;
          tmp[k] = val2[u];
          dtmp[u] = (k == 0 ? 0 : dtmp[v1[k].second] + 1);
          if(dtmp[u] == nodes[i].dep_start.back().size()) nodes[i].dep_start.back().push_back(k);
          in[u][dep][0] = k;
        }
        nodes[i].sum.back() = binary_indexed_tree_range_add_range_sum<Val>(tmp);
      }
    }
    // vからの距離が[ldist, rdist)の頂点の数
    int count(int v, int ldist, int rdist){
      int ret = 0;
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        int dist = (u == v ? 0 : in[v][d][3]);
        int M = nodes[u].dep_start.back().size(), L = std::max(0, ldist - dist), R = std::min(rdist - dist, M);
        if(L < R) ret += (R == M ? nodes[u].sum.back().M : nodes[u].dep_start.back()[R]) - nodes[u].dep_start.back()[L];
        if(u != v){
          int kc = in[v][d][2];
          M = nodes[u].dep_start[kc].size(), L = std::max(0, L - 1), R = std::min(rdist - dist - 1, M);
          if(L < R) ret -= (R == M ? nodes[u].sum[kc].M : nodes[u].dep_start[kc][R]) - nodes[u].dep_start[kc][L];
        }
      }
      return ret;
    }
    // vからの距離が[ldist, rdist)の頂点にxを加算
    void update(int v, int ldist, int rdist, Val x){
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        int dist = (u == v ? 0 : in[v][d][3]);
        int M = nodes[u].dep_start.back().size(), L = std::max(0, ldist - dist), R = std::min(rdist - dist, M);
        if(L < R) nodes[u].sum.back().update(nodes[u].dep_start.back()[L], R == M ? nodes[u].sum.back().M : nodes[u].dep_start.back()[R], x);
        if(u != v){
          int kc = in[v][d][2];
          M = nodes[u].dep_start[kc].size(), L = std::max(0, L - 1), R = std::min(rdist - dist - 1, M);
          if(L < R) nodes[u].sum[kc].update(nodes[u].dep_start[kc][L], R == M ? nodes[u].sum[kc].M : nodes[u].dep_start[kc][R], x);
        }
      }
    }
    // vの値
    Val get(int v){
      Val ret = 0;
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        ret += nodes[u].sum.back().query(in[v][d][0], in[v][d][0] + 1);
        if(u != v) ret -= nodes[u].sum[in[v][d][2]].query(in[v][d][1], in[v][d][1] + 1);
      }
      return ret;
    }
  };
  template<typename edge, typename Val>
  struct contour_sum_weighted{
    using weight = typename edge::weight;
    int n;
    contour_sum_weighted(){}
    struct node{
      vec<accumulate1d<Val>> sum;
      vec<vec<weight>> dep_list;
    };
    std::vector<std::vector<std::array<int, 4>>> in; // [頂点番号][深さ][全体でのin / 部分木でのin / 部分木番号 / 重心との距離]
    std::vector<node> nodes;
    fast_lca fl;
    centroid_decomposition_structure_weighted<edge> cd;

    void build(tree<edge> t, const std::vector<Val> &val){
      n = t.n;
      cd = centroid_decomposition_structure_weighted<edge>(t);
      fl = fast_lca(cd.traw.g, cd.r);
      std::vector<weight> dsum(n, 0);
      
      auto dfs = [&](auto &&dfs, int cur, int par) -> void{
        for(auto &e : cd.traw[cur]){
          if(e.t == par) continue;
          dsum[e.t] = dsum[cur] + e.w;
          dfs(dfs, e.t, cur);
        }
      };
      dfs(dfs, cd.r, -1);
      auto dist_weighted = [&](int a, int b){
        return dsum[a] + dsum[b] - 2 * dsum[fl.lca(a, b)];
      };
      nodes.resize(n);
      in.resize(n);
      std::vector<Val> val2(n);
      for(int i = 0; i < n; i++){
        in[i].resize(cd.dep[i] + 1);
        val2[cd.encode(i)] = val[i];
      }
      for(int i = 0; i < n; i++){
        auto v1 = cd.shortest_raw(i);
        auto v2 = cd.shortest_raw2(i);
        nodes[i].sum.resize(v2.size() + 1);
        nodes[i].dep_list.resize(v2.size() + 1);
        int dep = cd.dep[i];
        for(int j = 0; j < v2.size(); j++){
          std::vector<Val> tmp(v2[j].size());
          for(int k = 0; k < v2[j].size(); k++){
            int u = std::get<1>(v2[j][k]);
            tmp[k] = val2[u];
            nodes[i].dep_list[j].push_back(std::get<0>(v2[j][k]));
            in[u][dep][1] = k;
            in[u][dep][2] = j;
            in[u][dep][3] = dist_weighted(u, i);
          }
          nodes[i].sum[j] = accumulate1d<Val>(tmp);
        }
        std::vector<Val> tmp(v1.size());
        for(int k = 0; k < v1.size(); k++){
          int u = std::get<1>(v1[k]);
          tmp[k] = val2[u];
          nodes[i].dep_list.back().push_back(std::get<0>(v1[k]));
          in[u][dep][0] = k;
        }
        nodes[i].sum.back() = accumulate1d<Val>(tmp);
      }
    }
    // vからの距離が[ldist, rdist)の頂点の数
    int count(int v, weight ldist, weight rdist){
      int ret = 0;
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        weight dist = (u == v ? 0 : in[v][d][3]);
        ret += std::lower_bound(nodes[u].dep_list.back().begin(), nodes[u].dep_list.back().end(), rdist - dist) - std::lower_bound(nodes[u].dep_list.back().begin(), nodes[u].dep_list.back().end(), ldist - dist);
        if(u != v){
          int kc = in[v][d][2];
          ret -= std::lower_bound(nodes[u].dep_list[kc].begin(), nodes[u].dep_list[kc].end(), rdist - dist) - std::lower_bound(nodes[u].dep_list[kc].begin(), nodes[u].dep_list[kc].end(), ldist - dist);
        }
      }
      return ret;
    }
    // vからの距離が[ldist, rdist)の頂点の和
    Val query(int v, weight ldist, weight rdist){
      Val ret = 0;
      v = cd.encode(v);
      for(int d = cd.dep[v], u = v; d >= 0; d--, u = cd.par[u]){
        weight dist = (u == v ? 0 : in[v][d][3]);
        int L = std::lower_bound(nodes[u].dep_list.back().begin(), nodes[u].dep_list.back().end(), ldist - dist) - nodes[u].dep_list.back().begin();
        int R = std::lower_bound(nodes[u].dep_list.back().begin(), nodes[u].dep_list.back().end(), rdist - dist) - nodes[u].dep_list.back().begin();
        ret += nodes[u].sum.back().query(L, R);
        if(u != v){
          int kc = in[v][d][2];
          int L = std::lower_bound(nodes[u].dep_list[kc].begin(), nodes[u].dep_list[kc].end(), ldist - dist) - nodes[u].dep_list[kc].begin();
          int R = std::lower_bound(nodes[u].dep_list[kc].begin(), nodes[u].dep_list[kc].end(), rdist - dist) - nodes[u].dep_list[kc].begin();
          ret -= nodes[u].sum[kc].query(L, R);
        }
      }
      return ret;
    }
  };
};
#endif