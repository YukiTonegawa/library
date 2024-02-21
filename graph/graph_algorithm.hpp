#ifndef _GRAPH_ALGORITHM_H_
#define _GRAPH_ALGORITHM_H_
#include <vector>
#include <queue>
#include <numeric>
#include <algorithm>
#include <limits>
#include "edge.hpp"

namespace graph_algorithm{
  template<typename T>
  using vec = std::vector<T>;
  // O(V + E)
  // 辺の重みが1
  template<typename edge>
  struct bfs_shortest_path{
  private:
    using weight = typename edge::weight;
    using dist_p = std::pair<weight, int>;
    vec<vec<edge>> &g;
  public:
    bfs_shortest_path(vec<vec<edge>> &g): g(g){}
    static constexpr weight inf = std::numeric_limits<weight>::max() / 2;
    static constexpr weight minf = std::numeric_limits<weight>::min() / 2;
    vec<weight> dist;
    vec<edge> par;
    void build(int s){
      int n = g.size();
      if(dist.empty()){
        dist.resize(n, inf);
        par.resize(n, edge{});
      }else{
        std::fill(dist.begin(), dist.end(), inf);
        std::fill(par.begin(), par.end(), edge{});
      }
      std::queue<dist_p> que;
      dist[s] = edge::z();
      que.push(dist_p(edge::z(), s));
      while(!que.empty()){
        auto [w, v] = que.front();
        que.pop();
        if(dist[v] < w) continue;
        for(edge &e: g[v]){
          assert(e.wei() == 1);
          weight d = dist[v] + e.wei();
          int to = e.to();
          if(dist[to] > d){
            dist[to] = d;
            par[to] = e;
            que.push(dist_p(d, to));
          }
        }
      }
    }
    vec<edge> get_path(int v){
      assert(!dist.empty());
      vec<edge> ret;
      while(par[v].from() != -1) ret.push_back(par[v]), v = par[v].from();
      std::reverse(ret.begin(), ret.end());
      return ret;
    }
    weight operator [](int v){return dist[v];}
  };

  // O(V + E)
  // 辺の重みが0か1
  template<typename edge>
  struct zero_one_bfs_shortest_path{
  private:
    using weight = typename edge::weight;
    vec<vec<edge>> &g;
  public:
    zero_one_bfs_shortest_path(vec<vec<edge>> &g): g(g){}
    static constexpr weight inf = std::numeric_limits<weight>::max() / 2;
    static constexpr weight minf = std::numeric_limits<weight>::min() / 2;
    vec<weight> dist;
    vec<edge> par;
    void build(int s){
      int n = g.size();
      if(dist.empty()){
        dist.resize(n, inf);
        par.resize(n, edge{});
      }else{
        std::fill(dist.begin(), dist.end(), inf);
        std::fill(par.begin(), par.end(), edge{});
      }
      std::queue<int> que0, que1;
      dist[s] = edge::z();
      weight dcur = dist[s];
      que0.push(s);
      while(!que0.empty() || !que1.empty()){
        if(que0.empty()){
          std::swap(que0, que1);
          dcur++;
        }
        int v = que0.front();
        que0.pop();
        if(dist[v] < dcur) continue;
        for(edge &e: g[v]){
          weight w = e.wei();
          assert(w == 0 || w == 1);
          weight d = dist[v] + w;
          if(dist[e.t] > d){
            dist[e.t] = d;
            par[e.t] = e;
            if(w == 0) que0.push(e.t);
            else que1.push(e.t);
          }
        }
      }
    }
    vec<edge> get_path(int v){
      assert(!dist.empty());
      vec<edge> ret;
      while(par[v].from() != -1) ret.push_back(par[v]), v = par[v].from();
      std::reverse(ret.begin(), ret.end());
      return ret;
    }
    weight operator [](int v){return dist[v];}
  };

  // O((V + E)logV)
  // 辺の重みが非負(負の閉路がなければ一応動く)
  template<typename edge>
  struct dijkstra{
  private:
    using weight = typename edge::weight;
    using dist_p = std::pair<weight, int>;
    vec<vec<edge>> &g;
  public:
    dijkstra(vec<vec<edge>> &g): g(g){}
    static constexpr weight inf = std::numeric_limits<weight>::max() / 2;
    static constexpr weight minf = std::numeric_limits<weight>::min() / 2;
    vec<weight> dist;
    vec<edge> par;
    void build(int s){
      int n = g.size();
      if(dist.empty()){
        dist.resize(n, inf);
        par.resize(n, edge{});
      }else{
        std::fill(dist.begin(), dist.end(), inf);
        std::fill(par.begin(), par.end(), edge{});
      }
      std::priority_queue<dist_p, vec<dist_p>, std::greater<dist_p>> que;
      dist[s] = edge::z();
      que.push(dist_p(edge::z(), s));
      while(!que.empty()){
        auto [w, v] = que.top();
        que.pop();
        if(dist[v] < w) continue;
        for(edge &e: g[v]){
          weight d = dist[v] + e.wei();
          int to = e.to();
          if(dist[to] > d){
            dist[to] = d;
            par[to] = e;
            que.push(dist_p(d, to));
          }
        }
      }
    }
    vec<edge> get_path(int v){
      assert(!dist.empty());
      vec<edge> ret;
      while(par[v].from() != -1) ret.push_back(par[v]), v = par[v].from();
      std::reverse(ret.begin(), ret.end());
      return ret;
    }
    weight operator [](int v){return dist[v];}
  };

  // O(VE)
  // inf: 到達不可, minf: 負の閉路
  template<typename edge>
  struct bellman_ford{
  private:
    using weight = typename edge::weight;
    using dist_p = std::pair<weight, int>;
    vec<vec<edge>> &g;
  public:
    bellman_ford(vec<vec<edge>> &g): g(g){}
    static constexpr weight inf = std::numeric_limits<weight>::max() / 2;
    static constexpr weight minf = std::numeric_limits<weight>::min() / 2;
    vec<weight> dist;
    vec<edge> par;

    void build(int s){
      int n = g.size();
      if(dist.empty()){
        dist.resize(n, inf);
        par.resize(n);
      }else{
        std::fill(dist.begin(), dist.end(), inf);
        std::fill(par.begin(), par.end(), edge{});
      }
      dist[s] = edge::z();
      for(int lp = 0; ; lp++){
        bool update = false;
        for(int i = 0; i < n; i++){
          if(dist[i] == inf) continue;
          for(edge e : g[i]){
            weight &dto = dist[e.to()];
            if(dto == minf){
              if(dto != minf) update = true;
              dto = minf;
            }else if(dto == inf || dto > dist[i] + e.wei()){
              dto = (lp > n ? minf : dist[i] + e.wei());
              par[e.to()] = e;
              update = true;
            }
          }
        }
        if(!update) break;
      }
    }
    vec<edge> get_path(int v){
      assert(!dist.empty());
      vec<edge> ret;
      while(par[v].from() != -1) ret.push_back(par[v]), v = par[v].from();
      std::reverse(ret.begin(), ret.end());
      return ret;
    }
    weight operator [](int v){return dist[v];}
  };

  // O(V^3)
  template<typename edge>
  struct warshall_floyd{
  private:
    using weight = typename edge::weight;
    vec<vec<edge>> &g;
  public:
    warshall_floyd(vec<vec<edge>> &g): g(g){}
    static constexpr weight inf = std::numeric_limits<weight>::max() / 2;
    static constexpr weight minf = std::numeric_limits<weight>::min() / 2;
    vec<vec<weight>> dist;
    void build(){
      int n = g.size();
      dist.resize(n, vec<weight>(n, inf));
      for(int i = 0; i < n; i++){
        dist[i][i] = 0;
        for(edge &e : g[i]){
          dist[i][e.to()] = std::min(dist[i][e.to()], e.wei());
        }
      }
      for(int k = 0; k < n; k++){
        for(int s = 0; s < n; s++){
          for(int t = 0; t < n; t++){
            dist[s][t] = std::min(dist[s][t], dist[s][k] + dist[k][t]);
          }
        }
      }
    }
    vec<weight>& operator [](int v){return dist[v];}
  };
};

namespace graph_algorithm{
  // {連結成分, DAG}
  template<typename edge>
  std::pair<vec<int>, vec<vec<int>>> scc(vec<vec<edge>> &g){
    int n = g.size();
    vec<int> v(n), cmp(n, 0);
    vec<vec<int>> rg(n), V;
    auto scc_dfs = [&](auto &&scc_dfs, int cur, int &sz)->void{
      cmp[cur] = -1;
      for(edge &e : g[cur]){
        int to = e.to();
        rg[to].push_back(cur);
        if(cmp[to] == 0) scc_dfs(scc_dfs, to, sz);
      }
      v[sz++] = cur;
    };
    auto scc_rdfs = [&](auto &&scc_rdfs, int cur, const int k)->void{
      cmp[cur] = k;
      V[k].push_back(cur);
      for(int to : rg[cur]) if(cmp[to] == -1) scc_rdfs(scc_rdfs, to, k);
    };
    for(int i = 0, j = 0; i < n; i++) if(!cmp[i]) scc_dfs(scc_dfs, i, j);
    for(int i = (int)v.size() - 1, j = 0; i >= 0; i--){
      if(cmp[v[i]] == -1){
        V.push_back(vec<int>());
        scc_rdfs(scc_rdfs, v[i], j++);
      }
    }
    return {cmp, V};
  }
  // {連結成分, 森}
  template<typename edge>
  std::pair<vec<int>, vec<vec<int>>> two_edge_connected(vec<vec<edge>> &g){
    int n = g.size();
    vec<int> v(n), cmp(n, 0);
    vec<vec<int>> V;
    vec<vec<bool>> edge_used(n);
    auto tec_dfs = [&](auto &&tec_dfs, int cur, int &sz)->void{
      cmp[cur] = -1;
      for(int i = 0; i < g[cur].size(); i++){
        int to = g[cur][i].to();
        if(cmp[to] == 0) edge_used[cur][i] = true, tec_dfs(tec_dfs, to, sz);
      }
      v[sz++] = cur;
    };
    auto tec_rdfs = [&](auto &&tec_rdfs, int cur, const int k)->void{
      cmp[cur] = k;
      V[k].push_back(cur);
      for(int i = 0; i < g[cur].size(); i++){
        int to = g[cur][i].to();
        if(cmp[to] == -1 && !edge_used[cur][i]) tec_rdfs(tec_rdfs, to, k);
      }
    };
    for(int i = 0; i < n; i++) edge_used[i].resize(g[i].size(), 0);
    for(int i = 0, j = 0; i < n; i++) if(!cmp[i]) tec_dfs(tec_dfs, i, j);
    for(int i = (int)v.size() - 1, j = 0; i >= 0; i--){
      if(cmp[v[i]] == -1){
        V.push_back(vec<int>());
        tec_rdfs(tec_rdfs, v[i], j++);
      }
    }
    return {cmp, V};
  }
  // 二重頂点連結成分分解
  // {間接点フラグ, 各連結成分が含む頂点}
  template<typename edge>
  std::pair<vec<bool>, vec<vec<int>>> bcc(vec<vec<edge>> &g){
    int n = g.size();
    vec<vec<int>> V;
    vec<int> child(n, 0), dep(n, -1), low(n);
    vec<bool> used(n, false), is_articulation(n, false);
    vec<edge> tmp_edge;
    auto bcc_dfs = [&](auto &&bcc_dfs, int cur, int par, int d)->void{
      if(par != -1) child[par]++;
      dep[cur] = low[cur] = d;
      for(edge &e : g[cur]){
        int to = e.to();
        if(to == par) continue;
        if(dep[to] < dep[cur]) tmp_edge.push_back(e);
        if(dep[e.to()] == -1){
          bcc_dfs(bcc_dfs, to, cur, d + 1);
          if(low[to] >= dep[cur]){
            is_articulation[cur] = true;
            V.push_back(vec<int>());
            bool is_ok = false;
            while(!tmp_edge.empty() && !is_ok){
              edge e = tmp_edge.back();
              tmp_edge.pop_back();
              if(e.from() == cur && e.to() == to) is_ok = true;
              if(!used[e.to()]) V.back().push_back(e.to()), used[e.to()] = true;
              if(!used[e.from()]) V.back().push_back(e.from()), used[e.from()] = true;
            }
            for(int v : V.back()) used[v] = false;
          }
          low[cur] = std::min(low[cur], low[to]);
        }else low[cur] = std::min(low[cur], dep[to]);
      }
    };
    for(int i = 0; i < n; i++){
      if(dep[i] != -1) continue;
      int vsz_pre = V.size();
      bcc_dfs(bcc_dfs, i, -1, 0);
      is_articulation[i] = (child[i] > 1);
      if(V.size() == vsz_pre) V.push_back(vec<int>{i});// 孤立点
    }
    return {is_articulation, V};
  }
  template<typename edge>
  std::tuple<vec<bool>, vec<int>, vec<vec<simple_edge<int>>>> block_cut_tree(vec<vec<edge>> &g){
    auto [is_articulation, V] = bcc<edge>(g);
    int n = g.size();
    vec<int> cmp(n, -1);
    int m = V.size(), a = m;
    for(int i = 0; i < m; i++){
      for(int v : V[i]){
        if(is_articulation[v]) cmp[v] = a++;
        else cmp[v] = i;
      }
    }
    vec<vec<simple_edge<int>>> G(a);
    for(int i = 0; i < m; i++){
      for(int v : V[i]){
        if(is_articulation[v]){
          G[i].push_back({i, cmp[v]});
          G[cmp[v]].push_back({cmp[v], i});
        }
      }
    }
    return {is_articulation, cmp, G};
  }
};
namespace graph_algorithm{
  // 終了時にinが0でない要素がある -> 閉路が存在する
  // 閉路があるなら空のvectorを返す
  template<typename edge>
  vec<int> topological_sort(vec<vec<edge>> &g){
    int n = g.size();
    std::queue<int> que;
    vec<int> in(n, 0), ret;
    for(int i = 0; i < n; i++) for(edge e : g[i]) in[e.to()]++;
    for(int i = 0; i < n; i++) if(!in[i]) que.push(i);
    while(!que.empty()){
      int p = que.front();
      que.pop();
      ret.push_back(p);
      for(edge &e : g[p]){
        int to = e.to();
        if(!(--in[to])) que.push(to);
      }
    }
    for(int i = 0; i < n; i++) if(in[i] != 0) return {};
    return ret;
  }

  // プリム法, 連結なら始点sは関係ない
  template<typename edge>
  vec<edge> undirected_mst(vec<vec<edge>> &g, int s = 0){
    int n = g.size();
    assert(s < n);
    static vec<bool> V(n, 0);
    vec<edge> ret;
    using pde = std::pair<typename edge::weight, edge>;
    std::priority_queue<pde, vec<pde>, std::function<bool(pde, pde)>> que([](pde a, pde b){
      return a.first > b.first;
    });
    V[s] = true;
    for(edge &e : g[s]) que.push(pde{e.wei(), e});
    while(!que.empty()){
      auto [d, e] = que.top();
      que.pop();
      if(V[e.to()]) continue;
      V[e.to()] = true;
      ret.push_back(e);
      for(edge &ec : g[e.to()]) if(!V[ec.to()]) que.push({ec.wei(), ec});
    }
    for(edge &e : ret) V[e.to()] = V[e.from()] = false;
    return ret;
  }
  /*
  // プリム法より早い
  template<typename edge>
  vec<edge> undirected_mst_kruskal(vec<vec<edge>> &g){
    int n = g.size();
    union_find uf(n);
    vec<edge> E, res;
    for(int i = 0; i < n; i++){
      for(auto e : g[i]){
        if(i < e.to()) E.push_back(e);
      }
    }
    std::sort(E.begin(), E.end(), [](edge &a, edge &b){return a.wei() < b.wei();});
    for(auto e : E){
      if(!uf.same(e.from(), e.to())){
        res.push_back(e);
        uf.unite(e.from(), e.to());
      }
    }
    return res;
  }
  template<typename edge>
  vec<edge> undirected_mst_kruskal(int n, vec<edge> E){
    union_find uf(n);
    vec<edge> res;
    std::sort(E.begin(), E.end(), [](edge &a, edge &b){return a.wei() < b.wei();});
    for(auto e : E){
      if(!uf.same(e.from(), e.to())){
        res.push_back(e);
        uf.unite(e.from(), e.to());
      }
    }
    return res;
  }
  */
  // rを根とするbfs木O(V + E)
  template<typename edge>
  vec<vec<edge>> bfs_tree(vec<vec<edge>> &g, int r){
    int n = g.size();
    std::queue<int> que;
    vec<bool> used(n, false);
    que.push(r);
    vec<vec<edge>> ret(n);
    used[r] = true;
    while(!que.empty()){
      int v = que.front();
      que.pop();
      for(edge &e : g[v]){
        int to = e.to();
        if(used[to]) continue;
        used[to] = true;
        ret[v].push_back(e);
        que.push(to);
      }
    }
    return ret;
  }
  // rを根とするbfs木, 最短経路的 O((V + E)logV)
  // {木, 重みのテーブル}
  template<typename edge>
  std::pair<vec<vec<edge>>, vec<typename edge::weight>> bfs_tree_shortest(vec<vec<edge>> &g, int r){
    int n = g.size();
    using weight = typename edge::weight;
    using pdv = std::pair<weight, int>;
    static constexpr weight inf = std::numeric_limits<weight>::max() / 2;
    std::priority_queue<pdv, vec<pdv>, std::greater<pdv>> que;
    vec<weight> dist(n, inf);
    dist[r] = edge::z();
    que.push({edge::z(), r});

    vec<vec<edge>> pedge(n);
    vec<vec<edge>> ret(n);

    while(!que.empty()){
      auto [d, v] = que.top();
      que.pop();
      if(dist[v] < d) continue;
      for(edge &e : g[v]){
        int to = e.to();
        weight nxtd = d + e.wei();
        if(dist[to] > nxtd){
          dist[to] = nxtd;
          if(!pedge[to].empty()) pedge[to].pop_back();
          pedge[to].push_back(e);
          que.push({nxtd, to});
        }
      }
    }
    for(int i = 0; i < n; i++){
      if(!pedge[i].empty()){
        edge e = pedge[i][0];
        ret[e.s].push_back(e);
      }
    }
    return {ret, dist};
  }
  // g[i]の辺を{同じcmpへの辺, 異なるcmpへの辺}に並び替える, O(V + E)
  template<typename edge>
  void cmp_edge_arrange(const vec<int> &cmp, vec<vec<edge>> &g){
    int n = g.size();
    for(int i = 0; i < n; i++){
      int m = g[i].size();
      int l = 0, r = m - 1;
      while(l < r){
        while(l < m && cmp[i] == cmp[g[i][l].to()]) l++;
        while(0 < r && cmp[i] == cmp[g[i][r].to()]) r--;
        if(l < r) std::swap(g[i][l], g[i][r]);
      }
    }
  }
};
#endif
