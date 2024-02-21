#ifndef _GRAPH_ALGORITHM_CYCLE_H_
#define _GRAPH_ALGORITHM_CYCLE_H_
#include <vector>
#include <algorithm>
#include <queue>
#include "edge.hpp"
//#include "tree.hpp"
//#include "../data_structure/segment_tree/dual_segment_tree.hpp"

namespace graph_algorithm{
  template<typename T>
  using vec = std::vector<T>;
  // O(V + E)
  // ショートカットが無い閉路を1つ返す, 閉路が無い場合は空のvector
  // 自己辺はサイクルと見なさない(含めたい場合, e.from = e.toの辺を探す)
  // 多重辺があっても良い
  // 無向辺の場合辺にidをつけて逆流しないようにする
  template<typename edge>
  void cycle_noshortcut(int cur, vec<vec<edge>> &g, vec<bool> &use, vec<int> &in, vec<edge> &E, vec<edge> &res, int last_edge_id){
    if(in[cur] >= 0){
      int s = 0;
      while(E[s].from() != cur) s++;// cur->curの閉路
      std::queue<edge> shortcut;
      std::fill(use.begin(), use.end(), 0);
      while(res.empty()){
        use[cur] = true;
        int skip_in = -1, end_in = -1;
        edge skip_e = g[cur][0], end_e = g[cur][0];
        for(edge &e : g[cur]){
          int to = e.to();
          if(in[to] == -1 || (last_edge_id == e.id())) continue;
          if(in[to] > in[cur] && skip_in < in[to]) skip_in = in[to], skip_e = e;
          if(in[to] < in[cur] && use[to] && end_in < in[to]) end_in = in[to], end_e = e;
        }
        if(end_in != -1){
          int to = end_e.to();
          while(E[s].from() != to) s++;
          while(E.back().to() != cur) E.pop_back();
          while(!shortcut.empty() && in[shortcut.front().to()] <= in[to]) shortcut.pop();
          E.push_back(end_e);
          while(s != E.size()){
            int srank = (shortcut.empty() ? g.size() : in[shortcut.front().from()]);
            while(s != E.size() && in[E[s].from()] < srank) res.push_back(E[s++]);
            if(!shortcut.empty()){
              int trank = in[shortcut.front().to()];
              res.push_back(shortcut.front());
              shortcut.pop();
              while(in[E[s].from()] < trank) s++;
            }
          }
          return;
        }
        if(in[cur] + 1 != skip_in) shortcut.push(skip_e);
        cur = skip_e.to();
        last_edge_id = skip_e.id();
      }
      return;
    }
    in[cur] = E.size();
    for(edge &e : g[cur]){
      if(e.to() == cur) continue;
      if(!use[e.to()] && res.empty() && (e.id() == -1 || e.id() != last_edge_id)){
        E.push_back(e);
        cycle_noshortcut<edge>(e.to(), g, use, in, E, res, e.id());
        E.pop_back();
      }
    }
    use[cur] = true;
    in[cur] = -1;
  }
  // accept_self_loop := trueなら自己辺も長さ1のサイクルとみなす
  template<typename edge>
  vec<edge> cycle_noshortcut(vec<vec<edge>> &g, bool accept_self_loop){
    int n = g.size();
    if(accept_self_loop){
      for(int i = 0; i < n; i++){
        for(auto &e : g[i]){
          if(e.to() == e.from()){
            return {e}; // 自己辺をサイクルと見なすか
          }
        }
      }
    }
    vec<bool> use(n, false);
    vec<int> in(n, -1);
    vec<edge> res, E;
    for(int i = 0; i < n; i++){
      if(use[i]) continue;
      cycle_noshortcut<edge>(i, g, use, in, E, res, -1);
      if(!res.empty()) break;
    }
    return res;
  }

  // 全体の最小コスト閉路
  // 重みあり  O(V(V + E)logV)
  // 重みなし  O(V(V + E))
  template<typename edge>
  std::pair<typename edge::weight, vec<edge>> cycle_mincost(const vec<vec<edge>> &g, bool directed, bool weighted, bool accept_self_loop){
    using weight = typename edge::weight;
    using pdv = std::pair<weight, int>;
    static constexpr weight inf = std::numeric_limits<weight>::max() / 2;
    auto chmin = [&](weight &a, weight b)->bool{
      if(a > b){
        a = b;
        return true;
      }
      return false;
    };
    int n = g.size();
    vec<vec<edge>> cpg = g;
    for(int v = 0; v < n; v++) std::sort(cpg[v].begin(), cpg[v].end(), [](edge &a, edge &b){return a.to() < b.to();});

    std::priority_queue<pdv, vec<pdv>, std::greater<pdv>> q1;
    std::queue<pdv> q2;
    std::function<pdv()> get_top = [&]()->pdv{pdv ret = q1.top();q1.pop();return ret;};
    std::function<void(pdv)> push = [&](pdv a)->void{q1.push(a);};
    std::function<bool()> empty = [&]()->bool{return q1.empty();};

    if(!weighted){
      get_top = [&]()->pdv{pdv ret = q2.front();q2.pop();return ret;};
      push = [&](pdv a)->void{q2.push(a);};
      empty = [&]()->bool{return q2.empty();};
    }

    weight ans = inf;
    vec<edge> E;
    for(int v = n - 1; v >= 0; v--){
      vec<weight> dist(v + 1, inf);
      dist[v] = edge::z();
      push({edge::z(), v});
      vec<edge> par(v), unused;
      while(!empty()){
        auto [d, u] = get_top();
        if(dist[u] < d) continue;
        for(int i = 0; i < cpg[u].size(); i++){
          int to = cpg[u][i].to();
          if(to > v){
            cpg[u].pop_back();
            continue;
          }
          weight nxtd = d + cpg[u][i].wei();
          if(chmin(dist[to], nxtd)){
            par[to] = cpg[u][i];
            push({nxtd, to});
          }else if(par[u].id() != cpg[u][i].id()) unused.push_back(cpg[u][i]);
        }
      }
      weight pre_ans = ans;
      edge tmp_edge;
      if(!directed){
        for(edge &e : unused){
          if(dist[e.from()] == inf || dist[e.to()] == inf) continue;
          if(e.from() == e.to()){
            if(accept_self_loop && e.to() == v && chmin(ans, e.wei())){
              tmp_edge = e; // 自己辺をサイクルとみなす場合
            }
            continue;
          }
          if(chmin(ans, e.wei() + dist[e.from()] + dist[e.to()])) tmp_edge = e;
        }
      }else{
        for(edge &e : unused){
          if(dist[e.from()] == inf) continue;
          if(accept_self_loop && e.to() == v && chmin(ans, e.wei() + dist[e.from()])){
            tmp_edge = e; // 自己辺をサイクルとみなす場合
          }
          if(e.from() != v && e.to() == v && chmin(ans, e.wei() + dist[e.from()])) tmp_edge = e;
        }
      }
      if(ans < pre_ans){
        E.clear();
        int a = tmp_edge.from();
        while(a != v){
          E.push_back(par[a]);
          a = par[a].from();
        }
        std::reverse(E.begin(), E.end());
        E.push_back(tmp_edge);
        a = tmp_edge.to();
        while(a != v){
          E.push_back(par[a].reverse());
          a = par[a].from();
        }
      }
    }
    return {ans, E};
  }
  // vを含む最小コストの閉路
  // 重みあり O((V + E)logV)
  // 重みなし O(V + E)
  template<typename edge>
  std::pair<typename edge::weight, vec<edge>> cycle_mincost_specific_vertex(vec<vec<edge>> &g, int v, bool directed, bool weighted, bool accept_self_loop){
    using weight = typename edge::weight;
    using pdv = std::pair<weight, int>;
    static constexpr weight inf = std::numeric_limits<weight>::max() / 2;

    auto chmin = [&](weight &a, weight b)->bool{
      if(a > b){
        a = b;
        return true;
      }
      return false;
    };
    int n = g.size();

    std::priority_queue<pdv, vec<pdv>, std::greater<pdv>> q1;
    std::queue<pdv> q2;
    std::function<pdv()> get_top = [&]()->pdv{pdv ret = q1.top();q1.pop();return ret;};
    std::function<void(pdv)> push = [&](pdv a)->void{q1.push(a);};
    std::function<bool()> empty = [&]()->bool{return q1.empty();};

    if(!weighted){
      get_top = [&]()->pdv{pdv ret = q2.front();q2.pop();return ret;};
      push = [&](pdv a)->void{q2.push(a);};
      empty = [&]()->bool{return q2.empty();};
    }
    weight ans = inf;
    vec<edge> E;

    vec<weight> dist(n, inf);
    dist[v] = edge::z();
    push({edge::z(), v});
    vec<edge> par(n), unused;
    // 無向辺の場合lca(from, to)がvであることがvを含むサイクルになる必要十分条件
    // d1_ancestor[i] := iから辿ってどの深さ1のノードに辿り着くか
    vec<int> d1_ancestor(n, -1);

    while(!empty()){
      auto [d, u] = get_top();
      if(dist[u] < d) continue;
      if(par[u].from() == v) d1_ancestor[u] = u;
      else d1_ancestor[u] = d1_ancestor[par[u].from()];
      for(int i = 0; i < g[u].size(); i++){
        int to = g[u][i].to();
        weight nxtd = d + g[u][i].wei();
        if(chmin(dist[to], nxtd)){
          par[to] = g[u][i];
          push({nxtd, to});
        }else if(par[u].id() != g[u][i].id()) unused.push_back(g[u][i]); //無向辺の逆流対策
      }
    }
    edge tmp_edge;
    if(!directed){
      for(edge &e : unused){
        if(dist[e.from()] == inf || dist[e.to()] == inf) continue;
        if(e.from() == e.to()){
          if(accept_self_loop && e.to() == v && chmin(ans, e.wei())){
            tmp_edge = e; // 自己辺をサイクルと見なす場合
          }
          continue;
        }
        if(d1_ancestor[e.from()] == d1_ancestor[e.to()]) continue;
        if(chmin(ans, e.wei() + dist[e.from()] + dist[e.to()])) tmp_edge = e;
      }
    }else{
      for(edge &e : unused){
        if(dist[e.from()] == inf) continue;
        if(accept_self_loop && e.to() == v && chmin(ans, e.wei() + dist[e.from()])){
          tmp_edge = e; // 自己辺をサイクルと見なす場合
        }
        if(e.from() != v && e.to() == v && chmin(ans, e.wei() + dist[e.from()])) tmp_edge = e;
      }
    }
    if(ans != inf){
      int a = tmp_edge.from();
      while(a != v){
        E.push_back(par[a]);
        a = par[a].from();
      }
      std::reverse(E.begin(), E.end());
      E.push_back(tmp_edge);
      a = tmp_edge.to();
      while(a != v){
        E.push_back(par[a].reverse());
        a = par[a].from();
      }
    }
    return {ans, E};
  }
  // 全ての辺を1度だけ使うウォークが存在するか
  // オイラー路が存在するための必要十分条件
  // 無向グラフ : 奇次数の頂点が2個か0個
  // 有向グラフ : 入次数と出次数の異なる頂点が0個または(2個かつその頂点の間に辺を1本足すと0個にできる)
  // 上の条件を満たさない場合は空のvectorを返す
  // 無向グラフのときは辺にユニークなidを割り振って逆流しないようにする
  template<typename edge>
  vec<edge> eulerian_trail_directed(vec<vec<edge>> g){
    int n = g.size();
    vec<int> in(n, 0), out(n);
    for(int i = 0; i < n; i++){
      out[i] = g[i].size();
      for(auto &e : g[i]) in[e.t]++;
    }
    int s = -1;
    for(int i = 0; i < n; i++){
      if(s == -1 && out[i]) s = i;
      if(out[i] > in[i]) s = i;
    }
    std::vector<edge> res;
    auto dfs = [&](auto &&dfs, int v) -> void {
      while(out[v]){
        auto e = g[v].back();
        g[v].pop_back();
        out[v]--;
        in[e.t]--;
        dfs(dfs, e.t);
        res.push_back(e);
      }
    };
    dfs(dfs, s);
    std::reverse(res.begin(), res.end());
    for(int i = 0; i < n; i++) if(in[i] || out[i]) return {};
    for(int i = 1; i < res.size(); i++) if(res[i - 1].t != res[i].s) return {};
    return res;
  }


  template<typename edge>
  vec<edge> eulerian_trail_undirected(vec<vec<edge>> g){
    int n = g.size(), degsum = 0;
    vec<int> deg(n, 0);
    for(int i = 0; i < n; i++){
      deg[i] = g[i].size();
      degsum += deg[i];
    }
    int s = -1;
    for(int i = 0; i < n; i++){
      if(s == -1 && deg[i]) s = i;
      if(deg[i] & 1) s = i;
    }
    std::vector<bool> used_edge(degsum, 0);
    std::vector<edge> res;
    auto dfs = [&](auto &&dfs, int v) -> void {
      while(deg[v]){
        auto e = g[v].back();
        g[v].pop_back();
        while(used_edge[e.i]){
          assert(!g[v].empty()); // 同じラベルの辺が行き, 帰りの2本より多く存在すると引っかかる可能性がある
          e = g[v].back();
          g[v].pop_back();
        }
        deg[v]--;
        deg[e.t]--;
        used_edge[e.i] = 1;
        dfs(dfs, e.t);
        res.push_back(e);
      }
    };
    dfs(dfs, s);
    std::reverse(res.begin(), res.end());
    for(int i = 0; i < n; i++) if(deg[i]) return {};
    for(int i = 1; i < res.size(); i++) if(res[i - 1].t != res[i].s) return {};
    return res;
  }
  
  /*
  // ! グラフが連結でない場合

  // res[i] := iを含む閉路 O(V^2 + Elog^2V)
  // 全ての辺に[0, E)のuniqueな値を振っておく
  template<typename edge>
  vec<vec<edge>> undirected_cycle_detection_vv(vec<vec<edge>> &g){
    int n = g.size();
    tree<edge> t = bfs_tree<edge>(g, 0);
    vec<bool> used_edge{0};
    for(int i = 0; i < n; i++){
      for(edge &e : t.g[i]){
        assert(0 <= e.id() && e.id() <= 1e7);
        while(used_edge.size() <= e.id()) used_edge.resize(used_edge.size() << 1);
        used_edge[e.id()] = true;
      }
    }
    t.hld_build();
    t.bfs_build();
    dual_segment_tree<range_set_get_any, edge*, edge*> seg(n, nullptr);
    for(int i = 0; i < n; i++){
      for(edge &e : g[i]){
        if(e.id() >= used_edge.size() || !used_edge[e.id()]){// bfs木に含まれない
          for(auto [l, r] : t.hld_p->unordered_path(e.from(), e.to())) seg.update(l, r, &e);
          if(e.id() < used_edge.size()) used_edge[e.id()] = true;
        }
      }
    }
    vec<vec<edge>> res(n);
    for(int i = 0; i < n; i++){
      if(!res[i].empty()) continue;
      edge *e = seg.get(t.hld_p->index_vertex(i));
      if(!e) continue;// iを含む閉路が無い
      vec<edge> E;
      int a = e->from(), b = e->to(), l = t.hld_p->lca(a, b);
      while(a != l){
        E.push_back(t.bfs_p->get_parent_edge(a));
        a = t.parent[a];
      }
      std::reverse(E.begin(), E.end());
      E.push_back(*e);
      while(b != l){
        E.push_back(t.bfs_p->get_parent_edge(b).reverse());
        b = t.parent[b];
      }
      for(edge &e : E) if(res[e.from()].empty()) res[e.from()] = E;
    }
    return res;
  }
  */
};
#endif