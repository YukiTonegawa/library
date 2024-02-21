#ifndef _BIPARTITE_MATCHING_H_
#define _BIPARTITE_MATCHING_H_
#include "max_flow.hpp"

struct HopcroftKarp {
  std::vector<std::vector<int>> g;
  std::vector<int> dist, match;
  std::vector<bool> used, vv;
  HopcroftKarp(int n, int m) : g(n), match(m, -1), used(n){}
  void add_edge(int u, int v) {
    g[u].push_back(v);
  }
  void bfs() {
    dist.assign(g.size(), -1);
    simple_queue<int> que;
    for(int i = 0; i < g.size(); i++) {
      if(!used[i]) {
        que.push(i);
        dist[i] = 0;
      }
    }
    while(!que.empty()) {
      int a = que.front();
      que.pop();
      for(auto &b : g[a]){
        int c = match[b];
        if(c >= 0 && dist[c] == -1) {
          dist[c] = dist[a] + 1;
          que.push(c);
        }
      }
    }
  }
  bool dfs(int a) {
    vv[a] = true;
    for(auto &b : g[a]) {
      int c = match[b];
      if(c < 0 || (!vv[c] && dist[c] == dist[a] + 1 && dfs(c))) {
        match[b] = a;
        used[a] = true;
        return (true);
      }
    }
    return (false);
  }
  int bipartite_matching(){
    int ret = 0;
    while(true) {
      bfs();
      vv.assign(g.size(), false);
      int flow = 0;
      for(int i = 0; i < g.size(); i++) {
        if(!used[i] && dfs(i)) ++flow;
      }
      if(flow == 0) return (ret);
      ret += flow;
    }
  }
};
struct bipartite_matching{
  int n, m, s, t, edgenum;
  HopcroftKarp g;
  bipartite_matching(int n, int m): n(n), m(m), s(n + m), t(s + 1), edgenum(0), g(n, m){}
  // 左側のiと右側のjに辺を張る　
  void add_edge(int i, int j){
    assert(0 <= i && i < n);
    assert(0 <= j && j < m);
    g.add_edge(i, j);
  }
  // {最大流, v} v[i] := 左側のiに対する右側のペア　　存在しない場合は-1
  std::pair<int, std::vector<int>> flow(){
    int f = g.bipartite_matching();
    std::vector<int> p(n, -1);
    for(int i = 0; i < m; i++){
      if(g.match[i] != -1) p[g.match[i]] = i;
    }
    return {f, g.match};
  }
  // {最大流, E}
  std::pair<int, std::vector<std::pair<int, int>>> flow_edge(){
    int f = g.bipartite_matching();
    std::vector<std::pair<int, int>> E;
    for(int i = 0; i < m; i++){
      if(g.match[i] != -1) E.push_back({g.match[i], i});
    }
    return {f, E};
  }
};
#endif
