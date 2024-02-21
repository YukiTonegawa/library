#ifndef _CONVERT_DIRECTED_FUNCTIONAL_GRAPH_H_
#define _CONVERT_DIRECTED_FUNCTIONAL_GRAPH_H_

#include <vector>
#include <cassert>
#include <queue>

// n頂点n辺の無向グラフをfunctional graphに変換
std::vector<int> convert_directed_functional_graph(const std::vector<std::vector<int>> &g){
  int n = g.size();
  std::vector<int> next(n, -1);
  std::vector<int> deg(n);
  std::vector<bool> used(n, false);
  std::queue<int> q;
  for(int i = 0; i < n ; i++){
    deg[i] = g[i].size();
    if(deg[i] == 1) q.push(i);
  }
  // 足部分を処理
  while(!q.empty()){
    int v = q.front();
    q.pop();
    assert(deg[v] == 1);
    deg[v] = 0;
    used[v] = 1;
    for(int to : g[v]){
      if(used[to]) continue;
      next[v] = to;
      deg[to]--;
      if(deg[to] == 1) q.push(to);
    }
  }
  auto f = [&](auto &&f, int cur, int par){
    if(next[cur] != -1) return;
    for(int to : g[cur]){
      if(deg[to] < 2 || to == par) continue;
      next[cur] = to;
      f(f, to, cur);
      return;
    }
    next[cur] = par;
  };
  // ループ部分を処理
  for(int i = 0; i < n; i++){
    if(deg[i] >= 2){
      f(f, i, -1);
    }
  }
  return next;
}
#endif