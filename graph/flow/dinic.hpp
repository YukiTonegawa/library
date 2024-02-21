#ifndef _DINIC_H_
#define _DINIC_H_
#include <vector>
#include <cassert>
#include <queue>
#include <limits>
#include <algorithm>

template<typename T>
struct Dinic{
  int n;
  struct edge{
    int to, rev;
    T cap, initial;
    edge(int to, int rev, T cap):to(to), rev(rev), cap(cap), initial(cap){}
  };
  std::vector<std::vector<edge>> G;
  std::vector<int> level, itr;
  Dinic(int n):n(n), G(n, vector<edge>()), level(n), itr(n){}

  void add_edge(int from, int to, T cap){
    G[from].push_back(edge(to, (int)G[to].size(), cap));
    G[to].push_back(edge(from, (int)G[from].size()-1, 0));
  }
  void bfs(int s){
    level.assign(n, -1);
    std::queue<int> q;
    level[s] = 0;
    q.push(s);
    while(!q.empty()){
      int v = q.front();q.pop();
      for(auto &e:G[v]){
        if(e.cap > 0 && level[e.to] < 0){
          level[e.to] = level[v] + 1;
          q.push(e.to);
        }
      }
    }
  }
  T dfs(int v, int t, T f){
    if(v==t) return f;
    for(int &i = itr[v];i<(int)G[v].size();i++){
      edge &e = G[v][i];
      if(e.cap>0&&level[v] < level[e.to]){
        T d = dfs(e.to, t, min(f, e.cap));
        if(d > 0){
          e.cap -= d;
          G[e.to][e.rev].cap += d;
          return d;
        }
      }
    }
    return 0;
  }
  T max_flow(int s, int t) {
    T ret = 0, f;
    while(bfs(s), level[t] >= 0){
      itr.assign(n,0);
      while ((f = dfs(s, t, std::numeric_limits<T>::max())) > 0) ret += f;
    }
    return ret;
  }
};
#endif