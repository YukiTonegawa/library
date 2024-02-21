#ifndef ROLLBACK_UNION_FIND_H_
#define ROLLBACK_UNION_FIND_H_
#include <vector>
#include <tuple>
#include <numeric>
#include <cassert>
#include <queue>

struct rollback_union_find{
  int n, pos, cc;
  std::vector<int> par, sz, history;
  rollback_union_find(int n = 0): n(n), pos(0), cc(n), par(n, -1), sz(n, 1), history(2 * n){}
  int find(int u){
    while(par[u] != -1) u = par[u];
    return u;
  }
  int size(int u){
    return sz[find(u)];
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
  bool unite(int u, int v){
    u = find(u), v = find(v);
    if(history.size() <= pos + 1) history.resize(2 * (pos + 1));
    if(u == v){
      history[pos++] = n + u;
      return false;
    }
    if(sz[v] > sz[u]) std::swap(u, v);
    cc--;
    sz[u] += sz[v];
    par[v] = u;
    history[pos++] = u;
    history[pos++] = v;
    return true;
  }
  // 戻り値: {橋だったか, 頂点c, 頂点p}
  std::tuple<int, int> rollback(){
    assert(pos);
    int c = history[--pos];
    if(c >= n) return {c - n, c - n};
    int p = history[--pos];
    cc++;
    par[c] = -1;
    sz[p] -= sz[c];
    return {c, p};
  }
  // 連結成分の数
  int count_cc(){
    return cc;
  }
};

struct rollback_union_find_enumerate{
  int n, pos, cc;
  std::vector<int> par, sz, history;
  std::vector<std::vector<int>> ch;
  rollback_union_find_enumerate(int n): n(n), pos(0), cc(n), par(n, -1), sz(n, 1), history(2 * n), ch(n){}
  int find(int u){
    while(par[u] != -1) u = par[u];
    return u;
  }
  int size(int u){
    return sz[find(u)];
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
  bool unite(int u, int v){
    u = find(u), v = find(v);
    if(history.size() <= pos + 1) history.resize(2 * (pos + 1));
    if(u == v){
      history[pos++] = n + u;
      return false;
    }
    if(sz[v] > sz[u]) std::swap(u, v);
    cc--;
    sz[u] += sz[v];
    par[v] = u;
    ch[u].push_back(v);
    history[pos++] = u;
    history[pos++] = v;
    return true;
  }
  // 戻り値: {橋だったか, 頂点c, 頂点p}
  std::tuple<int, int> rollback(){
    assert(pos);
    int c = history[--pos];
    if(c >= n) return {c - n, c - n};
    int p = history[--pos];
    cc++;
    par[c] = -1;
    sz[p] -= sz[c];
    ch[p].pop_back();
    return {c, p};
  }
  // 連結成分の数
  int count_cc(){
    return cc;
  }
  // uを含む連結成分の全要素
  std::vector<int> enumerate(int u){
    std::vector<int> ret;
    std::queue<int> q;
    q.push(find(u));
    while(!q.empty()){
      int v = q.front();
      q.pop();
      ret.push_back(v);
      for(int c:ch[v]) q.push(c);
    }
    return ret;
  }
};
#endif