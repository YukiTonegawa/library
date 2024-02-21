#ifndef _UNION_FIND_H_
#define _UNION_FIND_H_
#include <vector>
#include <iostream>
#include <numeric>

struct union_find{
  std::vector<int> sz;
  int cc;
  union_find(int n): sz(n, -1), cc(n){}
  int find(int u){
    if(sz[u] < 0) return u;
    return sz[u] = find(sz[u]);
  }
  int size(int u){
    return -sz[find(u)];
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
  bool unite(int u, int v){
    u = find(u), v = find(v);
    if(u == v) return false;
    if(sz[v] < sz[u]) std::swap(u, v);
    cc--;
    sz[u] += sz[v];
    sz[v] = u;
    return true;
  }
  // 連結成分の数
  int count_cc(){
    return cc;
  }
};
#endif