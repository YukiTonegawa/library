#ifndef _SPARSE_UNION_FIND_H_
#define _SPARSE_UNION_FIND_H_
#include <unordered_map>
#include <vector>
#include <algorithm>
template<typename Idx>
struct sparse_union_find{
  Idx N, cc;
  std::unordered_map<Idx, int> mp;
  std::vector<int> sz, par, rev;
  sparse_union_find(Idx N): N(N), cc(N){}
  int encode(Idx x){
    auto itr = mp.find(x);
    if(itr != mp.end()){
      return itr->second;
    }else{
      int id = mp.size();
      sz.push_back(1);
      par.push_back(id);
      rev.push_back(x);
      mp.emplace(x, id);
      return id;
    }
  }
  Idx decode(int x){
    return rev[x];
  }
  int find(int u){
    if(par[u] == u) return u;
    return par[u] = find(par[u]);
  }
  // uを含む連結成分のサイズ
  int size(int u){
    return sz[find(u)];
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
  void unite(int u, int v){
    u = find(u);
    v = find(v);
    if(u == v) return;
    cc--;
    if(sz[v] > sz[u]) std::swap(u, v);
    sz[u] += sz[v];
    par[v] = u;
  }
  // 連結成分の数
  int count_cc(){
    return cc;
  }
};
#endif