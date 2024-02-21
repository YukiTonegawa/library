#ifndef _INCREMENTAL_TREE_H_
#define _INCREMENTAL_TREE_H_
#include <vector>
#include <array>
#include <algorithm>
struct incremental_tree{
  static constexpr int max_nodecount_log = 20;
  static constexpr int max_nodecount_log_half = max_nodecount_log / 2;
  std::vector<int> dep;
  std::vector<std::array<int, 10>> par;
  // 0を根とする1頂点の木を作る 
  incremental_tree(){
    dep.push_back(0);
    par.push_back({});
    par.back().fill(0);
  }
  int size(){
    return dep.size();
  }
  // pを親とする頂点を作る(頂点番号を返す)
  int make_node(int p){
    int id = dep.size();
    dep.push_back(dep[p] + 1);
    par.push_back({});
    par[id][0] = p;
    for(int i = 0; i + 1 < max_nodecount_log_half; i++) par[id][i + 1] = par[par[par[par[id][i]][i]][i]][i];
    return id;
  }
  // 親(0の親は0とする)
  int parent(int v){
    return par[v][0];
  }
  int depth(int v){
    return dep[v];
  }
  // depth(v) < kの場合0
  int la(int v, int k){
    for(int i = max_nodecount_log_half - 1; i >= 0; i--) while(k >= (1 << (i << 1))) v = par[v][i], k -= 1 << (i << 1);
    return v;
  }
  int lca(int u, int v){
    int depdiff = depth(u) - depth(v);
    if(depdiff < 0){
      std::swap(u, v);
      depdiff *= -1;
    }
    u = la(u, depdiff);
    if(u == v) return u;
    for(int i = max_nodecount_log_half - 1; i >= 0 && u != v; i--) while(par[u][i] != par[v][i]) u = par[u][i], v = par[v][i];
    return par[u][0];
  }
  int dist(int u, int v){
    int l = lca(u, v);
    return dep[u] + dep[v] - 2 * dep[l];
  }
  // u->v単純パスのk番目
  // dist(u, v) < kの場合-1
  int jump_on_tree(int u, int v, int k){
    int l = lca(u, v), dlu = dep[u] - dep[l];
    if(dlu > k) return la(u, k);
    k = dep[v] - dep[l] - k + dlu;
    if(k < 0) return -1;
    return la(v, k);
  }
};
/*
struct incremental_tree_centroid{
  static constexpr int max_nodecount_log = 20;
  static constexpr int max_nodecount_log_half = max_nodecount_log / 2;
  std::vector<int> dep;
  std::vector<std::array<int, 10>> par;
  // 0を根とする1頂点の木を作る 
  incremental_tree(){
    dep.push_back(0);
    par.push_back({});
    par.back().fill(0);
  }
  int size(){
    return dep.size();
  }
  // pを親とする頂点を作る(頂点番号を返す)
  int make_node(int p){
    int id = dep.size();
    dep.push_back(dep[p] + 1);
    par.push_back({});
    par[id][0] = p;
    for(int i = 0; i + 1 < max_nodecount_log_half; i++) par[id][i + 1] = par[par[par[par[id][i]][i]][i]][i];
    return id;
  }
  // 親(0の親は0とする)
  int parent(int v){
    return par[v][0];
  }
  int depth(int v){
    return dep[v];
  }
  // depth(v) < kの場合0
  int la(int v, int k){
    for(int i = max_nodecount_log_half - 1; i >= 0; i--) while(k >= (1 << (i << 1))) v = par[v][i], k -= 1 << (i << 1);
    return v;
  }
  int lca(int u, int v){
    int depdiff = depth(u) - depth(v);
    if(depdiff < 0){
      std::swap(u, v);
      depdiff *= -1;
    }
    u = la(u, depdiff);
    if(u == v) return u;
    for(int i = max_nodecount_log_half - 1; i >= 0 && u != v; i--) while(par[u][i] != par[v][i]) u = par[u][i], v = par[v][i];
    return par[u][0];
  }
  int dist(int u, int v){
    int l = lca(u, v);
    return dep[u] + dep[v] - 2 * dep[l];
  }
  // u->v単純パスのk番目
  // dist(u, v) < kの場合-1
  int jump_on_tree(int u, int v, int k){
    int l = lca(u, v), dlu = dep[u] - dep[l];
    if(dlu > k) return la(u, k);
    k = dep[v] - dep[l] - k + dlu;
    if(k < 0) return -1;
    return la(v, k);
  }
};
*/
#endif