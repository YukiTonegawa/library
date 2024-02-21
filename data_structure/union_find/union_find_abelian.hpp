#ifndef _UNION_FIND_COMMUTATIVE_H_
#define _UNION_FIND_COMMUTATIVE_H_
#include "../../algebraic_structure/abelian_group.hpp"
#include <numeric>
#include <vector>
#include <stack>
#include <queue>
#include <iostream>

// rollbackなし経路圧縮
template<typename abelian_group>
struct union_find_abelian{
  using Val = typename abelian_group::Val;
  using Lazy = typename abelian_group::Lazy;
  static constexpr auto id        = abelian_group::id;
  static constexpr auto id_lazy   = abelian_group::id_lazy;
  static constexpr auto inv       = abelian_group::inv;
  static constexpr auto inv_lazy  = abelian_group::inv_lazy;
  static constexpr auto propagate = abelian_group::propagate;
  static constexpr auto merge     = abelian_group::merge;
  static constexpr auto apply     = abelian_group::apply;
private:
  int n;
  std::vector<int> par, sz;
  std::vector<Val> val, sum;
  std::vector<Lazy> lazy;
public:
  union_find_abelian(){}
  union_find_abelian(int n): n(n), par(n, -1), sz(n, 1), val(n, id()), sum(val), lazy(n, id_lazy()){}
  union_find_abelian(const std::vector<Val> &v): n(v.size()), par(n, -1), sz(n, 1), val(v), sum(v), lazy(n, id_lazy()){}
  
  int find(int u){
    if(par[u] == -1) return u;
    int p = find(par[u]);
    if(p != par[u]) lazy[u] = propagate(lazy[u], lazy[par[u]]);
    return par[u] = p;
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
  void unite(int u, int v){
    u = find(u), v = find(v);
    if(u == v) return;
    if(sz[v] > sz[u]) std::swap(u, v);
    lazy[v] = propagate(lazy[v], inv_lazy(lazy[u]));
    sum[u] = merge(sum[u], sum[v]);
    par[v] = u;
    sz[u] += sz[v];
  }
  // aを含む連結成分のサイズ
  int size(int a){
    return sz[find(a)];
  }
  // uの値をxに更新
  void update(int u, Val x){
    sum[u] = merge(sum[u], x);
    val[u] = merge(val[u], x);
    while(par[u] != -1){
      sum[par[u]] = merge(sum[par[u]], x);
      u = par[u];
    }
  }
  // uの値
  Val get(int u){
    Val res = val[u];
    Lazy lz = lazy[u];
    while(par[u] != -1){
      lz = propagate(lz, lazy[par[u]]);
      u = par[u];
    }
    return apply(res, lz, 1);
  }
  void set(int u, Val x){
    update(u, merge(x, inv(get(u))));
  }
  // uを含む連結成分の和
  Val query_cc(int u){
    u = find(u);
    return apply(sum[u], lazy[u], 0, sz[u]);
  }
  // uを含む連結成分にxを作用させる
  void update_cc(int u, Lazy x){
    u = find(u);
    lazy[u] = propagate(lazy[u], x);
  }
};
template<typename abelian_group>
struct union_find_abelian_rollback{
  using Val = typename abelian_group::Val;
  using Lazy = typename abelian_group::Lazy;
  static constexpr auto id        = abelian_group::id;
  static constexpr auto id_lazy   = abelian_group::id_lazy;
  static constexpr auto inv       = abelian_group::inv;
  static constexpr auto inv_lazy  = abelian_group::inv_lazy;
  static constexpr auto propagate = abelian_group::propagate;
  static constexpr auto merge     = abelian_group::merge;
  static constexpr auto apply     = abelian_group::apply;
private:
  int n, cc;
  std::vector<int> par, sz, dup;
  std::vector<std::vector<int>> ch;
  std::vector<Val> val, sum;
  std::vector<Lazy> lazy;
  std::vector<int> history;
public:
  union_find_abelian_rollback(){}
  union_find_abelian_rollback(int n): n(n), cc(n), par(n, -1), sz(n, 1), dup(n, 0), ch(n), val(n, id()), sum(val), lazy(n, id_lazy()){}
  union_find_abelian_rollback(const std::vector<Val> &v): n(v.size()), cc(n), par(n, -1), sz(n, 1), dup(n, 0), ch(n), val(v), sum(v), lazy(n, id_lazy()){}
  int find(int u){
    while(par[u] != -1) u = par[u];
    return u;
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
  void unite(int u, int v){
    u = find(u), v = find(v);
    if(u == v){
      dup[u]++;
      history.push_back(u + n);
      return;
    }
    if(sz[v] > sz[u]) std::swap(u, v);
    cc--;
    lazy[v] = propagate(lazy[v], inv_lazy(lazy[u]));
    sum[u] = merge(sum[u], sum[v]);
    par[v] = u;
    ch[u].push_back(v);
    sz[u] += sz[v];
    dup[u] += dup[v];
    history.push_back(u);
    history.push_back(v);
  }
  // aを含む連結成分のサイズ
  int size(int a){
    return sz[find(a)];
  }
  // uの値をxに更新
  void update(int u, Val x){
    sum[u] = merge(sum[u], x);
    val[u] = merge(val[u], x);
    while(par[u] != -1){
      sum[par[u]] = merge(sum[par[u]], x);
      u = par[u];
    }
  }
  // uの値
  Val get(int u){
    Val res = val[u];
    Lazy lz = lazy[u];
    while(par[u] != -1){
      lz = propagate(lz, lazy[par[u]]);
      u = par[u];
    }
    return apply(res, lz, 1);
  }
  void set(int u, Val x){
    update(u, merge(x, inv(get(u))));
  }
  // uを含む連結成分の和
  Val query_cc(int u){
    u = find(u);
    return apply(sum[u], lazy[u], 0, sz[u]);
  }
  // uを含む連結成分にxを作用させる
  void update_cc(int u, Lazy x){
    u = find(u);
    lazy[u] = propagate(lazy[u], x);
  }
  // uniteを一回戻す
  bool rollback(){
    assert(!history.empty());
    int c = history.back();
    history.pop_back();
    if(c >= n){
      dup[c - n]--;
      return false;
    }
    int p = history.back();
    history.pop_back();
    cc++;
    par[c] = -1;
    sz[p] -= sz[c];
    dup[p] -= dup[c];
    ch[p].pop_back();
    sum[p] = merge(sum[p], inv(sum[c]));
    lazy[c] = propagate(lazy[c], lazy[p]);
    return true;
  }
  // 戻り値: 頂点uを含む連結成分の全要素
  std::vector<int> enumerate(int u){
    std::vector<int> ret;
    std::queue<int> q;
    q.push(find(u));
    while(!q.empty()){
      int v = q.front();
      q.pop();
      ret.push_back(v);
      for(int c : ch[v]) q.push(c);
    }
    return ret;
  }
  // uを含む連結成分が木か
  int is_tree(int u){
    return dup[find(u)] == 0;
  }
  // uを含む連結成分からあと何本辺を取れば木になるか
  int duplicated(int u){
    return dup[find(u)];
  }
  // 連結成分の数
  int count_cc(){
    return cc;
  }
};
#endif