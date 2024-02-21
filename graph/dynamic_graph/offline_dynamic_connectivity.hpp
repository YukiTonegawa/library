#ifndef _OFFLINE_DYNAMIC_CONNECTIVITY_H_
#define _OFFLINE_DYNAMIC_CONNECTIVITY_H_
#include <vector>
#include <unordered_map>
#include <cassert>
#include <functional>
#include <algorithm>
#include "../../data_structure/persistent/rollback_union_find.hpp"

// UF := コンストラクタuf(N), 関数unite(a, b), rollback(a, bを持つ
// F: (uf, query_id) -> void
template<typename UF, typename F = std::function<void(UF&, int)>>
struct offline_dynamic_connectivity{
public:
  UF uf;
private:
  static constexpr int inf_time = 1 << 30;
  int n;
  std::unordered_multimap<long long, int> E; // {a, b, inesrt_time}
  std::vector<std::pair<int, F>> Q;
  struct Update{int a, b, link_time, cut_time;};
  std::vector<Update> U{Update{0, 0, 0, 0}};

  void solve(int l, int r, int &j, int &k){
    assert(l != r);
    if(r - l == 1){
      while(k < Q.size() && Q[k].first == l) Q[k++].second(uf, j++);
      return;
    }
    int mid = (l + r) / 2;
    int cnt = 0;
    for(int i = mid; i < r; i++) if(U[i].link_time <= l) cnt++, uf.unite(U[i].a, U[i].b);
    solve(l, mid, j, k);
    while(cnt) cnt--, uf.rollback();
    for(int i = l; i <= mid; i++) if(r <= U[i].cut_time) cnt++, uf.unite(U[i].a, U[i].b);
    solve(mid, r, j, k);
    while(cnt) cnt--, uf.rollback();
  }
public:
  offline_dynamic_connectivity(int n): uf(n), n(n), E(n){}

  // 辺a-bを追加する, すでにある場合は多重辺が追加される
  void link(int a, int b){
    if(a > b) std::swap(a, b);
    if(a == b) return;
    int link_time = U.size();
    E.emplace(((long long)a << 30) + b, link_time);
    U.push_back(Update{a, b, link_time, inf_time});
  }
  // 辺a-bを追加する, すでにある場合は何もしない
  void link_unique(int a, int b){
    if(a > b) std::swap(a, b);
    if(a == b) return;
    auto itr = E.find(((long long)a << 30) + b);
    if(itr != E.end()) return;
    int link_time = U.size();
    E.emplace(((long long)a << 30) + b, link_time);
    U.push_back(Update{a, b, link_time, inf_time});
  }
  // 辺a-bを消す, ない場合は何もしない
  void cut(int a, int b){
    if(a > b) std::swap(a, b);
    if(a == b) return;
    auto itr = E.find(((long long)a << 30) + b);
    if(itr != E.end()){
      int link_time = itr->second;
      int cut_time = U.size();
      U[link_time].cut_time = cut_time;
      U.push_back({a, b, link_time, cut_time});
      E.erase(itr);
    }
  }
  void query(F f){
    int query_time = (int)U.size() - 1;
    Q.push_back(std::make_pair(query_time, f));
  }
  void solve(){
    int j = 0, k = 0;
    solve(0, U.size(), j, k);
  }
};
#endif
