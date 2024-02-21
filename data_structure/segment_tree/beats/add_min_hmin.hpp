#ifndef _BEAATS_ADD_MIN_HMIN_H_
#define _BEAATS_ADD_MIN_HMIN_H_
#include <limits>
#include <algorithm>
#include "../lazy_segment_tree.hpp"

template<typename T>
struct monoid_add_min_hmin{
  using Val = std::pair<T, T>;
  using Lazy = Val;
  static constexpr T inf = std::numeric_limits<T>::max() / 2;
  // {hmin, min}
  static Val id(){return {inf, inf};}
  // {pre, add}
  static Lazy id_lazy(){return {inf, 0};}
  static Val merge(Val a, Val b){return {std::min(a.first, b.first), std::min(a.second, b.second)};}
  static Val apply(Val a, Lazy b, int l, int r){return {std::min(a.first, a.second + b.first), a.second + b.second};}
  static Lazy propagate(Lazy a, Lazy b){return {std::min(a.first, a.second + b.first), a.second + b.second};}
};

// 区間加算-区間historic_minのmin
// O(logN)
template<typename Val>
struct beats_add_min_hmin{
  lazy_segment_tree<monoid_add_min_hmin<Val>> seg;
  beats_add_min_hmin(const std::vector<Val> &v){
    int n = v.size();
    std::vector<std::pair<Val, Val>> tmp(n);
    for(int i = 0; i < n; i++) tmp[i] = {v[i], v[i]};
    seg = lazy_segment_tree<monoid_add_min_hmin<Val>>(tmp);
  }
  void update_add(int l, int r, Val x){
    seg.update(l, r, {x, x});
  }
  Val get(int k){
    return seg.get(k).second;
  }
  Val query_min(int l, int r){
    return seg.query(l, r).second;
  }
  Val get_hmin(int k){
    return seg.get(k).first;
  }
  Val query_min_hmin(int l, int r){
    return seg.query(l, r).first;
  }
};
#endif