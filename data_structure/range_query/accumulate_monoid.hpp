#ifndef _ACCUMULATE_MONOID_H_
#define _ACCUMULATE_MONOID_H_
#include <vector>
#include <algorithm>
template<typename monoid>
struct accumulate_monoid{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  std::vector<Val> sum;
  accumulate_monoid(){}
  accumulate_monoid(const std::vector<Val> &v): sum(v){
    for(int i = 1; i < v.size(); i++) sum[i] = merge(sum[i - 1], v[i]);
  }
  // [0, r)のモノイド積, 範囲外の部分は全て単位元
  Val query(int r){
    r = std::min(r, (int)sum.size());
    if(r <= 0) return id();
    return sum[r - 1];
  }
  void push_back(Val x){
    Val y = (sum.empty() ? id() : sum.back());
    sum.push_back(merge(y, x));
  }
  void pop_back(){
    assert(!sum.empty());
    sum.pop_back();
  }
};
template<typename monoid>
struct accumulate_monoid_de{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  std::vector<Val> sum_l, sum_r;
  accumulate_monoid_de(){}
  accumulate_monoid_de(const std::vector<Val> &v): sum_l(v), sum_r(v){
    for(int i = 1; i < v.size(); i++) sum_l[i] = merge(sum_l[i - 1], v[i]);
    for(int i = (int)v.size() - 2; i >= 0; i--) sum_r[i] = merge(v[i], sum_r[i + 1]);
  }
  // [0, r)のモノイド積, 範囲外の部分は全て単位元
  Val query_front(int r){
    r = std::min(r, (int)sum_l.size());
    if(r <= 0) return id();
    return sum_l[r - 1];
  }
  // [l, n)のモノイド積, 範囲外の部分は全て単位元
  Val query_back(int l){
    l = std::max(0, l);
    if(l >= sum_r.size()) return id();
    return sum_r[l];
  }
};
#endif
