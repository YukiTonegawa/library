#ifndef _SYNCHRONOUS_H_
#define _SYNCHRONOUS_H_
#include <limits>
#include <algorithm>
#include <cassert>
#include "../segment_tree/dual_segment_tree.hpp"
#include "../segment_tree/binary_indexed_tree_range_add.hpp"

template<typename Val>
struct syncronous_chmax{
  using Cl = clamp_function_score<Val, Val>;
  static constexpr Val inf = Cl::inf;
  static constexpr Val minf = Cl::minf;
private:
  int n;
  dual_segment_tree_monoid<range_clamp_score<Val, Val>> seg;
  binary_indexed_tree_range_add<Val> bit;
public:
  syncronous_chmax(): n(0){}
  // n * 0
  syncronous_chmax(int n): n(n), seg(n), bit(n){
    std::vector<Cl> tmp(n);
    for(int i = 0; i < n; i++) tmp[i] = Cl(0, 0, 0);
    seg = dual_segment_tree_monoid<range_clamp_score<Val, Val>>(tmp);
  }
  syncronous_chmax(std::vector<Val> A, const std::vector<Val> &B): n(A.size()), bit(B){
    assert(B.size() == n);
    std::vector<Cl> tmp(n);
    for(int i = 0; i < n; i++){
      A[i] -= B[i];
      tmp[i] = Cl(0, A[i], A[i]);
    }
    seg = dual_segment_tree_monoid<range_clamp_score<Val, Val>>(tmp);
  }
  // ak <- x
  void set_A(int k, Val x){
    Val y = get_B(k);
    seg.set(k, {x - y, minf, inf});
    Val f = bit.query(k, k + 1);
    bit.update(k, y - f);
  }
  // bk <- y
  void set_B(int k, Val y){
    Val x = get_A(k);
    seg.set(k, {x - y, minf, inf});
    Val f = bit.query(k, k + 1);
    bit.update(k, y - f);
  }
  // ak
  Val get_A(int k){
    auto x = seg.get(k);
    return bit.query(k, k + 1) + x.lower + x.score_sum;
  }
  // bk
  Val get_B(int k){
    return bit.query(k, k + 1) + seg.get(k).score_sum;
  }
  // [l, r)のai <- max(ai + x, bi + y)
  void update_chmax_A(int l, int r, Val x, Val y){
    seg.update(l, r, {x, y, inf});
  }
  // [l, r)のbi <- max(ai + x, bi + y)
  void update_chmax_B(int l, int r, Val x, Val y){
    seg.update(l, r, {-y, minf, -x});
    bit.update(l, r, y);
  }
  // A[l, r)にx, B[l, r)にyを足す
  void update_add(int l, int r, Val x, Val y){
    seg.update(l, r, {x - y, minf, inf});
    bit.update(l, r, y);
  }
  // A[0, n)にx, B[0, n)にyを足す
  void update_add_all(Val x, Val y){
    seg.update_all({x - y, minf, inf});
    bit.update(0, n, y);
  }
};

#endif
