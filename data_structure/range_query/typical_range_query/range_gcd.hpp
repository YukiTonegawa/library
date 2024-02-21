#ifndef _RANGE_GCD_H_
#define _RANGE_GCD_H_
#include <vector>
#include <cassert>
#include "../../segment_tree/segment_tree.hpp"
#include "../../segment_tree/binary_indexed_tree.hpp"

// Val :　符号付き整数
template<typename Val>
struct range_gcd{
  int n;
  segment_tree<range_gcd<Val>> seg;
  binary_indexed_tree<Val> bit;
  range_gcd(int _n): n(_n), bit(n), seg(n){}
  range_gcd(const std::vector<Val> &v): n(v.size()){
    std::vector<Val> tmp(n, 0);
    if(n){
      tmp[0] = v[0];
      for(int i = 1; i < n; i++) tmp[i] = v[i] - v[i - 1];
    }
    bit = binary_indexed_tree<Val>(tmp);
    seg = segment_tree<range_gcd<Val>>(tmp);
  }
  // 負の数xは常にabs(x)となる
  void set(int i, Val x){
    assert(0 <= i && i < n);
    update(i, i + 1, abs(x) - get(i));
  }
  // 負の数xは常にabs(x)となる
  Val get(int i){
    assert(0 <= i && i < n);
    return abs(bit.query(i + 1));
  }
  // [l, r)にxを足す, 負の数xは常にabs(x)となる
  void update(int l, int r, Val x){
    assert(0 <= l && l < r && r <= n);
    bit.update(l, x);
    if(r < n) bit.update(r, -x);
    seg.update(l, x);
    seg.update(r, -x);
  }
  // [l, r)のgcd, 0は単位元
  Val query(int l, int r){
    assert(0 <= l && l < r && r <= n);
    Val xl = bit.query(l + 1);
    if(l + 1 == r) return xl;
    return range_gcd<Val>::__gcd(xl, seg.query(l + 1, r));
  }
  // gcd[l, r]がx未満になる最小のrとgcd, ない場合は{-1, gcd[l, n)}
  std::pair<int, Val> bisect_left(int l, Val x){
    assert(0 <= l && l < n);
    Val xl = get(l);
    if(xl < x) return {l, x};
    if(l + 1 == n) return {-1, x};
    int r = seg.bisect_from_left(l + 1, [&](Val z){return range_gcd<Val>::__gcd(xl, z) < x;});
    return {r, query(l, r == -1 ? n : r + 1)};
  }
  // gcd[l, r]がx未満になる最大のlとgcd, ない場合は{-1, gcd[0, r]}
  std::pair<int, Val> bisect_right(int r, Val x){
    assert(0 <= r && r < n);
    Val xr = get(r);
    if(xr < x) return {r, x};
    if(r == 0) return {-1, x};
    int l = seg.bisect_from_right(r, [&](Val z){return range_gcd<Val>::__gcd(xr, z) < x;});
    if(l > 0) l--;
    return {l, query(l == -1 ? 0 : l, r + 1)};
  }
  // gcd[l, r]が減少するO(log値)個の{r, gcd[l, r]}を列挙
  std::vector<std::pair<int, Val>> bisect_left_enumerate(int l){
    assert(0 <= l && l < n);
    Val g = get(l);
    std::vector<std::pair<int, Val>> ret{{l, g}};
    while(l < n){
      int r = seg.bisect_from_left(l + 1, [&](Val z){return range_gcd<Val>__gcd(g, z) < g;});
      if(r == -1) return ret;
      l = r;
      g = range_gcd<Val>__gcd(g, seg.get(r));
      ret.push_back({l, g});
    }
    return ret;
  }
  // gcd[l, r]が減少するO(log値)個の{l, gcd[l, r]}を列挙
  std::vector<std::pair<int, Val>> bisect_right_enumerate(int r){
    assert(0 <= r && r < n);
    Val g = get(r);
    std::vector<std::pair<int, Val>> ret{{r, g}};
    while(0 < r){
      int l = seg.bisect_from_right(r, [&](Val z){return range_gcd<Val>__gcd(g, z) < g;});
      if(l > 0) l--;
      if(l == -1) return ret;
      r = l;
      g = range_gcd<Val>__gcd(g, seg.get(l + 1));
      ret.push_back({r, g});
    }
    return ret;
  }
};
#endif