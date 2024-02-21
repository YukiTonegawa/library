#ifndef _SEGMENT_TREE_H_
#define _SEGMENT_TREE_H_
#include <cassert>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "../../algebraic_structure/monoid.hpp"
template<typename monoid>
struct segment_tree{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  int N, M;
  std::vector<Val> sum;
  int ceil_pow2(int y){
    int x = 0;
    while ((1U << x) < (unsigned int)(y)) x++;
    return x;
  };
  segment_tree(){}
  segment_tree(int n): N(n), M(1 << ceil_pow2(N)), sum(2 * M - 1, id()){}
  segment_tree(const std::vector<Val> &v): N(v.size()), M(1 << ceil_pow2(N)), sum(2 * M - 1, id()){
    std::copy(v.begin(), v.end(), sum.begin() + M - 1);
    for(int i = M - 2; i >= 0; i--){
      sum[i] = merge(sum[i * 2 + 1], sum[i * 2 + 2]);
    }
  }
  int size(){return N;}
  void set(int k, Val x){
    assert(0 <= k && k < N);
    k += M - 1;
    sum[k] = x;
    while(k){
      k = (k - 1) >> 1;
      sum[k] = merge(sum[k * 2 + 1], sum[k * 2 + 2]);
    }
  }
  Val get(int k){
    assert(0 <= k && k < N);
    return sum[M - 1 + k];
  }
  Val query(int l, int r){
    l = std::max(l, 0), r = std::min(r, N);
    assert(l <= r);
    l += M, r += M;
    Val L = id(), R = L;
    while(l < r){
      if(l & 1) L = merge(L, sum[(l++) - 1]);
      if(r & 1) R = merge(sum[(--r) - 1], R);
      l >>= 1;
      r >>= 1;
    }
    return merge(L, R);
  }
  Val query_all(){
    return sum[0];
  }
  // f(sum[l, r])がtrueになる最左のr. ない場合は-1
  template<typename F>
  int bisect_from_left(int l, const F &f){
    assert(0 <= l);
    assert(!f(id()));
    if(l >= N) return -1;
    l += M;
    Val ret = id();
    do{
      while(l % 2 == 0) l >>= 1;
      if(f(merge(ret, sum[l - 1]))){
        while(l < M){
          l = 2 * l;
          if(!f(merge(ret, sum[l - 1]))){
            ret = merge(ret, sum[l - 1]);
            l++;
          }
        }
        return l - M;
      }
      ret = merge(ret, sum[l - 1]);
      l++;
    }while((l & -l) != l);
    return -1;
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合は-1
  template<typename F>
  int bisect_from_right(int r, const F &f){
    assert(0 <= r && r < N);
    assert(!f(id()));
    r++;
    r += M;
    Val ret = id();
    do{
      r--;
      while(r > 1 && (r % 2)) r >>= 1;
      if(f(merge(sum[r - 1], ret))){
        while(r < M){
          r = 2 * r + 1;
          if(!f(merge(sum[r - 1], ret))){
            ret = merge(sum[r - 1], ret);
            r--;
          }
        }
        return r - M;
      }
      ret = merge(sum[r - 1], ret);
    }while((r & -r) != r);
    return -1;
  }
  std::vector<Val> to_list(){
    return std::vector<Val>(sum.begin() + M - 1, sum.begin() + M - 1 + N);
  }
};
#endif
