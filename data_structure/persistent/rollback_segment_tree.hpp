#ifndef _ROLLBACK_SEGMENT_TREE_H_
#define _ROLLBACK_SEGMENT_TREE_H_
#include "rollback_array.hpp"
#include "../../algebraic_structure/monoid.hpp"
#include <vector>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <stack>

template<typename monoid>
struct rollback_segment_tree{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  int N, M;
  std::vector<Val> sum;
  rollback_array<Val> hist;
  rollback_segment_tree(){}
  rollback_segment_tree(int n):
  N(n), M(ceil_pow(N, 2)), sum(2 * M - 1, id()){}
  rollback_segment_tree(const std::vector<Val> v):
  N(v.size()), M(ceil_pow(N, 2)), sum(2 * M - 1, id()){
    std::copy(v.begin(), v.end(), sum.begin() + M - 1);
    for(int i = M - 2; i >= 0; i--){
      sum[i] = merge(sum[i * 2 + 1], sum[i * 2 + 2]);
    }
  }
  void set(int k, Val x){
    assert(0 <= k && k < N);
    hist.set(k, x);
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
    assert(0 <= r && r <= N);
    assert(!f(id()));
    if(r == 0) return -1;
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
        return r + 1 - M;
      }
    }while((r & -r) != r);
    return -1;
  }
  void bookmark(){
    hist.bookmark();
  }
  void rollback(){
    for(auto [k, x] : hist.rollback_bookmark_diff()) set(k, x);
  }
};
#endif
