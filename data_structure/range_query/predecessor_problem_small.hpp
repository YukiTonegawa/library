#ifndef _PREDECESSOR_PROBLEM_SMALL_H_
#define _PREDECESSOR_PROBLEM_SMALL_H_
#include <array>
#include <vector>
#include "../bit_sequence/bit_operation.hpp"
// [0, max_elem)
template<int max_elem>
struct predecessor_problem_small{
  using ull = unsigned long long;
  static constexpr ull inf = ~(ull)0;
  static constexpr int bitlen = 64;
  static constexpr int bitlen_mod = 63;
  static constexpr int bitlen_div = 6;
  static constexpr int block = (max_elem + bitlen - 1) / bitlen;
  std::array<ull, block> v;
  predecessor_problem_small(bool f = 0){v.fill(f ? inf : (ull)0);}
  // O(1)
  void set(int k, bool f){
    bool g = (v[k >> bitlen_div] >> (k & bitlen_mod)) & 1;
    if(f != g) v[k >> bitlen_div] ^= (ull)1 << (k & bitlen_mod);
  }
  // O(1)
  bool get(int k){
    return (v[k >> bitlen_div] >> (k & bitlen_mod)) & 1;
  }
  // O(max_elem / bitlen), ない場合は-1
  int find_next1(int k){
    int a = k >> bitlen_div, b = k & bitlen_mod;
    int res = find_next_64bit(v[a], b);
    if(res != -1) return (a << bitlen_div) + res < max_elem ? (a << bitlen_div) + res : -1;
    while(++a < block){
      if(v[a]){
        res = (a << bitlen_div) + __builtin_ctzll(v[a]);
        return res < max_elem ? res : -1;
      }
    }
    return -1;
  }
  // O(max_elem / bitlen), ない場合は-1
  int find_prev1(int k){
    int a = k >> bitlen_div, b = k & bitlen_mod;
    int res = find_prev_64bit(v[a], b);
    if(res != -1) return (a << bitlen_div) + res < max_elem ? (a << bitlen_div) + res : -1;
    while(--a >= 0){
      if(v[a]){
        return (a << bitlen_div) + 63 - __builtin_clzll(v[a]);
      }
    }
    return -1;
  }
};
#endif
