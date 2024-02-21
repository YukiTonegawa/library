#ifndef _NAZOKI_H_
#define _NAZOKI_H_
#include <vector>
#include <cassert>
#include <array>
#include <algorithm>
#include <string>
// [0, 64^D)
template<int D>
struct nazotree{
  static_assert(1 <= D && D <= 4); // 64, 4096, 262144, 16777216
  using ull = unsigned long long;
  static constexpr int bitlen = 64;
  static constexpr int bitlen_mod = 63;
  static constexpr int bitlen_div = 6;
  static constexpr int node_count = (1 << (6 * D)) / (bitlen - 1); // 1 + 64 + 64^2 ... 64^(D-1)
  std::array<ull, node_count> flag;
  nazotree(){flag.fill(0);}
  nazotree(const std::string &s, char one = '1'){
    flag.fill(0);
    int n = std::min((int)s.size(), 1 << (6 * D));
    for(int i = 0; i < n; i++){
      if(s[i] == one){
        int j = (i + node_count - 1);
        flag[j >> bitlen_div] |= 1ULL << (j & bitlen_mod);
      }
    }
    for(int i = node_count - 1; i > 0; i--){
      if(flag[i]) flag[(i - 1) >> bitlen_div] |= 1ULL << ((i - 1) & bitlen_mod);
    }
  }
  // kを追加. すでにある場合は何もしない
  void insert(int k){
    k += node_count;
    while(k--){
      int dir = k & bitlen_mod;
      k >>= bitlen_div;
      flag[k] |= 1ULL << dir;
    }
  }
  // kを削除. すでにない場合は何もしない
  void erase(int k){
    k += node_count;
    bool f = true;
    while(k--){
      int dir = k & bitlen_mod;
      k >>= bitlen_div;
      flag[k] &= ~((ull)f << dir);
      f = !flag[k];
    }
  }
  bool find(int k){
    k += node_count;
    int p0 = (k - 1) >> bitlen_div, d0 = (k - 1) & bitlen_mod;
    return (flag[p0] >> d0) & 1;
  }
  // k以上の最小要素, 存在しない場合-1
  int successor(int k){
    if(find(k)) return k;
    k += node_count;
    while(k--){
      int dir = k & bitlen_mod;
      k >>= bitlen_div;
      ull f = (flag[k] >> dir) >> 1; // dirより大きい1が存在する
      if(f){
        k = (k << bitlen_div) + dir + 2 + __builtin_ctzll(f);
        while(k < node_count) k = __builtin_ctzll(flag[k]) + 1 + (k << bitlen_div);
        return k - node_count;
      }
    }
    return -1;
  }
  // k以下の最大要素, 存在しない場合-1
  int predecessor(int k){
    if(find(k)) return k;
    k += node_count;
    while(k--){
      int dir = k & bitlen_mod;
      k >>= bitlen_div;
      ull f = flag[k] & ((1ULL << dir) - 1); // dirより小さい1が存在する
      if(f){
        k = (k << bitlen_div) + 64 - __builtin_clzll(f);
        while(k < node_count) k = 64 - __builtin_clzll(flag[k]) + (k << bitlen_div);
        return k - node_count;
      }
    }
    return -1;
  }
};

#endif