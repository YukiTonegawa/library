#ifndef _BINARY_INDEXED_TREE_SET_H_
#define _BINARY_INDEXED_TREE_SET_H_
#include "../data_structure/segment_tree/binary_indexed_tree.hpp"
#include <vector>

struct binary_indexed_tree_set{
  using ull = unsigned long long;
  static constexpr int bitlen = 64;
  static constexpr int bitlen_mod = 63;
  static constexpr int bitlen_shift = 6;
  int n, m;
  std::vector<ull> val;
  binary_indexed_tree<int> bit;
  binary_indexed_tree_set(): n(0), m((n + bitlen_mod) >> bitlen_shift), val(m, 0), bit(m){}
  binary_indexed_tree_set(int n): n(n), m((n + bitlen_mod) >> bitlen_shift), val(m, 0), bit(m){}
  // xを追加(trueを返す), すでにある場合は何もしない(falseを返す)
  bool insert(int x){
    assert(0 <= x && x < n);
    int i = x >> bitlen_shift;
    int j = x & bitlen_mod;
    if((val[i] >> j) & 1) return false;
    val[i] ^= ull(1) << j;
    bit.update(i, 1);
    return true;
  }
  // xを削除(trueを返す), すでにない場合は何もしない(falseを返す)
  bool erase(int x){
    assert(0 <= x && x < n);
    int i = x >> bitlen_shift;
    int j = x & bitlen_mod;
    if(!((val[i] >> j) & 1)) return false;
    val[i] ^= ull(1) << j;
    bit.update(i, -1);
    return true;
  }
  // xがあるか
  bool find(int x){
    assert(0 <= x && x < n);
    int i = x >> bitlen_shift;
    int j = x & bitlen_mod;
    return (val[i] >> j) & 1;
  }
  // [0, r)の要素数
  int rank1(int r){
    assert(0 <= r && r <= n);
    int i = r >> bitlen_shift;
    int j = r & bitlen_mod;
    return bit.query(i) + (!j ? 0 : __builtin_popcountll(val[i] << (bitlen - j)));
  }
  // [0, r)の存在しない要素の数
  int rank0(int r){
    return r - rank1(r);
  }
  // select, findprev, findnext
};

#endif