#ifndef _DYNAMIC_BIT_H_
#define _DYNAMIC_BIT_H_

#include "bit_operation.hpp"
struct dynamic_bit64{
  using ull = unsigned long long;
  ull x;
  dynamic_bit64(ull _x): x(_x){}
  // 32bitずつに分割
  std::pair<ull, ull> split_half(){
    return {x >> 32, x & mask_0_32};
  }
  // 先頭k-bitと残りの間にfを挿入
  void insert(int k, bool f){
    ull y = ((ull)1 << k) - 1;
    ull left = x & ~y, right = x & y;
    x = (left << 1) | ((ull)f << k) | right;
  }
  // k番目を削除
  bool erase(int k){
    bool ret = (x >> k) & 1;
    if(k == 63){
      x ^= (ull)ret << k;
    }else{
      ull y = ((ull)1 << (k + 1)) - 1;
      ull left = (x & ~y) >> 1, right = (x & (y >> 1));
      x = left | right;
    }
    return ret;
  }
  // k番目にfをセット
  void set(int k, bool f){
    if(f) x |= ((ull)1 << k);
    else x &= ~((ull)1 << k);
  }
  // k番目の値を取得
  bool get(int k){
    return (x >> k) & 1;
  }
  int popcount(){
    return __builtin_popcountll(x >> 32) + __builtin_popcountll(x & mask_0_32);
  }
  bool access(int k){
    return (x >> k) & 1;
  }
  // rank [0, k) 下位k-bitのpop
  int rank(int k){
    if(k == 64) return popcount();
    return __builtin_popcountll(x & ((1ULL << k) - 1));
  }
  // k番目の1, 無い場合壊れる
  int select1(int k){
    return select_64bit(x, k);
  }
  // k番目の0, 無い場合壊れる
  int select0(int k){
    return select_64bit(~x, k);
  }
};

#endif