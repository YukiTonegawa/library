#ifndef _BIT_OPERATION_H_
#define _BIT_OPERATION_H_
#include <stdint.h>
#include <vector>
#include <array>
static constexpr __uint128_t mask_0_64 = ((__uint128_t)1 << 64) - 1;
//static constexpr uint64_t mask_32_64 = 0xFFFFFFFF00000000;
static constexpr uint64_t mask_0_32  = 0x00000000FFFFFFFF;
static constexpr uint64_t mask_48_64 = 0xFFFF000000000000;
static constexpr uint64_t mask_32_48 = 0x0000FFFF00000000;
static constexpr uint64_t mask_16_32 = 0x00000000FFFF0000;
static constexpr uint64_t mask_0_16  = 0x000000000000FFFF;
static constexpr int TABLE_SIZE_LOG = 16, TABLE_SIZE = 1 << TABLE_SIZE_LOG;

using __table =  std::vector<std::array<int8_t, TABLE_SIZE_LOG>>;
using __table_p =  std::vector<std::array<std::pair<int8_t, int8_t>, TABLE_SIZE_LOG>>;
__table select_build(){
  __table res(TABLE_SIZE);
  for(int i = 0; i < TABLE_SIZE; i++){
    res[i].fill(-1);
    int pcnt = 0;
    for(int j = 0; j < TABLE_SIZE_LOG; j++) if((i >> j) & 1) res[i][pcnt++] = j;
  }
  return res;
}
// k-bit目以降(kも含む)に初めて現れる1の位置, 無い場合は-1
int find_next_32bit(uint32_t x, int k){
  uint32_t b = x >> k;
  if(!b) return -1;
  return k + __builtin_ctz(b);
}
// k-bit目以降(kも含む)に初めて現れる1の位置, 無い場合は-1
int find_next_64bit(uint64_t x, int k){
  uint64_t b = x >> k;
  if(!b) return -1;
  return k + __builtin_ctzll(b);
}
// 0 <= k <= 63
// k-bit目以前(kも含む)に初めて現れる1の位置, 無い場合は-1
int find_prev_64bit(uint64_t x, int k){
  uint64_t b = x << (63 - k);
  if(!b) return -1;
  return k - __builtin_clzll(b);
}
// k番目(0-indexed)の1の場所(0-indexed)を返す. 無い場合壊れる
int select_32bit(uint32_t x, int k){
  static __table table = select_build();
  int r = __builtin_popcount(x & mask_0_16);
  if(r > k) return table[x & mask_0_16][k];
  return 16 + table[(x & mask_16_32) >> 16][k - r];
}
// k番目(0-indexed)の1の場所(0-indexed)を返す. 無い場合壊れる
int select_64bit(uint64_t x, int k){
  static __table table = select_build();
  int r = __builtin_popcount(x & mask_0_32);
  if(r > k){
    int rr = __builtin_popcount(x & mask_0_16);
    if(rr > k) return table[x & mask_0_16][k];
    else return 16 + table[(x & mask_16_32) >> 16][k - rr];
  }else{
    k -= r;
    int lr = __builtin_popcountll(x & mask_32_48);
    if(lr > k) return 32 + table[(x & mask_32_48) >> 32][k];
    else return 48 + table[(x & mask_48_64) >> 48][k - lr];
  }
}
// 先頭k_bit(0 <= k <= 32)の1の数
int rank_32bit(uint32_t x, int k){
  return k == 32 ? __builtin_popcount(x) : __builtin_popcount(x & ((1ULL << k) - 1));
}
// 先頭k_bit(0 <= k <= 64)の1の数
int rank_64bit(uint64_t x, int k){
  return k == 64 ? __builtin_popcountll(x) : __builtin_popcountll(x & ((1ULL << k) - 1));
}
// 128bit
int pop_count_128bit(__uint128_t x){
  return __builtin_popcountll(x >> 64) + __builtin_popcountll(x & mask_0_64);
}
int rank_128bit(__uint128_t x, int k){
  if(k == 128) return pop_count_128bit(x);
  if(k < 64) return __builtin_popcountll((x & mask_0_64) & ((1ULL << k) - 1));
  k -= 64;
  return __builtin_popcountll(x & mask_0_64) + __builtin_popcountll((x >> 64)  & ((1ULL << k) - 1));
}
// k番目の1, 無い場合壊れる
int select1_128bit(__uint128_t x, int k){
  int left_pop = __builtin_popcountll(x & mask_0_64);
  if(left_pop > k) return select_64bit(x & mask_0_64, k);
  return 64 + select_64bit(x >> 64, k - left_pop);
}
// k番目の0, 無い場合壊れる
int select0_128bit(__uint128_t x, int k){
  __uint128_t y = ~x;
  int left_unpop = __builtin_popcountll(y & mask_0_64);
  if(left_unpop > k) return select_64bit(y & mask_0_64, k);
  return 64 + select_64bit(y >> 64, k - left_unpop);
}
// k-bit目以降(kも含む)に初めて現れる1の位置, 無い場合は-1
int find_next_128bit(__uint128_t x, int k){
  __uint128_t b = x >> k;
  if(!b) return -1;
  // 末尾64bitが0
  if(!(b & mask_0_64)) return k + 64 + __builtin_ctzll(b >> 64);
  return k + __builtin_ctzll(b);
}
// 0 <= k <= 63
// k-bit目以前(kも含む)に初めて現れる1の位置, 無い場合は-1
int find_prev_128bit(__uint128_t x, int k){
  __uint128_t b = x << (127 - k);
  if(!b) return -1;
  // 先頭64bitが0
  if(!(b >> 64)) return k - 64 - __builtin_clzll(b);
  return k - __builtin_clzll(b >> 64);
}
#endif
