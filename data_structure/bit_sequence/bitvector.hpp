#ifndef _BITVECTOR_H_
#define _BITVECTOR_H_
#include <iostream>
#include <vector>
#include <cassert>
#include <array>

//rank O(1), select O(lgN), memory 2Nbit
struct bitvector_memory{
  static constexpr int s = 32;
  int n;
  std::vector<int> RS{0}, table;
  bitvector_memory(): n(0){}
  bitvector_memory(const std::vector<bool> &S): n(S.size()){
    int pop = 0, m = 0;
    for(int i = 0, t = 0; i < n; i++, t++){
      if(S[i]) pop++, m += 1 << t;
      if(t == s - 1 || i == (n - 1)){
        RS.push_back(pop);
        table.push_back(m);
        t = -1, m = 0;
      }
    }
  }
  // bit[k]
  bool access(int k){
    assert(k < n);
    return (table[k / s] >> (k % s)) & 1;
  }
  // count 1, i < k
  int rank1(int k){
    assert(k <= n);
    int ret = RS[k / s];
    if(k % s) ret += __builtin_popcount(table[k / s] << (s - (k % s)));
    return ret;
  }
  // count 0, i < k
  int rank0(int k){
    return k - rank1(k);
  }
  // kth 1, 0-indexed
  int select1(int k){
    int L = 0, R = n + 1;
    while(R - L > 1){
      int mid = (L + R) >> 1;
      if(rank1(mid) > k) R = mid;
      else L = mid;
    }
    return L;
  }
  // kth 0, 0-indexed
  int select0(int k){
    int L = 0, R = n + 1;
    while(R - L > 1){
      int mid = (L + R) >> 1;
      if(rank0(mid) > k) R = mid;
      else L = mid;
    }
    return L;
  }
  // leftmost1, k < i
  int succ1(int k){
    return select1(rank1(k + 1));
  }
  // leftmost0, k < i
  int succ0(int k){
    return select0(rank0(k + 1));
  }
  // rightmost1, i < k
  int pred1(int k){
    int r = rank1(k);
    return r == 0 ? -1 : select1(r - 1);
  }
  // rightmost0, i < k
  int pred0(int k){
    int r = rank0(k);
    return r == 0 ? -1 : select0(r - 1);
  }
};
// rank O(1), select O(1), memory (2 + 32)Nbit
struct bitvector_fast_select{
  static constexpr int s = 32;
  int n;
  std::vector<int> pos1, pos0;
  std::vector<int> RS{0}, table;
  bitvector_fast_select(){}
  bitvector_fast_select(const std::vector<bool> &v): n(v.size()){
    for(int i = 0; i < n; i++){
      if(v[i]) pos1.push_back(i);
      else pos0.push_back(i);
    }
    int pop = 0, m = 0;
    for(int i = 0, t = 0; i < n; i++, t++){
      if(v[i]) pop++, m += 1 << t;
      if(t == s - 1 || i == (n - 1)){
        RS.push_back(pop);
        table.push_back(m);
        t = -1, m = 0;
      }
    }
  }
  // bit[k]
  bool access(int k){
    assert(k < n);
    return (table[k / s] >> (k % s)) & 1;
  }
  // count 1, i < k
  int rank1(int k){
    assert(k <= n);
    int ret = RS[k / s];
    if(k % s) ret += __builtin_popcount(table[k / s] << (s - (k % s)));
    return ret;
  }
  // count 0, i < k
  int rank0(int k){
    return k - rank1(k);
  }
  int select1(int k){
    if(k >= pos1.size()) return n;
    return pos1[k];
  }
  int select0(int k){
    if(k >= pos0.size()) return n;
    return pos0[k];
  }
  // leftmost1, k < i
  int succ1(int k){
    return select1(rank1(k + 1));
  }
  // leftmost0, k < i
  int succ0(int k){
    return select0(rank0(k + 1));
  }
  // rightmost1, i < k
  int pred1(int k){
    int r = rank1(k);
    return r == 0 ? -1 : select1(r - 1);
  }
  // rightmost0, i < k
  int pred0(int k){
    int r = rank0(k);
    return r == 0 ? -1 : select0(r - 1);
  }
};
/*
https://codeforces.com/contest/1746/my バグ
// rank O(1), select O(1), memory (2 + 4 + (0 ~ 4))N bit
struct bitvector{
private:
  constexpr static int s = 16, th_l = 255, th_m = 64, s2 = 32;
  int n;
  // rank
  std::vector<int> RS{0}, table;

  // select
  struct block{
    int t:2, idx:30;
    block(int t, int idx): t(t), idx(idx){}
  };
  std::array<std::vector<block>, 2> B;
  std::array<std::vector<std::array<int, s>>, 2> large;
  std::array<std::vector<std::array<uint8_t, s>>, 2> medium;
  std::array<std::vector<int>, 2> block_idx;

  static uint64_t flip(uint64_t x){
    constexpr static uint64_t y = 0xffffffffffffffffULL;
    return x ^ y;
  }
  static int find_kth_set_bit(uint64_t mask, int k){
    constexpr static uint64_t m1 = 0x5555555555555555ULL; // even bits
    constexpr static uint64_t m2 = 0x3333333333333333ULL; // even 2-bit groups
    constexpr static uint64_t m4 = 0x0f0f0f0f0f0f0f0fULL; // even nibbles
    constexpr static uint64_t m8 = 0x00ff00ff00ff00ffULL; // even bytes
    int t, i = k, r = 0;
    uint64_t c1 = mask;
    uint64_t c2 = c1 - ((c1 >> 1) & m1);
    uint64_t c4 = ((c2 >> 2) & m2) + (c2 & m2);
    uint64_t c8 = ((c4 >> 4) + c4) & m4;
    uint64_t c16 = ((c8 >> 8) + c8) & m8;
    uint64_t c32 = (c16 >> 16) + c16;
    int c64 = (int)(((c32 >> 32) + c32) & 0x7f);
    t = (c32    ) & 0x3f; if (i >= t) { r += 32; i -= t; }
    t = (c16>> r) & 0x1f; if (i >= t) { r += 16; i -= t; }
    t = (c8 >> r) & 0x0f; if (i >= t) { r +=  8; i -= t; }
    t = (c4 >> r) & 0x07; if (i >= t) { r +=  4; i -= t; }
    t = (c2 >> r) & 0x03; if (i >= t) { r +=  2; i -= t; }
    t = (c1 >> r) & 0x01; if (i >= t) { r +=  1;         }
    if (k >= c64) r = -1;
    return r;
  }
  // get[l, r], len <= 64
  uint64_t get64(int l, int r){
    int lb = l / s2, rb = r / s2;
    l %= s2, r %= s2;
    if(rb - lb == 2){
      uint64_t ans_l = table[lb] >> l;
      uint32_t ans_r = table[rb] << (s2 - 1 - r);
      return ans_l + (table[lb + 1] << (s2 - l)) + ((uint64_t)ans_r << (r - l + s2));
    }else if(rb - lb == 1){
      uint64_t ans_l = table[lb] >> l;
      uint32_t ans_r = table[rb] << (s2 - 1 - r);
      return ans_l + ((uint64_t)ans_r << (r - l));
    }else{
      assert(lb == rb);
      uint32_t ans = table[lb];
      ans <<= (s2 - 1 - r);
      ans >>= (s2 - (r - l + 1));
      return ans;
    }
  }
  int select(int k, bool b){
    int idx_l = k / s;
    if(idx_l >= block_idx[b].size()) return n;
    int ans = block_idx[b][idx_l], block_type = B[b][idx_l].t, idx_s = B[b][idx_l].idx;
    if(block_type == 2) ans += large[b][idx_s][k % s];
    else if(block_type == 1) ans += medium[b][idx_s][k % s];
    else{
      int L = block_idx[b][idx_l], R = B[b][idx_l].idx;
      uint64_t x = b ? get64(L, R) : flip(get64(L, R));
      int tmp = find_kth_set_bit(x, k % s);
      ans = (tmp == -1 || ans + tmp >= n ? n : ans + tmp);
    }
    return ans;
  }
public:
  bitvector(){}
  bitvector(const std::vector<bool> &v): n(v.size()){
    // rank
    int pop = 0, m = 0;
    for(int i = 0, t = 0; i < n; i++, t++){
      if(v[i]) pop++, m += 1 << t;
      if(t == s2 - 1 || i == (n - 1)){
        RS.push_back(pop);
        table.push_back(m);
        t = -1, m = 0;
      }
    }
    // select
    // 区間が[l, r]の0/1のブロックを作成
    auto make_block = [&](int l, int r, bool b){
      block_idx[b].push_back(l);
      if(r - l >= th_l){
        B[b].push_back(block(2, large[b].size()));
        large[b].push_back(std::array<int, s>());
        int k = 0;
        for(int i = l; i <= r; i++) if(v[i] == b) large[b].back()[k++] = i - l;
        while(k < s) large[b].back()[k++] = r - l + 1;
      }else if(r - l >= th_m){
        B[b].push_back(block(1, medium[b].size()));
        medium[b].push_back(std::array<uint8_t, s>());
        int k = 0;
        for(int i = l; i <= r; i++) if(v[i] == b) medium[b].back()[k++] = i - l;
        while(k < s) large[b].back()[k++] = r - l + 1;// [0, 255]
      }else{
        B[b].push_back(block(0, r));
      }
    };
    for(int b = 0; b < 2; b++){
      for(int i = 0, k = 0, l = 0; i < n; i++){
        if(v[i] == b){
          k++;
          if(k % s == 0) make_block(l, i, b), l = i + 1;
        }
        if(i == n - 1 && (block_idx[b].empty() || block_idx[b].back() != l)) make_block(l, i, b);
      }
    }
  }
  int size(){
    return n;
  }
  // bit[k]
  bool access(int k){
    assert(k < n);
    return (table[k / s2] >> (k % s2)) & 1;
  }
  // count 1, i < k
  int rank1(int k){
    assert(0 <= k && k <= n);
    int ret = RS[k / s2];
    if(k % s2) ret += __builtin_popcount(table[k / s2] << (s2 - (k % s2)));
    return ret;
  }
  // count 0, i < k
  int rank0(int k){
    return k - rank1(k);
  }
  int select1(int k){
    return select(k, 1);
  }
  int select0(int k){
    return select(k, 0);
  }
};
*/

template<int R>
struct bitvector_arbitrary_radix{
  static constexpr int s = 32, sdiv = 5, smod = 31;
  using Z = uint32_t;//s bit
  int N;
  std::vector<std::array<int, R>> B;
  std::vector<std::array<Z, R>> S;

  bitvector_arbitrary_radix(): N(0){}
  // init O(N + NR/s)
  bitvector_arbitrary_radix(const std::vector<uint8_t> &v): N(v.size()){
    int M = (N + s - 1) / s;
    std::array<int, R> pop;
    std::array<Z, R> pop_small;
    pop.fill(0);
    pop_small.fill(0);
    B.resize(M + 1, pop);
    S.resize(M, pop_small);
    for(int i = 0, t = 0, sz = 0; i < N; i++, t++){
      int x = v[i];
      assert(0 <= x && x < R);
      pop[x]++, pop_small[x] |= (Z(1) << t);
      if(t == s - 1 || i == N - 1){
        for(int j = 0; j < R; j++){
          if(j) pop[j] += pop[j - 1], pop_small[j] |= pop_small[j - 1];
          B[sz + 1][j] = pop[j] + B[sz][j];
          S[sz][j] = pop_small[j];
        }
        pop.fill(0);
        pop_small.fill(0);
        t = -1;
        sz++;
      }
    }
  }
  // r未満のcの数
  int rank(int r, int c){
    int rq = r >> sdiv, rm = r & smod;
    int ret = B[rq][c] - (c ? B[rq][c - 1] : 0);
    if(rm) ret += __builtin_popcount((S[rq][c] ^ (c ? S[rq][c - 1] : 0)) << (s - rm));
    return ret;
  }
  // r未満のc未満の数
  int rank_lower(int r, int c){
    if(c == 0) return 0;
    int rq = r >> sdiv, rm = r & smod;
    int ret = B[rq][--c];
    if(rm) ret += __builtin_popcount(S[rq][c] << (s - rm));
    return ret;
  }
};
#endif