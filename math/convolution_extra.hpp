#ifndef _CONVOLUTION_EXTRA_H_
#define _CONVOLUTION_EXTRA_H_
#include "convolution.hpp"
// @param 答えが__uint128_tに収まる
std::vector<__uint128_t> convolution_128bit(const std::vector<unsigned long long> &a, const std::vector<unsigned long long> &b){
  int n = int(a.size()), m = int(b.size());
  if(!n || !m) return {};
  if(std::min(n, m) <= 60){
    std::vector<__uint128_t> ans(n + m - 1, 0);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        ans[i + j] += (__uint128_t)a[i] * b[j];
      }
    }
    return ans;
  }
  static constexpr unsigned long long MOD1 = 754974721;  // 2^24
  static constexpr unsigned long long MOD2 = 998244353;  // 2^23
  static constexpr unsigned long long MOD3 = 167772161;  // 2^25
  static constexpr unsigned long long MOD4 = 469762049;  // 2^26
  static constexpr unsigned long long MOD5 = 1224736769; // 2^24
  static constexpr unsigned long long M1M2 = MOD1 * MOD2;
  static constexpr __uint128_t M1M2M3 = (__uint128_t)M1M2 * MOD3;
  static constexpr __uint128_t M1M2M3M4 = M1M2M3 * MOD4;
  static constexpr unsigned long long IM1_M2 = inv_gcd(MOD1, MOD2).second;
  static constexpr unsigned long long IM1M2_M3 = inv_gcd(MOD1 * MOD2, MOD3).second;
  static constexpr unsigned long long IM1M2M3_M4 = inv_gcd((long long)(M1M2M3 % MOD4), MOD4).second;
  static constexpr unsigned long long IM1M2M3M4_M5 = inv_gcd((long long)(M1M2M3M4 % MOD5), MOD5).second;
  auto c1 = convolution_mod<MOD1>(a, b);
  auto c2 = convolution_mod<MOD2>(a, b);
  auto c3 = convolution_mod<MOD3>(a, b);
  auto c4 = convolution_mod<MOD4>(a, b);
  auto c5 = convolution_mod<MOD5>(a, b);
  std::vector<__uint128_t> c(n + m - 1);
  for(int i = 0; i < n + m - 1; i++){
    __uint128_t x = c1[i], y, t;
    t = ((c2[i] < x ? MOD2 + c2[i] - x : c2[i] - x) * IM1_M2) % MOD2;
    x += t * MOD1; // x < 2^60
    y = x % MOD3;
    t = ((c3[i] < y ? MOD3 + c3[i] - y : c3[i] - y) * IM1M2_M3) % MOD3;
    x += t * M1M2; // x < 2^90
    y = x % MOD4;
    t = ((c4[i] < y ? MOD4 + c4[i] - y : c4[i] - y) * IM1M2M3_M4) % MOD4;
    x += t * M1M2M3; // x < 2^120
    y = x % MOD5;
    t = ((c5[i] < y ? MOD5 + c5[i] - y : c5[i] - y) * IM1M2M3M4_M5) % MOD5;
    // 上の解の集合を x_5 + k * M1M2M3M4M5, 0 <= x_5 < M1M2M3M4M5とすると
    // 求めたい解が2^128未満にただ1つ存在することがわかっており, M1M2M3M4M5 > 2^128であるためk = 0
    /// t * M1M2M3M4 <= M1M2M3M4M5 - M1M2M3M4
    //             x_4 < M1M2M3M4
    // であるため, x_4 + t * M1M2M3M4 = x_5 < 2^128でオーバーフローしない
    x += t * M1M2M3M4;
    c[i] = x;
  }
  return c;
}
// @param 答えが__uint128_tに収まる
std::vector<__uint128_t> square_128bit(const std::vector<unsigned long long> &a){
  int n = int(a.size());
  if(!n) return {};
  if(n <= 60){
    std::vector<__uint128_t> ans(2 * n - 1, 0);
    for(int i = 0; i < n; i++){
      for(int j = i + 1; j < n; j++){
        ans[i + j] += __uint128_t(2) * a[i] * a[j];
      }
      ans[2 * i] += __uint128_t(a[i]) * a[i];
    }
    return ans;
  }
  static constexpr unsigned long long MOD1 = 754974721;  // 2^24
  static constexpr unsigned long long MOD2 = 998244353;  // 2^23
  static constexpr unsigned long long MOD3 = 167772161;  // 2^25
  static constexpr unsigned long long MOD4 = 469762049;  // 2^26
  static constexpr unsigned long long MOD5 = 1224736769; // 2^24
  static constexpr unsigned long long M1M2 = MOD1 * MOD2;
  static constexpr __uint128_t M1M2M3 = (__uint128_t)M1M2 * MOD3;
  static constexpr __uint128_t M1M2M3M4 = M1M2M3 * MOD4;
  static constexpr unsigned long long IM1_M2 = inv_gcd(MOD1, MOD2).second;
  static constexpr unsigned long long IM1M2_M3 = inv_gcd(MOD1 * MOD2, MOD3).second;
  static constexpr unsigned long long IM1M2M3_M4 = inv_gcd((long long)(M1M2M3 % MOD4), MOD4).second;
  static constexpr unsigned long long IM1M2M3M4_M5 = inv_gcd((long long)(M1M2M3M4 % MOD5), MOD5).second;
  auto c1 = square_mod<MOD1>(a);
  auto c2 = square_mod<MOD2>(a);
  auto c3 = square_mod<MOD3>(a);
  auto c4 = square_mod<MOD4>(a);
  auto c5 = square_mod<MOD5>(a);
  std::vector<__uint128_t> c(2 * n - 1);
  for(int i = 0; i < 2 * n - 1; i++){
    __uint128_t x = c1[i], y, t;
    t = ((c2[i] < x ? MOD2 + c2[i] - x : c2[i] - x) * IM1_M2) % MOD2;
    x += t * MOD1; // x < 2^60
    y = x % MOD3;
    t = ((c3[i] < y ? MOD3 + c3[i] - y : c3[i] - y) * IM1M2_M3) % MOD3;
    x += t * M1M2; // x < 2^90
    y = x % MOD4;
    t = ((c4[i] < y ? MOD4 + c4[i] - y : c4[i] - y) * IM1M2M3_M4) % MOD4;
    x += t * M1M2M3; // x < 2^120
    y = x % MOD5;
    t = ((c5[i] < y ? MOD5 + c5[i] - y : c5[i] - y) * IM1M2M3M4_M5) % MOD5;
    // 上の解の集合を x_5 + k * M1M2M3M4M5, 0 <= x_5 < M1M2M3M4M5とすると
    // 求めたい解が2^128未満にただ1つ存在することがわかっており, M1M2M3M4M5 > 2^128であるためk = 0
    /// t * M1M2M3M4 <= M1M2M3M4M5 - M1M2M3M4
    //             x_4 < M1M2M3M4
    // であるため, x_4 + t * M1M2M3M4 = x_5 < 2^128でオーバーフローしない
    x += t * M1M2M3M4;
    c[i] = x;
  }
  return c;
}
template<typename mint>
std::vector<mint> convolution_mod_large(const std::vector<mint> &a, const std::vector<mint> &b){
  if(a.empty() || b.empty()) return {};
  static constexpr int block_size = 1 << 22; // 998244353
  auto split = [](const std::vector<mint> &_a){
    std::vector<std::vector<mint>> res;
    int n = _a.size();
    for(int i = 0; i < n; i++){
      if(i % block_size == 0) res.push_back(std::vector<mint>{});
      res.back().push_back(_a[i]);
    }
    return res;
  };
  auto a2 = split(a);
  auto b2 = split(b);
  int n = (int)a2.size();
  int m = (int)b2.size();
  for(int i = 0; i < n; i++){
    a2[i].resize(block_size << 1);
    butterfly(a2[i]);
  }
  for(int i = 0; i < m; i++){
    b2[i].resize(block_size << 1);
    butterfly(b2[i]);
  }
  std::vector<mint> tmp(block_size << 1);
  std::vector<mint> res((n + m) * block_size);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      for(int k = 0; k < (block_size << 1); k++){
        tmp[k] = a2[i][k] * b2[j][k];
      }
      butterfly_inv(tmp);
      int offset = (i + j) * block_size;
      for(int k = 0; k < (block_size << 1); k++) res[offset + k] += tmp[k];
    }
  }
  mint iz = mint(block_size << 1).inv();
  int res_size = (int)a.size() + (int)b.size() - 1;
  res.resize(res_size);
  for(int i = 0; i < res_size; i++) res[i] *= iz;
  return res;
}
#endif