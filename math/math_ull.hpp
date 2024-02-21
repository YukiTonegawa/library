#ifndef _MATH_ULL_H_
#define _MATH_ULL_H_
#include <vector>
#include <limits>
#include <numeric>
#include <cassert>
namespace math_ull{
  using ull = unsigned long long;
  static constexpr ull inf = std::numeric_limits<ull>::max();

  // min(a + b, inf)
  ull add_ull(ull a, ull b){
    return inf - a < b ? inf : a + b;
  }
  // min(a * b, inf)
  ull mul_ull(ull a, ull b){
    __uint128_t res = (__uint128_t)a * b;
    return res > inf ? inf : res;
  }
  // min(a^k, inf)
  ull pow_ull(ull a, ull k){
    ull res = 1;
    while(k){
      if(k & 1) res = mul_ull(res, a);
      a = mul_ull(a, a);
    }
    return res;
  }
  // min(n!, inf)
  ull fac_ull(ull n){
    assert(n >= 0);
    if(n > 20) return inf;
    static std::vector<ull> res;
    if(res.empty()){
      res.resize(21, 1);
      for(int i = 0; i < 21; i++){
        for(int j = 1; j <= i; j++) res[i] *= j;
      }
    }
    return res[n];
  }
  // min(nCk, inf)
  // 連続するd数の積はdで必ず割り切れる
  ull comb_ull(ull n, ull k){
    if(n < 0 || k < 0 || n < k) return 0;
    if(n - k < k) k = n - k;
    __uint128_t res = 1;
    for(ull d = 1; d <= k; d++){
      res *= n--;
      res /= d;
      if(res > inf) return inf;
    }
    return res;
  }
  // min(nPk, inf)
  ull perm_ull(ull n, ull k){
    if(n < 0 || k < 0 || n < k) return 0;
    __uint128_t res = 1;
    for(ull x = n; x > n - k; x--){
      res *= x;
      if(res > inf) return inf;
    }
    return res;
  }
};
#endif