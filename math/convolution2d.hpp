#ifndef _CONVOLUTION_2d_H_
#define _CONVOLUTION_2d_H_
#include "convolution.hpp"

template<typename mint> 
void butterfly2d(std::vector<std::vector<mint>> &a){
  int n = a.size();
  int m = a.empty() ? 0 : a[0].size();
  std::vector<std::vector<mint>> b(m, std::vector<mint>(n));
  for(int i = 0; i < n; i++) butterfly<mint>(a[i]);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++) b[j][i] = a[i][j];
  }
  for(int i = 0; i < m; i++) butterfly<mint>(b[i]);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++) a[i][j] = b[j][i];
  }
}
template<typename mint> 
void butterfly_inv2d(std::vector<std::vector<mint>> &a){
  int n = a.size();
  int m = a.empty() ? 0 : a[0].size();
  std::vector<std::vector<mint>> b(m, std::vector<mint>(n));
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++) b[j][i] = a[i][j];
  }
  for(int i = 0; i < m; i++) butterfly_inv<mint>(b[i]);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++) a[i][j] = b[j][i];
  }
  for(int i = 0; i < n; i++) butterfly_inv<mint>(a[i]);
}

// n * m配列の2Dfft
template<typename mint>
std::vector<std::vector<mint>> square_mod2d(std::vector<std::vector<mint>> a){
  int n = a.size();
  int m = a.empty() ? 0 : a[0].size();
  if(!n || !m) return {};
  int N = 1 << ceil_pow2(2 * n - 1);
  int M = 1 << ceil_pow2(2 * m - 1);
  a.resize(N);
  for(int i = 0; i < N; i++){
    assert(a[i].size() <= m);
    a[i].resize(M, 0);
  }
  butterfly2d<mint>(a);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < M; j++) a[i][j] *= a[i][j];
  }
  butterfly_inv2d<mint>(a);
  mint iz = (mint(N) * M).inv();
  a.resize(2 * n - 1);
  for(int i = 0; i < 2 * n - 1; i++){
    a[i].resize(2 * m - 1);
    for(int j = 0; j < 2 * m - 1; j++) a[i][j] *= iz;
  }
  return a;
}

template <unsigned int mod = 998244353, class T>
std::vector<std::vector<T>> square_mod2d(const std::vector<std::vector<T>> &a){
  int n = a.size();
  int m = a.empty() ? 0 : a[0].size();
  if(!n || !m) return {};
  using mint = static_modint<mod>;
  std::vector<std::vector<mint>> a2(n, std::vector<mint>(m));
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++) a2[i][j] = mint(a[i][j]);
  }
  auto c2 = square_mod2d<mint>(std::move(a2));
  std::vector<std::vector<T>> c(2 * n - 1, std::vector<T>(2 * m - 1));
  for(int i = 0; i < 2 * n - 1; i++){
    for(int j = 0; j < 2 * m - 1; j++) c[i][j] = c2[i][j].val();
  }
  return c;
}

template<typename mint>
std::vector<std::vector<mint>> square_int_mod2d(const std::vector<std::vector<mint>> & a){
  if(mint::mod() == 998244353) return square_mod2d<mint>(a);
  int n = a.size();
  int m = a.empty() ? 0 : a[0].size();
  if(!n || !m) return {};
  static constexpr long long MOD1 = 167772161, MOD2 = 469762049, MOD3 = 1224736769;
  static constexpr long long M1M2 = MOD1 * MOD2, ix = inv_gcd(MOD1, MOD2).second, i3 = inv_gcd(MOD1 * MOD2, MOD3).second;
  std::vector<std::vector<int>> a2(n, std::vector<int>(m));
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++) a2[i][j] = a[i][j].val();
  }
  auto c1 = square_mod2d<MOD1>(a2);
  auto c2 = square_mod2d<MOD2>(a2);
  auto c3 = square_mod2d<MOD3>(a2);
  std::vector<std::vector<mint>> c(2 * n - 1, std::vector<mint>(2 * m - 1));
  int M1M2m = M1M2 % mint::mod();
  for(int i = 0; i < 2 * n - 1; i++){
    for(int j = 0; j < 2 * m - 1; j++){
      long long v = (((long long)c2[i][j] - c1[i][j]) * ix) % MOD2;
      if(v < 0) v += MOD2;
      long long xxv = (long long)c1[i][j] + MOD1 * v;
      v = (((long long)c3[i][j] - (xxv % MOD3)) * i3) % MOD3;
      if(v < 0) v += MOD3;
      c[i][j] = mint(xxv + M1M2m * v);
    }
  }
  return c;
}
#endif
