#ifndef _CONVOLUTION_H_
#define _CONVOLUTION_H_
#include "mod.hpp"
#include <algorithm>

template<typename mint>
void butterfly(std::vector<mint>& a){
  static constexpr int g = primitive_root<mint::mod()>;
  int n = int(a.size());
  int h = ceil_pow2(n);
  static bool first = true;
  static mint sum_e[30];  // sum_e[i] = ies[0] * ... * ies[i - 1] * es[i]
  if(first){
    first = false;
    mint es[30], ies[30];  // es[i]^(2^(2+i)) == 1
    int cnt2 = bsf(mint::mod() - 1);
    mint e = mint(g).pow((mint::mod() - 1) >> cnt2), ie = e.inv();
    for(int i = cnt2; i >= 2; i--){
      // e^(2^i) == 1
      es[i - 2] = e;
      ies[i - 2] = ie;
      e *= e;
      ie *= ie;
    }
    mint now = 1;
    for(int i = 0; i <= cnt2 - 2; i++){
      sum_e[i] = es[i] * now;
      now *= ies[i];
    }
  }
  for(int ph = 1; ph <= h; ph++){
    int w = 1 << (ph - 1), p = 1 << (h - ph);
    mint now = 1;
    for(int s = 0; s < w; s++){
      int offset = s << (h - ph + 1);
      for(int i = 0; i < p; i++){
        auto l = a[i + offset];
        auto r = a[i + offset + p] * now;
        a[i + offset] = l + r;
        a[i + offset + p] = l - r;
      }
      now *= sum_e[bsf(~(unsigned int)(s))];
    }
  }
}
template<typename mint>
void butterfly_inv(std::vector<mint>& a){
  static constexpr int g = primitive_root<mint::mod()>;
  int n = int(a.size());
  int h = ceil_pow2(n);
  static bool first = true;
  static mint sum_ie[30];  // sum_ie[i] = es[0] * ... * es[i - 1] * ies[i]
  if(first){
    first = false;
    mint es[30], ies[30];  // es[i]^(2^(2+i)) == 1
    int cnt2 = bsf(mint::mod() - 1);
    mint e = mint(g).pow((mint::mod() - 1) >> cnt2), ie = e.inv();
    for(int i = cnt2; i >= 2; i--){
      // e^(2^i) == 1
      es[i - 2] = e;
      ies[i - 2] = ie;
      e *= e;
      ie *= ie;
    }
    mint now = 1;
    for(int i = 0; i <= cnt2 - 2; i++){
      sum_ie[i] = ies[i] * now;
      now *= es[i];
    }
  }
  for(int ph = h; ph >= 1; ph--){
    int w = 1 << (ph - 1), p = 1 << (h - ph);
    mint inow = 1;
    for(int s = 0; s < w; s++){
      int offset = s << (h - ph + 1);
      for(int i = 0; i < p; i++){
        auto l = a[i + offset];
        auto r = a[i + offset + p];
        a[i + offset] = l + r;
        a[i + offset + p] = (unsigned long long)(mint::mod() + l.val() - r.val()) * inow.val();
      }
      inow *= sum_ie[bsf(~(unsigned int)(s))];
    }
  }
}
template<typename mint>
std::vector<mint> convolution_mod(std::vector<mint> a, std::vector<mint> b){
  int n = int(a.size()), m = int(b.size());
  if(!n || !m) return {};
  if(std::min(n, m) <= 60){
    if(n < m){
      std::swap(n, m);
      std::swap(a, b);
    }
    std::vector<mint> ans(n + m - 1);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        ans[i + j] += a[i] * b[j];
      }
    }
    return ans;
  }
  int z = 1 << ceil_pow2(n + m - 1);
  a.resize(z);
  butterfly(a);
  b.resize(z);
  butterfly(b);
  for(int i = 0; i < z; i++) a[i] *= b[i];
  butterfly_inv(a);
  a.resize(n + m - 1);
  mint iz = mint(z).inv();
  for(int i = 0; i < n + m - 1; i++) a[i] *= iz;
  return a;
}
template <unsigned int mod = 998244353, typename T>
std::vector<T> convolution_mod(const std::vector<T>& a, const std::vector<T>& b) {
  int n = int(a.size()), m = int(b.size());
  if(!n || !m) return {};
  using mint = static_modint<mod>;
  std::vector<mint> a2(n), b2(m);
  for(int i = 0; i < n; i++) a2[i] = mint(a[i]);
  for(int i = 0; i < m; i++) b2[i] = mint(b[i]);
  auto c2 = convolution_mod<mint>(std::move(a2), std::move(b2));
  std::vector<T> c(n + m - 1);
  for (int i = 0; i < n + m - 1; i++) c[i] = c2[i].val();
  return c;
}
template<typename mint>
std::vector<mint> square_mod(std::vector<mint> a){
  int n = int(a.size());
  if(!n) return {};
  if(n <= 60){
    std::vector<mint> ans(2 * n - 1);
    for(int i = 0; i < n; i++){
      for(int j = i + 1; j < n; j++){
        ans[i + j] += 2 * a[i] * a[j];
      }
      ans[2 * i] += a[i] * a[i];
    }
    return ans;
  }
  int z = 1 << ceil_pow2(2 * n - 1);
  a.resize(z);
  butterfly(a);
  for(int i = 0; i < z; i++) a[i] *= a[i];
  butterfly_inv(a);
  a.resize(2 * n - 1);
  mint iz = mint(z).inv();
  for (int i = 0; i < 2 * n - 1; i++) a[i] *= iz;
  return a;
}
template <unsigned int mod = 998244353, class T>
std::vector<T> square_mod(const std::vector<T> &a){
  int n = int(a.size());
  if(!n) return {};
  using mint = static_modint<mod>;
  std::vector<mint> a2(n);
  for(int i = 0; i < n; i++) a2[i] = mint(a[i]);
  auto c2 = square_mod<mint>(std::move(a2));
  std::vector<T> c(2 * n - 1);
  for(int i = 0; i < 2 * n - 1; i++) c[i] = c2[i].val();
  return c;
}

// ntt-friendlyでないmodで畳み込みを行う
template<typename mint>
std::vector<mint> convolution_int_mod(const std::vector<mint>& a, const std::vector<mint>& b){
  if(mint::mod() == 998244353) return convolution_mod<mint>(a, b);
  int n = int(a.size()), m = int(b.size());
  if(!n || !m) return {};
  if(std::min(n, m) <= 60){
    std::vector<mint> ans(n + m - 1, 0);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        ans[i + j] += a[i] * b[j];
      }
    }
    return ans;
  }
  static constexpr long long MOD1 = 167772161, MOD2 = 469762049, MOD3 = 1224736769;
  static constexpr long long M1M2 = MOD1 * MOD2, ix = inv_gcd(MOD1, MOD2).second, i3 = inv_gcd(MOD1 * MOD2, MOD3).second;
  std::vector<long long> a2(n), b2(m);
  for(int i = 0; i < n; i++) a2[i] = a[i].val();
  for(int i = 0; i < m; i++) b2[i] = b[i].val();
  auto c1 = convolution_mod<MOD1>(a2, b2);
  auto c2 = convolution_mod<MOD2>(a2, b2);
  auto c3 = convolution_mod<MOD3>(a2, b2);
  std::vector<mint> c(n + m - 1);
  int M1M2m = M1M2 % mint::mod();
  for(int i = 0; i < n + m - 1; i++){
    long long v = (((long long)c2[i] - c1[i]) * ix) % MOD2;
    if(v < 0) v += MOD2;
    long long xxv = c1[i] + MOD1 * v;
    v = ((c3[i] - (xxv % MOD3)) * i3) % MOD3;
    if(v < 0) v += MOD3;
    c[i] = mint(xxv + M1M2m * v);
  }
  return c;
}
// ntt-friendlyでないmodでa * aを計算する
template<typename mint>
std::vector<mint> square_int_mod(const std::vector<mint>& a){
  if(mint::mod() == 998244353) return square_mod<mint>(a);
  int n = int(a.size());
  if(!n) return {};
  if(n <= 60){
    std::vector<mint> ans(2 * n - 1);
    for(int i = 0; i < n; i++){
      for(int j = i + 1; j < n; j++){
        ans[i + j] += 2 * a[i] * a[j];
      }
      ans[2 * i] += a[i] * a[i];
    }
    return ans;
  }
  static constexpr long long MOD1 = 167772161, MOD2 = 469762049, MOD3 = 1224736769;
  static constexpr long long M1M2 = MOD1 * MOD2, ix = inv_gcd(MOD1, MOD2).second, i3 = inv_gcd(MOD1 * MOD2, MOD3).second;
  std::vector<long long> a2(n);
  for(int i = 0; i < n; i++) a2[i] = a[i].val();
  auto c1 = square_mod<MOD1>(a2);
  auto c2 = square_mod<MOD2>(a2);
  auto c3 = square_mod<MOD3>(a2);
  std::vector<mint> c(2 * n - 1);
  int M1M2m = M1M2 % mint::mod();
  for(int i = 0; i < 2 * n - 1; i++){
    long long v = (((long long)c2[i] - c1[i]) * ix) % MOD2;
    if(v < 0) v += MOD2;
    long long xxv = c1[i] + MOD1 * v;
    v = ((c3[i] - (xxv % MOD3)) * i3) % MOD3;
    if(v < 0) v += MOD3;
    c[i] = mint(xxv + M1M2m * v);
  }
  return c;
}
// @param 答えがlong longに収まる
std::vector<long long> convolution_ll(const std::vector<long long>& a, const std::vector<long long>& b){
  int n = int(a.size()), m = int(b.size());
  if (!n || !m) return {};
  if(std::min(n, m) <= 60){
    std::vector<long long> ans(n + m - 1, 0);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        ans[i + j] += a[i] * b[j];
      }
    }
    return ans;
  }
  // 2^24, 2^25, 2^26
  static constexpr unsigned long long MOD1 = 754974721, MOD2 = 167772161, MOD3 = 469762049;
  static constexpr unsigned long long M2M3 = MOD2 * MOD3, M1M3 = MOD1 * MOD3, M1M2 = MOD1 * MOD2, M1M2M3 = MOD1 * MOD2 * MOD3;
  static constexpr unsigned long long i1 = inv_gcd(MOD2 * MOD3, MOD1).second, i2 = inv_gcd(MOD1 * MOD3, MOD2).second, i3 = inv_gcd(MOD1 * MOD2, MOD3).second;
  auto c1 = convolution_mod<MOD1>(a, b);
  auto c2 = convolution_mod<MOD2>(a, b);
  auto c3 = convolution_mod<MOD3>(a, b);
  std::vector<long long> c(n + m - 1);
  for(int i = 0; i < n + m - 1; i++){
    unsigned long long x = 0;
    x += (c1[i] * i1) % MOD1 * M2M3;
    x += (c2[i] * i2) % MOD2 * M1M3;
    x += (c3[i] * i3) % MOD3 * M1M2;
    long long diff = c1[i] - safe_mod((long long)(x), (long long)(MOD1));
    if (diff < 0) diff += MOD1;
    static constexpr unsigned long long offset[5] = {0, 0, M1M2M3, 2 * M1M2M3, 3 * M1M2M3};
    x -= offset[diff % 5];
    c[i] = x;
  }
  return c;
}
// @param 答えがlong longに収まる
std::vector<long long> square_ll(const std::vector<long long>& a){
  int n = int(a.size());
  if (!n) return {};
  if(n <= 60){
    std::vector<long long> ans(2 * n - 1, 0);
    for(int i = 0; i < n; i++){
      for(int j = i + 1; j < n; j++){
        ans[i + j] += 2 * a[i] * a[j];
      }
      ans[2 * i] += a[i] * a[i];
    }
    return ans;
  }
  // 2^24, 2^25, 2^26
  static constexpr unsigned long long MOD1 = 754974721, MOD2 = 167772161, MOD3 = 469762049;
  static constexpr unsigned long long M2M3 = MOD2 * MOD3, M1M3 = MOD1 * MOD3, M1M2 = MOD1 * MOD2, M1M2M3 = MOD1 * MOD2 * MOD3;
  static constexpr unsigned long long i1 = inv_gcd(MOD2 * MOD3, MOD1).second, i2 = inv_gcd(MOD1 * MOD3, MOD2).second, i3 = inv_gcd(MOD1 * MOD2, MOD3).second;
  auto c1 = square_mod<MOD1>(a);
  auto c2 = square_mod<MOD2>(a);
  auto c3 = square_mod<MOD3>(a);
  std::vector<long long> c(2 * n - 1);
  for(int i = 0; i < 2 * n - 1; i++){
    unsigned long long x = 0;
    x += (c1[i] * i1) % MOD1 * M2M3;
    x += (c2[i] * i2) % MOD2 * M1M3;
    x += (c3[i] * i3) % MOD3 * M1M2;
    long long diff = c1[i] - safe_mod((long long)(x), (long long)(MOD1));
    if (diff < 0) diff += MOD1;
    static constexpr unsigned long long offset[5] = {0, 0, M1M2M3, 2 * M1M2M3, 3 * M1M2M3};
    x -= offset[diff % 5];
    c[i] = x;
  }
  return c;
}
#endif