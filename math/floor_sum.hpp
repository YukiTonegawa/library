#ifndef _FLOOR_SUM_H_
#define _FLOOR_SUM_H_

#include "integer.hpp"

template<typename T>
T __mod(T a, T mod){
  a %= mod;
  return a < 0 ? a + mod : a;
}
// ∑{i, 0 <= i < n} floor((ai + b) / c)
// n, a, b, c <= 10^9
long long floor_sum(long long n, long long a, long long b, long long c){
  long long ans = 0;
  if(a >= c){
    ans += (n - 1) * n * (a / c) / 2;
    a %= c;
  }
  if(b >= c){
    ans += n * (b / c);
    b %= c;
  }
  // 変形:
  // https://rsk0315.hatenablog.com/entry/2020/12/13/231307
  // https://atcoder.jp/contests/practice2/editorial/579
  long long y_max = a * n + b;
  if(y_max >= c) ans += floor_sum(y_max / c, c, y_max % c, a);
  return ans;
}
// ∑{i, 0 <= i < n} ceil((ai + b) / c)
// n, a, b, c <= 10^9
long long ceil_sum(long long n, long long a, long long b, long long c){
  return floor_sum(n, a, b + c - 1, c);
}

template<typename T>
struct __fs_struct{
  T f, g, h;
  __fs_struct(T _f = 0, T _g = 0, T _h = 0): f(_f), g(_g), h(_h){}
};
template<typename T>
__fs_struct<T> __fs_query(T n, long long a, long long b, long long c){
  if(!n) return {0, 0, 0};
  long long ac = a / c, bc = b / c;
  T nc2 = n * (n - 1) / 2;
  if(!a){
    return {(T)bc * n, bc * nc2, (T)bc * n * bc};
  }
  T f, g, h;
  __fs_struct<T> next;
  if(a >= c || b >= c){
    T sqsum = nc2 * (2 * n - 1) / 3;
    next = __fs_query<T>(n, a % c, b % c, c);
    f = next.f + nc2 * ac + n * bc;
    g = next.g + ac * sqsum + bc * nc2;
    h = next.h + ac * (2 * next.g + ac * sqsum) + bc * (2 * next.f + (T)bc * n + (T)2 * ac * nc2);
    return {f, g, h};
  }
  long long m = (a * (n - 1) + b) / c;
  next = __fs_query<T>(m, c, c - b - 1, a);
  f = (T)m * (n - 1) - next.f;
  g = ((T)m * n * (n - 1) - next.h - next.f) / 2;
  h = (T)(n - 1) * m * (m + 1) - 2 * next.g - 2 * next.f - f;
  return {f, g, h};
}
// ∑{i, 0 <= i < n} i * floor((ai + b) / c)
// n, a, b, c <= 10^9
template<typename T>
T floor_sum_slope(long long n, long long a, long long b, long long c){
  assert(0 <= n && 0 <= a && 0 <= b && 0 < c);
  return __fs_query<T>(n, a, b, c).g;
}
// ∑{i, 0 <= i < n} floor((ai + b) / c) ^ 2
// n, a, b, c <= 10^9
template<typename T>
T floor_sum_square(long long n, long long a, long long b, long long c){
  assert(0 <= n && 0 <= a && 0 <= b && 0 < c);
  return __fs_query<T>(n, a, b, c).h;
}
// min((ai + b) % mod)  0 <= i < r
long long min_linear_mod(long long l, long long r, long long a, long long b, long long mod){
  assert(0 < mod && l <= 0 && l < r);
  a = __mod(a, mod);
  b = __mod(b + a * l, mod);
  if(a == 0) return b;
  long long upper = floor_sum(r, a, b, mod), L = 0, R = b + 1;
  // 定数倍高速化
  long long tmp_val = b, tmp_idx = 0;
  for(int i = 0; i < 10; i++){
    long long sa = mod - tmp_val;
    tmp_idx += (sa + a - 1) / a;
    tmp_val = (tmp_idx * a + b) % mod;
    if(tmp_idx >= r) break;
    R = std::min(R, tmp_val + 1);
  }
  while(R - L > 1){
    long long mid = (L + R) / 2;
    if(upper == floor_sum(r, a, b - mid, mod)) L = mid;
    else R = mid;
  }
  return L;
}
// max((ai + b) % mod)  0 <= i < r
long long max_linear_mod(long long l, long long r, long long a, long long b, long long mod){
  return mod - 1 - min_linear_mod(l, r, -a, -(b + 1), mod);
}

// 0 <= i < nかつ0 <= i % mod < rとなるようなiの数
// 0 <= n, r
// 0 < mod
long long count_range_mod(long long n, long long r, long long mod){
  if(!n || !r) return 0;
  assert(0 < n && 0 < mod && 0 <= r);
  if(r >= mod) return n;
  n--;
  long long d = n / mod, ret = r * d;
  n %= mod;
  ret += std::min(n + 1, r);
  return ret;
}
// 0 <= i < n　かつ ai + bがmの倍数になるiの数
// 0 <= n
// 0 < m
long long count_multiple_linear(long long n, long long a, long long b, long long m){
  assert(0 < m);
  a = __mod(a, m);
  b = __mod(b, m);
  auto [g, x] = inv_ext_gcd(a, m);
  if((m - b) % g) return 0;
  long long k = m / g; // ak ≡ 0 (mod m)となる最小の自然数
  long long y = (((m - b) / g) * x) % k; // ay + b ≡ 0
  if(y >= n) return 0;
  return 1 + (n - 1 - y) / k;
}
// l <= i　かつ ai + bがmの倍数になる最小のi
// 存在しない場合は-1
// 0 <= l
// 0 < m
long long next_multiple_linear(long long l, long long a, long long b, long long m){
  assert(0 < m);
  a = __mod(a, m);
  b = __mod(b, m);
  auto [g, x] = inv_ext_gcd(a, m);
  if((m - b) % g) return -1;
  long long k = m / g; // ak ≡ 0 (mod m)となる最小の自然数
  long long y = (((m - b) / g) * x) % k; // ay + b ≡ 0
  if(y >= l) return y;
  return y + ((l - y + k - 1) / k) * k;
}
// 0 <= i < nについて, 0 <= (ai + b) % mod < rとなるようなiの数
// 0 < mod
// 0 <= n, r
long long count_range_mod_linear(long long n, long long r, long long a, long long b, long long mod){
  assert(0 <= n && 0 <= r);
  a = __mod(a, mod);
  b = __mod(b, mod);
  if(r >= mod) return n;
  long long x = floor_sum(n, a, b + mod, mod);
  long long y = floor_sum(n, a, b + mod - r, mod);
  return x - y;
}
// l <= iかつ 0 <= (ai + b) % mod < rとなるような最小の整数i
// 存在しないなら-1
// 0 <= r, l
// 0 < mod
long long next_range_mod_linear(long long l, long long r, long long a, long long b, long long mod){
  a = __mod(a, mod);
  b = __mod(b + a * l, mod);
  assert(0 <= l && 0 <= r && 0 < mod);
  if(r >= mod) return l;
  if(a == 0) return b >= r ? -1 : l;
  long long g = gcd(a, mod);
  if(b % g >= r) return -1;
  long long L = -1, R = mod / g;
  while(R - L > 1){
    long long M = (L + R) / 2;
    if(floor_sum(M + 1, a, b + mod, mod) == floor_sum(M + 1, a, b + mod - r, mod)) L = M;
    else R = M;
  }
  return l + R;
}
// 0 <= i < n, 0 <= j < m かつ ai + bjがcの倍数
long long count_multiple_linear2(long long n, long long m, long long a, long long b, long long c){
  if(c == 0) return n * m; // 常に倍数
  if(c < 0) a *= -1, b *= -1, c *= -1;
  a = __safe_mod(a, c);
  b = __safe_mod(b, c);
  n--, m--;
  long long g_all = gcd(a, gcd(b, c));
  a /= g_all;
  b /= g_all;
  c /= g_all;
  long long g = gcd(a, c);
  long long a2 = a / g;
  long long c2 = c / g;
  long long b2 = b * inv_ext_gcd(a2, c2).second % c2;
  long long q = n / c2, r = n % c2;
  return count_range_mod_linear(m / g + 1, r + 1, -b2, 0, c2) + count_range_mod_linear(m / g + 1, c2, -b2, 0, c2) * q;
}
// [xi + y, xi + z) (0 <= y < z <= x)と
// [aj + b, aj + c) (0 <= b < c <= a)が交差する最小の非負整数i
// 解がない場合は-1
long long find_intersection_of_periodic_segments(long long x, long long y, long long z, long long a, long long b, long long c){
  assert(x > 0 && a > 0);
  assert(0 <= y && y < z && z <= x);
  assert(0 <= b && b < c && c <= a);
  long long ag = a / gcd(x, a);
  long long L = 0, R = ag + 1;
  while(R - L > 1){
    long long M = (L + R) / 2;
    if(floor_sum(M, x, z - 1 + a - b, a) == floor_sum(M, x, y + a - c, a)) L = M;
    else R = M;
  }
  if(R == ag + 1) return -1;
  return L;
}
/*
// 0 <= a < cのとき、[ai + b, ci + d]にeの倍数(負の数も含む)が含まれないような自然数iの数(有限)
long long simple_math3(long long a, long long b, long long c, long long d, long long e){
  assert(a >= 0 && c >= 0 && a < c && a + b <= c + d);
  if(std::min(b, d) < 0){
    long long add = (abs(std::min(b, d) + e - 1) / e) * e;
    b += add;
    d += add;
  }
  long long k = (e - 2 - (d - b)) / (c - a);
  if(k <= 0) return 0;
  return k - floor_sum(k, c, c + d, e) + floor_sum(k, a, a + b - 1, e);
}
*/

// (ai + b) % modのprefix_minumumを与える等差数列
std::vector<std::tuple<long long, long long, long long>> prefix_minimum(long long a, long long b, long long mod){
  assert(0 < mod);
  a = __mod(a, mod);
  b = __mod(b, mod);
  if(a == 0 || b == 0) return {{0, 1, 1}};
  long long p = 0, q = 1, r = 1, s = 0, g = gcd(a, mod);
  a = (mod - a) / g;
  mod /= g;
  long long X = 0, Y = b;
  std::vector<std::tuple<long long, long long, long long>> ret;
  auto step = [&](){
    long long P = p + r, Q = q + s;
    long long R = mod * r - a * s;
    long long S = a * q - mod * p;
    if(a * Q <= mod * P){
      long long x = R / S;
      r += p * x;
      s += q * x;
    }else{
      std::swap(R, S);
      long long x = R / S;
      long long y = std::max(0LL, (g * R - Y + g * S - 1) / (g * S));
      long long c = g * (a * (q + s * y) - mod * (p + r * y));
      if(c == 0){
        s = mod;
        return ;
      }
      if(y <= x && 0 < c && c <= Y){
        long long k = Y / c;
        ret.push_back({X, k + 1, q + s * y});
        Y -= c * k;
        X += (q + s * y) * k;
        p += r * std::min(x, y + 1);
        q += s * std::min(x, y + 1);
      }else{
        p += r * x;
        q += s * x;
      }
    }
  };
  while(q < mod && s < mod) step();
  if(q < mod && a * q - mod * p > 0){
    long long c = g * (a * q - mod * p);
    long long k = Y / c;
    ret.push_back({X, k + 1, q});
  }
  return ret;
}
// min((ai + b) % mod)  0 <= i < r
long long min_linear_mod_fast(long long l, long long r, long long a, long long b, long long mod){
  assert(0 < mod && l <= 0 && l < r);
  a = __mod(a, mod);
  b = __mod(b + a * l, mod);
  r -= l;
  if(a == 0) return b;
  auto pm = prefix_minimum(a, b, mod);
  int idx = std::lower_bound(pm.begin(), pm.end(), std::tuple<long long, long long, long long>{r, 0, 0}) - pm.begin() - 1;
  auto [x0, k, d] = pm[idx];
  long long x = x0 + d * std::min(k - 1, (r - 1 - x0) / d);
  return (a * x + b) % mod;
}
// max((ai + b) % mod)  0 <= i < r
long long max_linear_mod_fast(long long l, long long r, long long a, long long b, long long mod){
  return mod - 1 - min_linear_mod_fast(l, r, -a, -(b + 1), mod);
}

// l <= iかつ 0 <= (ai + b) % mod < rとなるような最小の整数i
// 存在しないなら-1
// 0 <= r, l
// 0 < mod
long long next_range_mod_linear_fast(long long l, long long r, long long a, long long b, long long mod){
  a = __mod(a, mod);
  b = __mod(b + a * l, mod);
  assert(0 <= l && 0 <= r && 0 < mod);
  if(r >= mod) return l;
  if(a == 0) return b >= r ? -1 : l;
  auto pm = prefix_minimum(a, b, mod);
  for(auto [x0, k, d] : pm){
    long long x_max = x0 + d * (k - 1);
    if((a * x_max + b) % mod >= r) continue;
    long long y0 = (a * x0 + b) % mod;
    if(y0 < r) return l + x0;
    long long diff = mod - (d * a % mod); // xがd増加した時の減少量
    assert(diff != mod);
    long long ret = x0 + d * ((y0 - r + diff) / diff);
    assert(ret <= x_max);
    return l + ret;
  }
  return -1;
}

// [xi + y, xi + z) (0 <= y < z <= x)と
// [aj + b, aj + c) (0 <= b < c <= a)が交差する最小の非負整数i
// 解がない場合は-1
long long find_intersection_of_periodic_segments_fast(long long x, long long y, long long z, long long a, long long b, long long c){
  assert(x > 0 && a > 0);
  assert(0 <= y && y < z && z <= x);
  assert(0 <= b && b < c && c <= a);
  // [0, xi + z - 1 - b)と[0, xi + y - c)で含まれるaの倍数の数が異なる < - > (xi + (z - 1 - b)) % aがdiff未満になる最小のi
  long long diff = (z - 1 - b) - (y - c);
  assert(diff > 0);
  return next_range_mod_linear_fast(0, diff, x, z - 1 - b, a);
}
// ∑(0 <= i < n) (ai + b) % mod
long long sum_linear_mod(long long n, long long a, long long b, long long mod){
  a = __mod(a, mod);
  b = __mod(b, mod);
  // (ai + b) % mod = (ai + b) - mod * floor((ai + b) / mod)
  return arithmetic_sum<__int128_t>(a, b, 0, n) - (__int128_t)floor_sum(n, a, b, mod) * mod;
}
// !
// 数列Sのi項目 S_i := ai + bとしたとき
// ∑{0 <= i < n} (S_i % mod) ^ 2
__int128_t sum_linear_mod_square(long long n, long long a, long long b, long long mod){
  __int128_t x = square_sum<__int128_t>(a, 0, n) + arithmetic_sum<__int128_t>((__int128_t)2 * a * b, b * b, 0, n);
  __int128_t y = (__int128_t)mod * a * floor_sum_slope<__int128_t>(n, a, b, mod) + (__int128_t)mod * b * floor_sum(n, a, b, mod);
  __int128_t z = (__int128_t)mod * mod * floor_sum_square<__int128_t>(n, a, b, mod);
  return x - 2 * y + z;
}

// !
// 数列Sのi項目 S_i := ai + bとしたとき
// ∑{0 <= i < n} (S_i % mod) * S_i
__int128_t sum_linear_weighted_with_mod(long long n, long long a, long long b, long long mod){
  __int128_t x = square_sum<__int128_t>(a, 0, n) + arithmetic_sum<__int128_t>((__int128_t)2 * a * b, b * b, 0, n);
  __int128_t y = (__int128_t)mod * a * floor_sum_slope<__int128_t>(n, a, b, mod);
  __int128_t z = (__int128_t)mod * b * floor_sum(n, a, b, mod);
  return x - y - z;
}

// !
// i項目がai + bの等差数列の要素のうち, 以下の条件を満たすものの和
// 0 <= i < n
// 0 <= (ci + d) % mod < r
long long sum_linear_extracted_by_range_mod_linear(long long n, long long r, long long a, long long b, long long c, long long d, long long mod){
  __int128_t x = floor_sum_slope<__int128_t>(n, c, d + mod, mod) - floor_sum_slope<__int128_t>(n, c, d + mod - r, mod);
  long long y = floor_sum(n, c, d + mod, mod) - floor_sum(n, c, d + mod - r, mod);
  return a * x + b * y;
}
// !
// i項目がfloor((ai + b) / c)の数列の要素のうち, 以下の条件を満たすものの和
// 0 <= i < n
// 0 <= (ai + b) % mod < r
long long floor_sum_extracted_by_range_mod_linear(long long n, long long r, long long a, long long b, long long c, long long mod){
  __int128_t x = floor_sum_square<__int128_t>(n, mod * a, mod * b + c * mod, mod * c);
  __int128_t y = floor_sum_square<__int128_t>(n, mod * a, mod * b + c * (mod - r), mod * c);
  long long count = floor_sum(n, a, b + mod, mod) - floor_sum(n, a, b + mod - r, mod); // 条件を満たす項の数
  return (x - y + count) / 2;
}
#endif
