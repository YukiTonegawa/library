#ifndef _INTEGER_H_
#define _INTEGER_H_
#include <cassert>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <cmath>
#include <tuple>

// ak + b のl <= k < rにおける和
template<typename T = long long>
T arithmetic_sum(T a, T b, T l, T r){
  return a * (r * (r - 1) - l * (l - 1)) / 2 + b * (r - l);
}
// a * k^2の　l <= k < rにおける和
template<typename T = long long>
T square_sum(T a, T l, T r){
  return (r * (r - 1) * (2 * r - 1) - l * (l - 1) * (2 * l - 1)) / 6 * a;
}
// a^k * b のl <= k < rにおける和
// b * (a^r - a^l) / (a - 1)
template<typename T>
T geometric_sum(T a, T b, T l, T r){
  if(a == 1) return 0;
  return b * (a.pow(r) - a.pow(l)) / (a - 1);
}
// @param 0 <= b, a^bがオーバーフローしない
long long llpow(long long a, long long b){
  assert(0 <= b);
  long long ret = 1;
  while(b){
    if(b & 1) ret *= a;
    a *= a;
    b >>= 1;
  }
  return ret;
}
// a^b, オーバーフローする場合numeric_limits<ull>
// @param 0 <= a, b
unsigned long long ullpow_of(unsigned long long a, unsigned long long b){
  static constexpr unsigned long long inf = std::numeric_limits<unsigned long long>::max();
  unsigned long long ret = 1;
  while(b){
    if(b & 1){
      if(ret >= inf / a) return inf;
      ret *= a;
    }
    a = (a >= inf / a ? inf : a * a);
    b >>= 1;
  }
  return ret;
}

long long gcd(long long _a, long long _b){
  uint64_t a = abs(_a), b = abs(_b);
  if(a == 0) return b;
  if(b == 0) return a;
  int shift = __builtin_ctzll(a | b);
  a >>= __builtin_ctzll(a);
  do{
    b >>= __builtin_ctzll(b);
    if(a > b) std::swap(a, b);
    b -= a;
  }while(b);
  return (a << shift);
}
// 64bitに収まらない可能性あり
template<typename T>
long long lcm(T a, T b){
  return __int128_t(a) * b / gcd(a, b);
}

// min(lcm(a, b), inf)
template<typename T>
T lcm_limited(T a, T b, T inf){
  if(a >= inf || b >= inf) return inf;
  b /= gcd(a, b);
  return (inf / a) < b ? inf : a * b;
}

std::tuple<long long, long long, long long> extgcd(long long a, long long b){
  long long x, y;
  for(long long u = y = 1, v = x = 0; a;){
    long long q = b / a;
    std::swap(x -= q * u, u);
    std::swap(y -= q * v, v);
    std::swap(b -= q * a, a);
  }
  return {x, y, b};//return {x, y, gcd(a, b)} s.t. ax + by = gcd(a, b)
}
constexpr long long __safe_mod(long long x, long long m){
  x %= m;
  if (x < 0) x += m;
  return x;
}
// @param b `1 <= b`
// @return pair(g, x) s.t. g = gcd(a, b), xa = g (mod b), 0 <= x < b/g
constexpr std::pair<long long, long long> inv_ext_gcd(long long a, long long b) {
  a = __safe_mod(a, b);
  if (a == 0) return {b, 0};
  long long s = b, t = a;
  long long m0 = 0, m1 = 1;
  while (t){
    long long u = s / t;
    s -= t * u;
    m0 -= m1 * u;
    auto tmp = s;
    s = t;
    t = tmp;
    tmp = m0;
    m0 = m1;
    m1 = tmp;
  }
  if (m0 < 0) m0 += b / s;
  return {s, m0};
}
// {lcm, ans}
std::pair<long long, long long> crt(const std::vector<long long>& r, const std::vector<long long>& m){
  assert(r.size() == m.size());
  int n = int(r.size());
  // Contracts: 0 <= r0 < m0
  long long r0 = 0, m0 = 1;
  for(int i = 0; i < n; i++){
    assert(1 <= m[i]);
    long long r1 = __safe_mod(r[i], m[i]), m1 = m[i];
    if(m0 < m1){
      std::swap(r0, r1);
      std::swap(m0, m1);
    }
    if(m0 % m1 == 0){
      if(r0 % m1 != r1) return {0, 0};
      continue;
    }
    long long g, im;
    std::tie(g, im) = inv_ext_gcd(m0, m1);
    long long u1 = (m1 / g);
    if ((r1 - r0) % g) return {0, 0};
    long long x = (r1 - r0) / g % u1 * im % u1;
    r0 += x * m0;
    m0 *= u1;
    if (r0 < 0) r0 += m0;
  }
  return {r0, m0};
}

// ai + bj = cの li <= i < ri かつ lj <= j < rjを満たす整数解
// 制約: a, b, c, li, ri, lj, riの絶対値が10^18以下(負でもいい)
//      a, b != 0
long long count_solution(long long a, long long b, long long c, long long li, long long ri, long long lj, long long rj){
  if(li >= (ri--) || lj >= (rj--)) return 0;
  if(a < 0){
    a *= -1, li *= -1, ri *= -1;
    std::swap(li, ri);
  }
  if(b < 0){
    b *= -1, lj *= -1, rj *= -1;
    std::swap(lj, rj);
  }
  assert(a && b);
  auto [i, j, g] = extgcd(a, b);
  if(c % g != 0) return 0;
  a /= g, b /= g, c /= g;
  __int128_t i2 = (__int128_t)i * c;
  __int128_t j2 = (__int128_t)j * c;
  // 任意の整数kに対してa(i2 + bk) + b(j2 - ak) = cを満たす
  __int128_t lk1 = i2 >= li ? -(i2 - li) / b : (li - i2 + b - 1) / b;
  __int128_t rk1 = i2 >= ri ? -(i2 - ri + b - 1) / b : (ri - i2) / b;
  __int128_t lk2 = j2 >= rj ? (j2 - rj + a - 1) / a : -(rj - j2) / a;
  __int128_t rk2 = j2 >= lj ? (j2 - lj) / a : -(lj - j2 + a - 1) / a;
  lk1 = std::max(lk1, lk2);
  rk1 = std::min(rk1, rk2);
  if(lk1 > rk1) return 0;
  return rk1 - lk1 + 1;
}

// 2^k >= xとなる最小の2^k
uint64_t ceil_2_pow(uint64_t x){
  static uint64_t INF = std::numeric_limits<uint64_t>::max();
  uint64_t ret = uint64_t(1) << (63 - __builtin_clzll(x));
  if(ret == x) return x;
  assert(ret <= (INF >> 1));
  return ret << 1;
}
#endif