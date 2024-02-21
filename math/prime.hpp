#ifndef _PRIME_H_
#define _PRIME_H_
#include <vector>
#include <random>
#include <cassert>
#include <algorithm>
#include <unordered_map>
#include "math_data_structure/montgomery_reduction.hpp"
namespace prime_sieve{
  std::vector<int> primes, min_factor;// 素数, 各数を割り切る最小の素数
  // O(MAX_N loglog MAX_N)
  // [1, MAX_N]を扱えるように初期化
  void init(int MAX_N){
    min_factor.resize(MAX_N + 1, -1);
    for(int i = 2; i <= MAX_N; i++){
      if(min_factor[i] == -1){
        primes.push_back(i);
        min_factor[i] = i;
      }
      for(int p : primes){
        if((long long)p * i > MAX_N || p > min_factor[i]) break;
        min_factor[p * i] = p;
      }
    }
  }
  bool is_prime(int n){
    assert(n < min_factor.size());
    return n == min_factor[n];
  }
  // {{素因数, 数}}, O(log n)
  std::vector<std::pair<int, int>> factorize(int n){
    assert(n < min_factor.size());
    std::vector<std::pair<int, int>> ret;
    while(n > 1){
      int cnt = 0, f = min_factor[n];
      while(n % f == 0){
        n /= f;
        cnt++;
      }
      ret.push_back({f, cnt});
    }
    return ret;
  }
  // 約数列挙, O(√n)
  std::vector<int> divisor(int n){
    auto p = factorize(n);
    std::vector<std::vector<int>> x;
    for(int i = 0; i < p.size(); i++){
      x.push_back(std::vector<int>{1});
      for(int j = 0; j < p[i].second; j++) x[i].push_back(x[i][j] * p[i].first);
    }
    int l = 0, r = 1;
    std::vector<int> ret{1};
    for(int i = 0; i < x.size(); i++){
      for(auto e : x[i]){
        for(int j = l; j < r; j++){
          ret.push_back(ret[j] * e);
        }
      }
      l = r;
      r = ret.size();
    }
    return std::vector<int>(ret.begin() + l, ret.end());
  }
  // O(logN)
  unsigned long long totient_function(unsigned long long n){
    unsigned long long res = n;
    int prev = -1;
    while(n > 1){
      if(min_factor[n] > prev){
        res -= res / min_factor[n];
        prev = min_factor[n];
      }
      n /= min_factor[n];
    }
    return res;
  }
  int mobius_function(int x){
    int pcnt = 0;
    while(x > 1){
      int y = x / min_factor[x];
      if(min_factor[x] == min_factor[y]) return 0;
      x = y;
      pcnt++;
    }
    return pcnt % 2 == 0 ? 1 : -1;
  }
};

bool _miller_rabin_mr(unsigned long long n, const montgomery_reduction_64bit &mr){
  static std::vector<int> small_p{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
  static std::vector<unsigned long long> A{2, 325, 9375, 28178, 450775, 9780504, 1795265022};
  static std::vector<unsigned long long> B{2, 7, 61};
  if(n <= 1) return false;
  if(n <= 50){
    for(int l = n < 20 ? 0 : 8, r = n < 20 ? 8 : 15; l < r; l++) if(small_p[l] == n) return true;
    return false;
  }
  if(!(n & 1)) return false;
  unsigned long long d = n - 1;
  unsigned long long one = mr.one(), mone = mr.generate(n - 1);
  d >>= __builtin_ctzll(d);
  for(unsigned long long a : (n >> 32) ? A : B){
    if(a % n == 0) continue;
    unsigned long long d2s = d; // d * 2^s, 0 <= s <= (n - 1)が2で割れる回数
    unsigned long long y = mr.pow_safe(mr.generate(a), d);
    while(d2s != n - 1 && y != one && y != mone){
      y = mr.mul_safe(y, y);
      d2s <<= 1;
    }
    if(y != mone && !(d2s & 1)) return false;
  }
  return true;
}
bool miller_rabin_mr(unsigned long long n){
  if(n % 2 == 0) return n == 2 ? true : false;
  montgomery_reduction_64bit mr(n);
  return _miller_rabin_mr(n, mr);
}
// https://en.wikipedia.org/wiki/Binary_GCD_algorithm
unsigned long long binary_gcd(unsigned long long a, unsigned long long b){
  if(!a || !b) return !a ? b : a;
  int shift = __builtin_ctzll(a | b); // [1] gcd(2a', 2b') = 2 * gcd(a', b')
  a >>= __builtin_ctzll(a);
  do{
    // if b is odd
    // gcd(2a', b) = gcd(a', b), if a = 2a'(a is even)
    // gcd(a, b) = gcd(|a - b|, min(a, b)), if a is odd
    b >>= __builtin_ctzll(b); // make b odd
    if(a > b) std::swap(a, b);
    b -= a;
  }while(b); // gcd(a, 0) = a
  return a << shift; // [1]
}
unsigned long long generate_random_prime(unsigned long long min_n = 2, unsigned long long max_n = ~0ULL){
  std::random_device seed_gen;
  std::mt19937_64 engine(seed_gen());
  __uint128_t len = max_n - min_n + 1;
  // https://en.wikipedia.org/wiki/Prime_number_theorem
  while(true){
    unsigned long long a = engine() % len + min_n;
    if(miller_rabin_mr(a)){
      return a;
    }
  }
}
namespace rho_factorization{
  unsigned long long rho(unsigned long long n){
    static std::vector<int> small_p{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

    for(int sp : small_p) if(n % sp == 0) return sp; // n < 50

    montgomery_reduction_64bit mr(n);
    if(_miller_rabin_mr(n, mr)) return n;

    auto try_factorize = [n, mr](unsigned long long c){
      c = mr.generate(c);
      auto f = [mr, c](unsigned long long mx){
        return mr.add(mr.mul(mx, mx), c);
      };
      unsigned long long m = 1LL << ((64 - __builtin_clzll(n)) / 8);
      unsigned long long y = n, r = 1, q = 1;
      unsigned long long x, g, k, ys;
      do{
        x = y;
        y = mr.generate(y);
        for(int i = 0; i < r; i++) y = f(y);
        y = mr.reduce(y);

        k = 0;
        while(k < r && g == 1){
          q = mr.generate(q);
          y = mr.generate(y);
          ys = y;
          for(int i = 0; i < std::min(m, r - k); i++){
            y = f(y);
            unsigned long long z = mr.reduce(y);
            q = mr.mul(q, mr.generate(x > z ? x - z : z - x));
          }
          y = mr.reduce(y);
          g = binary_gcd(mr.reduce(q), n);
          k += m;
        }
        r <<= 1;
      }while(g == 1);
      if(g == n){
        do{
          ys = f(ys);
          unsigned long long z = mr.reduce(ys);
          g = binary_gcd(x > z ? x - z : z - x, n);
        }while(g == 1);
      }
      return g; // g == n when failure
    };
    unsigned long long c = 1, res = n;
    do{
      res = try_factorize(c);
      // c = generate_random_prime(2, n - 1);
      c = (c + 1) % n;
    }while(res == n);
    return res;
  }
  std::vector<unsigned long long> factorize(unsigned long long n){
    if(n <= 1) return {};
    unsigned long long x = rho(n);
    if(x == n) return {x};
    auto l = factorize(x);
    auto r = factorize(n / x);
    l.insert(l.end(), r.begin(), r.end());
    return l;
  }
  // {素数, 個数}の形で返す
  std::vector<std::pair<unsigned long long, int>> factorize2(unsigned long long n){
    auto p = factorize(n);
    sort(p.begin(), p.end());
    std::vector<std::pair<unsigned long long, int>> ret;
    for(int i : p){
      if(ret.empty() || ret.back().first != i) ret.push_back({i, 1});
      else ret.back().second++;
    }
    return ret;
  }
  // 素因数の集合(重複なし, ソート済)を返す
  std::vector<unsigned long long> prime_factor(unsigned long long n){
    auto p = factorize(n);
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    return p;
  }
  // 10^18以下の高度合成数 897612484786617600の約数が103680個なので全列挙して良さそう
  std::vector<unsigned long long> divisor(unsigned long long n){
    auto p = factorize(n);
    std::sort(p.begin(), p.end());

    std::vector<std::pair<unsigned long long, int>> x;

    for(int i = 0; i < p.size(); i++){
      if(!i || p[i] != p[i - 1]) x.push_back({p[i], 1});
      else x.back().second++;
    }
    int sz = 1;
    for(auto [p_cur, cnt] : x) sz *= cnt + 1;

    std::vector<unsigned long long> res(sz);
    res[0] = 1;
    int r_prev = 1, r = 1;
    for(auto [p_cur, cnt] : x){
      unsigned long long ppow = 1;
      for(int c = 0; c < cnt; c++){
        ppow *= p_cur;
        for(int i = 0; i < r_prev; i++){
          res[r++] = res[i] * ppow;
        }
      }
      r_prev = r;
    }
    return res;
  }
  int mobius_function(long long x){
    auto P = rho_factorization::factorize(x);
    for(long long p : P) if((x / p) % p == 0) return 0;
    return P.size() % 2 == 0 ? 1 : -1;
  }
  unsigned long long totient_function(unsigned long long n){
    unsigned long long res = n;
    auto prims = rho_factorization::prime_factor(n);
    for(auto p : prims) res -= res / p;
    return res;
  }
};
// p: 素数
unsigned long long find_primitive_root(unsigned long long p){
  static std::random_device seed_gen;
  static std::mt19937_64 engine(seed_gen());
  //assert(miller_rabin_mr(p));
  auto primes = rho_factorization::prime_factor(p - 1);
  while(true){
    bool f = true;
    unsigned long long a = engine() % (p - 1) + 1;
    for(unsigned long long pk : primes){
      // a ^ (p - 1) / pk ≡ 1 (mod p) -> no
      if(mod_pow_mr(a, (p - 1) / pk, p) == 1){
        f = false;
        break;
      }
    }
    if(f) return a;
  }
}

struct bigint_prime_factor{
  int sq, max_x;
  std::vector<int> low;
  std::unordered_map<int, int> high;
  bigint_prime_factor(int max_x): sq(sqrt(max_x)), max_x(max_x), low(sq + 1, 0){
    // 篩が初期化済か
    assert(prime_sieve::min_factor.size() > max_x);
  }
  // O(log(x))
  // x^kを掛ける
  void mul(int x, int k = 1){
    assert(0 < x && x <= max_x);
    while(x > 1){
      int p = prime_sieve::min_factor[x];
      if(p <= sq) low[p] += k;
      else{
        auto [itr, f] = high.emplace(x, k);
        if(!f) itr->second += k;
      }
      x /= p;
    }
  }
  // O(log(x))
  // x^kで割る(割り切れない場合エラー)
  void div(int x, int k = 1){
    assert(0 < x && x <= max_x);
    while(x > 1){
      int p = prime_sieve::min_factor[x];
      if(p <= sq){
        assert(low[p] >= k);
        low[p] -= k;
      }else{
        auto itr = high.find(p);
        assert(itr != high.end() && itr->second >= k);
        itr->second -= k;
        if(itr->second == 0) high.erase(itr);
      }
      x /= p;
    }
  }
  // 素数pで何回割れるか, O(1)
  // pが素数でない場合0
  int k_divisible_prime(int p){
    if(p > max_x || !prime_sieve::is_prime(p)) return 0;
    if(p > sq){
      auto itr = high.find(p);
      return itr == high.end() ? 0 : itr->second;
    }
    return low[p];
  }
  // xで何回割れるか, O(log(x))
  // x = 1のときinf回割れる
  // @param 0 < x <= 篩の最大値
  int k_divisible(int x){
    assert(x > 0);
    static constexpr int inf = std::numeric_limits<int>::max();
    int ans = inf;
    for(auto [p, num] : prime_sieve::factorize(x)){
      if(p > sq){
        auto itr = high.find(p);
        if(itr == high.end()) return 0;
        ans = std::min(ans, itr->second / num);
      }else{
        ans = std::min(ans, low[p] / num);
      }
    }
    return ans;
  }
  // rで何回割り切れるか, O(sqrt(r.max_x)以下の素数 + rの素因数の種類)
  int k_divisible(const bigint_prime_factor &r){
    static constexpr int inf = std::numeric_limits<int>::max();
    int ans = inf;
    
    for(int i = 0; prime_sieve::primes[i] <= r.sq; i++){
      int p = prime_sieve::primes[i];
      if(!r.low[p]) continue;
      if(p <= sq){
        if(low[p] < r.low[p]) return 0;
        ans = std::min(ans, low[p] / r.low[p]);
      }else{
        auto itr = high.find(p);
        if(itr->second < r.low[p]) return 0;
        ans = std::min(ans, itr->second / r.low[p]);
      }
    }
    for(auto [p, num] : r.high){
      assert(num);
      if(p <= sq){
        if(low[p] < num) return 0;
        ans = std::min(ans, low[p] / num);
      }else{
        auto itr = high.find(p);
        if(itr->second < num) return 0;
        ans = std::min(ans, itr->second / num);
      }
    }
    return ans;
  }
};
// res[i] := v[i]の位数
template<typename mint>
std::vector<long long> multiplicative_order_many(const std::vector<mint> &v){
  int n = v.size();
  std::vector<long long> res(n, 1);
  for(auto [p, q] : rho_factorization::factorize2(mint::mod() - 1)){
    long long x = mint(p).pow(q).val(), y = (mint::mod() - 1) / x;
    for(int i = 0; i < n; i++){
      long long z = x;
      for(int j = 0; j < q; j++){
        z /= p;
        if(v[i].pow(y * z).val() != 1){
          res[i] *= z * p;
          break;
        }
      }
    }
  }
  return res;
}
#endif