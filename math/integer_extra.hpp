#ifndef _INTEGER_EXTRA_H_
#define _INTEGER_EXTRA_H_
#include "integer.hpp"
#include "prime.hpp"
#include <unordered_map>
// m_iで割るとr_i余る情報から(全ての条件を満たすlcm(m_i)以下の数)%modを返す
// mが互いに素でないと壊れる
// O(n^2)
long long garner_coprime(std::vector<long long> r, std::vector<long long> m, long long mod){
  int n = r.size();
  assert(n == m.size());

  // ans = r_1 + x_1 * m_1 + x_2 * m_1 * m_2 ...
  // x_i = (r_i - (i - 1項までの和)) / (m_1 * m_2 ... m_i-1)  mod m_i
  // M_i := (m_1 * m_2 ... m_i - 1) mod m_i
  std::vector<long long> accum(n + 1, 0), M(n + 1, 1);
  m.push_back(mod);

  for(int i = 0; i < n; i++){
    auto [g, inv] = inv_ext_gcd(M[i], m[i]);
    assert(g == 1);
    long long x = (((__int128_t)r[i] - accum[i]) * inv) % m[i];
    if(x < 0) x += m[i];
    for(int j = i + 1; j <= n; j++){
      accum[j] = (accum[j] + (__int128_t)M[j] * x) % m[j];
      M[j] = ((__int128_t)M[j] * m[i]) % m[j];
    }
  }
  return accum[n];
}
// m_iで割るとr_i余る情報から(全ての条件を満たすlcm(m_i)以下の数)%modを返す
// 解が無い場合は-1
// O(n^2 + ∑m_i^(1/4)log(m_i))
long long garner(std::vector<long long> r, std::vector<long long> m, long long mod){
  int n = r.size();
  assert(n == m.size());
  std::unordered_map<long long, std::pair<long long, int>> P;
  for(int i = 0; i < n; i++){
    auto f = rho_factorization::factorize(m[i]);
    std::sort(f.begin(), f.end());
    f.push_back(0);
    long long p = -1, p_pow = 1;
    for(int j = 0; j < f.size(); j++){
      if(f[j] == p){
        p_pow *= p;
      }else{
        if(p != -1){
          auto itr = P.find(p);
          if(itr == P.end()){
            P.emplace(p, std::make_pair(p_pow, i));
          }else{
            auto [p_pow_max, idx_max] = itr->second;
            long long p_pow_gcd = std::min(p_pow_max, p_pow);
            if(r[idx_max] % p_pow_gcd != r[i] % p_pow_gcd) return -1; // 解なし, 例: 2で割って1余りつつ0余る数
            if(p_pow_max < p_pow){
              m[idx_max] /= p_pow_max;
              r[idx_max] %= m[idx_max];
              itr->second = std::make_pair(p_pow, i);
            }else{
              m[i] /= p_pow;
              r[i] %= m[i];
            }
          }
        }
        p = f[j];
        p_pow = p;
      }
    }
  }
  /*
  // rが全て0の場合答えは0になるが, 正整数が欲しい場合はlcm(m)が答え
  bool is_all_zero = true;
  long long _lcm = 1 % mod;
  for(int i = 0; i < n; i++){
    if(r[i]){
      is_all_zero = false;
      break;
    }
    _lcm = (__int128_t)_lcm * m[i] % mod;
  }
  if(is_all_zero) return _lcm;
  */
  return garner_coprime(r, m, mod);
}

// a^k >= xとなる最小のa^k
uint64_t ceil_pow(uint64_t x, uint64_t a){
  assert(a > 1);
  if(x == 0) return 1;
  static uint64_t INF = std::numeric_limits<uint64_t>::max();
  if(a == 2){
    uint64_t ret = uint64_t(1) << (63 - __builtin_clzll(x));
    if(ret == x) return x;
    assert(ret <= (INF >> 1));
    return ret << 1;
  }
  if(a > 10){
    uint64_t ret = 1;
    while(true){
      if(ret > x / a) break;
      ret *= a;
    }
    if(ret == x) return ret;
    assert(ret <= INF / a);
    return ret * a;
  }
  static std::vector<std::vector<uint64_t>> pow_table(11);
  if(pow_table[a].empty()){
    uint64_t tmp = 1;
    pow_table[a].push_back(1);
    while(true){
      if(tmp > INF / a) break;
      pow_table[a].push_back(tmp *= a);
    }
  }
  auto itr = std::lower_bound(pow_table[a].begin(), pow_table[a].end(), x);
  assert(itr != pow_table[a].end());
  return *itr;
}
// a^k <= xを満たす最大のa^k
uint64_t floor_pow(uint64_t x, uint64_t a){
  assert(x > 0 && a > 1);
  if(a == 2) return uint64_t(1) << (63 - __builtin_clzll(x));
  if(a > 10){
    uint64_t ret = 1;
    while(true){
      if(ret > x / a) return ret;
      ret *= a;
    }
  }
  static std::vector<std::vector<uint64_t>> pow_table(11);
  static uint64_t INF = std::numeric_limits<uint64_t>::max();
  if(pow_table[a].empty()){
    uint64_t tmp = 1;
    pow_table[a].push_back(1);
    while(true){
      if(tmp > INF / a) break;
      pow_table[a].push_back(tmp *= a);
    }
  }
  return *--std::upper_bound(pow_table[a].begin(), pow_table[a].end(), x);
}

// floor(x / x以下の自然数)としてあり得る値の昇順
template<typename T>
std::vector<T> enumerate_quotients(T x){
  assert(x >= 0);
  if(x == 0) return {};
  T sq = sqrtl(x);
  std::vector<T> ans(sq);
  std::iota(ans.begin(), ans.end(), 1);
  if(x / sq != ans.back()) ans.push_back(x / sq);
  for(T i = sq - 1; i >= 1; i--) ans.push_back(x / i);
  return ans;
}
// floor(x / y(x以下の自然数))としてあり得る値の昇順
// {商, yの最小値, yの最大値 + 1}
template<typename T>
std::vector<std::tuple<T, T, T>> enumerate_quotients3(T x){
  assert(x >= 0);
  if(x == 0) return {};
  T sq = sqrtl(x);
  std::vector<std::tuple<T, T, T>> ans(sq);
  for(T i = 1, pre1 = x; i <= sq; i++){
    T pre2 = x / (i + 1);
    ans[i - 1] = {i, pre2 + 1, pre1 + 1};
    std::swap(pre1, pre2);
  }
  if(x / sq == std::get<0>(ans.back())) sq--;
  for(T i = sq; i >= 1; i--){
    ans.push_back({x / i, i, i + 1});
  }
  return ans;
}
// enumerate_quotients(x)を呼んだ時にyは何番目の要素に含まれるか
// 含まれない場合は-1
template<typename T>
int index_quotients(T x, T y){
  if(y <= 0 || y > x) return -1;
  T z = x / y;
  if(z <= y) return z;
  T sq = sqrtl(x);
  if(x / sq == sq){
    return 2 * sq - 1 - y; // 2sq - 1型
  }else{
    return 2 * sq - y; // 2sq型
  }
}
// enumerate_quotients(x)を呼んだ時のk番目の要素
// ない場合は-1
template<typename T>
T kth_quotients(T x, int k){
  T sq = sqrtl(x);
  if(k < sq) return k + 1;
  int len = 2 * sq - (x / sq == sq);
  if(len <= k) return -1;
  return x / (len - k);
}
// enumerate_quotients3(x)を呼んだ時のk番目の要素
// ない場合は-1
template<typename T>
std::tuple<T, T, T> kth_quotients3(T x, int k){
  T sq = sqrtl(x);
  if(k < sq) return {k + 1, x / (k + 2) + 1, x / (k + 1) + 1};
  int len = 2 * sq - (x / sq == sq);
  if(len <= k) return {-1, -1, -1};
  len -= k;
  return {x / len, len, len + 1};
}
// 2変数enumerate_quotients
  // (x / k, y / k)ごとに{lk, rk}を列挙
  // {x / k, y / k} = {0, 0}となるようなものは列挙しない(k > max(x, y))
  template<typename T>
  std::vector<std::tuple<T, T>> enumerate_quotients_pair(T x, T y){
    if(x < y) std::swap(x, y);
    auto X = enumerate_quotients3(x);
    auto Y = enumerate_quotients3(y);
    X.push_back({-1, 0, 0});
    Y.push_back({-1, 0, 0});
    if(x != y) Y.insert(Y.begin(), {-1, y + 1, -1});
    int r = x + 1;
    std::vector<std::tuple<T, T>> res;
    int i = 0, j = 0;
    while(i + 1 < X.size() || j + 1 < Y.size()){
      auto [_, xl, xr] = X[i];
      auto [__, yl, yr] = Y[j];
      if(xl > yl){
        res.push_back({xl, r});
        r = xl;
        i++;
      }else if(xl < yl){
        res.push_back({yl, r});
        r = yl;
        j++;
      }else{
        res.push_back({xl, r});
        r = xl;
        i++, j++;
      }
    }
    return res;
  }

// x%1 + x%2 ... x%n
// @param 0 <= x, n
// O(√x)
long long mod_sum(long long x, long long n){
  assert(0 <= x && 0 <= n);
  if(x == 0) return 0;
  long long ans = x * n;
  if(n > x) n = x;
  // floor(x / i)の値でグループ分け
  long long sq = sqrtl(x);
  for(int i = 1; i <= sq; i++){
    long long l = x / (i + 1) + 1, r = std::min(n + 1, x / i + 1);
    if(l < r){
      long long sum_lr = (r * (r - 1) - l * (l - 1)) / 2;
      ans -= i * sum_lr;
    }
  }
  if(x / sq == sq) sq--;
  for(int i = std::min(n, sq); i >= 1; i--){
    ans -= i * (x / i);
  }
  return ans;
}

// a ^ k <= xを満たす最大のa
uint64_t kth_root_stable(uint64_t x, uint64_t k){
  if(!x) return 0;
  if(k == 1 || x == 1) return x;
  if(k >= 64) return 1;
  uint64_t l = 1, r = x;
  const static uint64_t threshold = std::numeric_limits<uint64_t>::max();
  while(r - l > 1){
    bool f = false;
    uint64_t mid = l + ((r - l) >> 1), z = 1;
    uint64_t lim = threshold / mid;
    for(int i = 0; i < k; i++){
      if(z > lim){
        f = true;
        break;
      }
      z *= mid;
    }
    if(f | (z > x)) r = mid;
    else l = mid;
  }
  return l;
}

// a^k <= x を満たす最大のa
uint64_t kth_root_fast(uint64_t x, uint64_t k){
  if(x <= 1) return (!k ? 1 : x);
  if(k <= 2) return (k <= 1 ? !k ? 1 : x : sqrtl(x));
  if(k >= 64||!(x >> k)) return 1;
  const static int sz[8] = {40, 31, 27, 24, 22, 21, 20, 19};
  static std::vector<std::vector<uint64_t>> table;
  if(table.empty()){
    table.resize(40);
    for(int i = 0; i < 40; i++){
      for(uint64_t j = 0; j < 8 && sz[j] > i; j++){
        table[i].push_back((!i ? 1 : table[i - 1][j]) *(j + 3));
      }
    }
  }
  if(k >= 19){
    int ans = 10;
    for(;ans > 2; ans--){
      if(sz[ans - 3]<k || table[k - 1][ans - 3] > x) continue;
      return ans;
    }
    return 2;
  }
  uint64_t r = (k != 3 ? pow(x, (1.0 + 1e-6) / k) : cbrt(x) + 1);
  const static uint64_t threshold = std::numeric_limits<uint64_t>::max();
  while(true){
    uint64_t lim = threshold / r, z = 1;
    for(int i = 0; i < k; i++, z *= r) if(z > lim) goto upper;
    if(z > x) upper: r--;
    else break;
  }
  return r;
}


std::vector<long long> v{1,2,4,6,12,24,36,48,60,120,180,240,360,720,840,1260,1680,2520,5040,7560,10080,15120,20160,25200,27720,45360,50400,55440,83160,110880,166320,221760,277200,332640,498960,554400,665280,720720,1081080,1441440,2162160,2882880,3603600,4324320,6486480,7207200,8648640,10810800,14414400,17297280,21621600,32432400,36756720,43243200,61261200,73513440,110270160,122522400,147026880,183783600,245044800,294053760,367567200,551350800,698377680,735134400,1102701600,1396755360,2095133040,2205403200,2327925600,2793510720,3491888400,4655851200,5587021440,6983776800,10475665200,13967553600,20951330400,27935107200,41902660800,48886437600,64250746560,73329656400,80313433200,97772875200,128501493120,146659312800,160626866400,240940299600,293318625600,321253732800,481880599200,642507465600,963761198400,1124388064800,1606268664000,1686582097200,1927522396800,2248776129600,3212537328000,3373164194400,4497552259200,6746328388800,8995104518400,9316358251200,13492656777600,18632716502400,26985313555200,27949074753600,32607253879200,46581791256000,48910880818800,55898149507200,65214507758400,93163582512000,97821761637600,130429015516800,195643523275200,260858031033600,288807105787200,391287046550400,577614211574400,782574093100800,866421317361600,1010824870255200,1444035528936000,1516237305382800,1732842634723200,2021649740510400,2888071057872000,3032474610765600,4043299481020800,6064949221531200,8086598962041600,10108248702552000,1212898443062400,18194847664593600,20216497405104000,24259796886124800,30324746107656000,36389695329187200,48519593772249600,60649492215312000,72779390658374400,74801040398884800,106858629141264000,112201560598327200,149602080797769600,224403121196654400,299204161595539200,374005201994424000,448806242393308800,673209363589963200,748010403988848000,897612484786617600};

std::vector<long long> ans{1,2,3,4,6,8,9,10,12,16,18,20,24,30,32,36,40,48,60,64,72,80,84,90,96,100,108,120,128,144,160,168,180,192,200,216,224,240,256,288,320,336,360,384,400,432,448,480,504,512,576,600,640,672,720,768,800,864,896,960,1008,1024,1152,1200,1280,1344,1440,1536,1600,1680,1728,1792,1920,2016,2048,2304,2400,2688,2880,3072,3360,3456,3584,3600,3840,4032,4096,4320,4608,4800,5040,5376,5760,6144,6720,6912,7168,7200,7680,8064,8192,8640,9216,10080,10368,10752,11520,12288,12960,13440,13824,14336,14400,15360,16128,16384,17280,18432,20160,20736,21504,23040,24576,25920,26880,27648,28672,28800,30720,32256,32768,34560,36864,40320,41472,43008,46080,48384,49152,51840,53760,55296,57600,61440,62208,64512,65536,69120,73728,80640,82944,86016,92160,96768,98304,103680
};

// 高度合成数(N以下の数の中で最も多い約数を持つ数)
// {N以下の高度合成数, その約数}
std::pair<long long, long long> highly_composite_number(long long N){
  assert(0 < N && N <= 1000000000000000000);
  int idx = upper_bound(v.begin(), v.end(), N) - v.begin() - 1;
  assert(idx != 0);
  return {v[idx], ans[idx]};
}


//---tips---
//φ(n) := [1, n]の中でnと互いに素な正整数の数
//φ(n) =  n * π{p:prime factor of n}(1 - 1/p)
//aとnが互いに素な場合, a^φ(n) ≡ 1 (mod n)
unsigned long long totient_function(unsigned long long n){
  unsigned long long res = n;
  auto prims = rho_factorization::prime_factor(n);
  for(auto p : prims) res -= res / p;
  return res;
}
// n以下のtotient_sum
std::vector<unsigned long long> totient_sum_table(unsigned long long n){
  std::vector<unsigned long long> res(n + 1);
  std::iota(res.begin() + 1, res.end(), 0);
  res[1] = 1;
  for(int i = 2; i <= n; i++){
    // prime
    if(res[i] == i - 1){
      for(int j = 2 * i; j <= n; j += i){
        res[j] -= res[j] / i;
      }
    }
  }
  for(int i = 2; i <= n; i++) res[i] += res[i - 1];
  return res;
}
// [1]
// totient_sum(n) = (i, j){1 <= i <= n, 1 <= j <= i}となるペアの数 - そのうちgcd(i, j) != 1であるものの数
//                = n * (n + 1) / 2 - ∑{2 <= g <= n} gcd(i, j) == gのペア
//                = n * (n + 1) / 2 - ∑{2 <= g <= n} totient_sum(n / g)

// この演算をメモ化再帰で行った時の n / g (切り捨て)の種類数について考える

// [2]
// floor(x / ab) = floor(floor(x / a) / b)
// 証明:
// x = qab + r (0 <= r < ab)とすると
// 左辺 = floor(q + r / ab) = q
// 右辺 = floor((qb + r') / b) (0 <= r' < b) = q

// [3]
// [2]よりnを正整数で0回以上除算した時に得られる商はO(√n)種類
__uint128_t totient_sum(unsigned long long n){
  static std::vector<unsigned long long> low;
  if(low.empty()) low = totient_sum_table(std::min(n, (unsigned long long)4000000));
  std::unordered_map<unsigned long long, __uint128_t> mp;

  unsigned long long memo_max = 0;
  auto tsum = [&](auto &&tsum, unsigned long long m)->__uint128_t{
    if(m < low.size()) return low[m];
    if(m <= memo_max) return mp.at(m);
    __uint128_t res = (__uint128_t)m * (m + 1) / 2;
    auto d = enumerate_quotients3(m);
    for(auto [q, a, b] : d){
      if(q == m) continue;
      res -= (b - a) * tsum(tsum, q);
    }
    mp.emplace(m, res);
    memo_max = m;
    return res;
  };
  return tsum(tsum, n);
}

// res[i][j] = gcd(i, j) = 1であるペアの数　(0 < x <= i, 0 < y <= j), O(nm)
std::vector<std::vector<int>> coprime_sum_table(int n, int m){
  assert(n <= 3000 && m <= 3000);
  std::vector<std::vector<int>> res(n + 1, std::vector<int>(m + 1, 1));
  for(int i = 0; i <= n; i++) res[i][0] = 0;
  for(int i = 0; i <= m; i++) res[0][i] = 0;
  for(int k = 2; k <= std::min(n, m); k++){
    if(res[k][k] == 0) continue;
    for(int i = k; i <= n; i += k){
      for(int j = k; j <= m; j += k){
        res[i][j] = 0;
      }
    }
  }
  for(int i = 1; i <= n; i++) for(int j = 2; j <= m; j++) res[i][j] += res[i][j - 1];
  for(int j = 1; j <= m; j++) for(int i = 2; i <= n; i++) res[i][j] += res[i - 1][j];
  return res;
}
// n, mが10^9で900msくらい
long long coprime_sum(int n, int m){
  static std::vector<std::vector<int>> low;
  if(low.empty()) low = coprime_sum_table(2000, 2000);
  std::unordered_map<long long, long long> mp;
  long long memo_max = 0;
  if(n < m) std::swap(n, m);
  auto csum = [&](auto &&csum, int x, int y)->long long {
    if(x < low.size()) return low[x][y];
    if(y == 0) return 0;
    long long xy = (long long)x * y;
    if(xy <= memo_max) return mp.at(xy);
    long long res = xy;
    auto d = enumerate_quotients_pair(x, y);
    for(auto [a, b] : d){
      if(a == 1) continue;
      res -= (b - a) * csum(csum, x / a, y / a);
    }
    mp.emplace(xy, res);
    memo_max = xy;
    return res;
  };
  return csum(csum, n, m);
}

/*
// mが固定の場合, 定数の前計算やmontgomery reductionによる高速化ができる
// 解が無いと壊れる
struct fixed_mod_garner{
  std::vector<unsigned long long> m, iM, accum;
  std::vector<std::vector<unsigned long long>> M; // M[i][j] := (m_1 ~ m_j-1の積) % m_i
  std::vector<montgomery_reduction_64bit> mr;

  void set_mod(std::vector<unsigned long long> _m){
    int n = _m.size();
    std::unordered_map<long long, std::pair<long long, int>> P;
    for(int i = 0; i < n; i++){
      auto f = rho_factorization::factorize(_m[i]);
      std::sort(f.begin(), f.end());
      f.push_back(0);
      long long p = -1, p_pow = 1;
      for(int j = 0; j < f.size(); j++){
        if(f[j] == p){
          p_pow *= p;
        }else{
          if(p != -1){
            auto itr = P.find(p);
            if(itr == P.end()){
              P.emplace(p, std::make_pair(p_pow, i));
            }else{
              auto [p_pow_max, idx_max] = itr->second;
              if(p_pow_max < p_pow){
                _m[idx_max] /= p_pow_max;
                itr->second = std::make_pair(p_pow, i);
              }else{
                _m[i] /= p_pow;
              }
            }
          }
          p = f[j];
          p_pow = p;
        }
      }
    }
    m = _m;
    mr.resize(n + 1);
    for(int i = 0; i < n; i++) mr[i].set_mod(m[i]);
    M.resize(n + 1, std::vector<unsigned long long>(n + 1));
    for(int i = 0; i < n; i++){
      M[i][0] = mr[i].generate(1);
      for(int j = 1; j <= i; j++){
        M[i][j] = mr[i].mul(M[i][j - 1], mr[i].generate(m[j - 1]));
      }
    }
    iM.resize(n);
    for(int i = 0; i < n; i++) iM[i] = mr[i].generate(inv_gcd(mr[i].reduce(M[i][i]), m[i]).second);
    accum.resize(n + 1, 0);
  }
  unsigned long long query(std::vector<unsigned long long> r, unsigned long long mod){
    int n = r.size();
    assert(n == m.size());
    std::fill(accum.begin(), accum.end(), 0);
    mr[n].set_mod(mod);
    M[n][0] = mr[n].generate(1);
    for(int i = 1; i <= n; i++) M[n][i] = mr[n].mul(M[n][i - 1], mr[n].generate(m[i - 1]));
    for(int i = 0; i < n; i++){
      r[i] = mr[i].generate(r[i] % m[i]);
      unsigned long long x = mr[i].reduce(mr[i].mul(mr[i].sub(r[i], accum[i]), iM[i]));
      for(int j = i + 1; j <= n; j++) accum[j] = mr[j].add(accum[j], mr[j].mul(M[j][i], mr[j].generate(x)));
    }
    return mr[n].reduce(accum[n]);
  }
};
*/
// [0, N]を扱える 
struct Mobius{
  std::vector<int> is_prime, mobius;
  Mobius(int N): is_prime(N + 1, 1), mobius(N + 1, 1){
    for(int i = 2; i <= N; i++){
      if(!is_prime[i]) continue;
      mobius[i] = -1;
      for(int j = 2; j <= N / i; j++){
        is_prime[i * j] = 0;
        mobius[i * j] *= -1;
      }
      for(int j = 1; j <= N / ((long long)i * i); j++) mobius[j * i * i] = 0;
    }
  }
  int operator [](int k)const{return mobius[k];}
};
// O(√n * loglogn)
long long square_free_slow(long long n){
  int sq = sqrtl(n);
  auto m = Mobius(sq);
  long long ans = 0;
  for(int i = 1; i <= sq; i++) ans += m[i] * (n / (i * i));
  return ans;
}
// [1, n]の無平方数
long long square_free(long long n){
  long long I = 1;
  while(n > (I * I * I * I * I)) I++;
  int D = sqrtl(n / I);
  Mobius m(D);
  long long s1 = 0;
  for(int i = 1; i <= D; i++) s1 += m[i] * (n / ((long long)i * i));
  std::vector<int> M_list{0};
  int M = 0;
  for(int i = 1; i <= D; i++){
    M += m[i];
    M_list.push_back(M);
  }
  std::vector<int> Mxi_list;
  long long Mxi_sum = 0;
  for(int i = I - 1; i > 0; i--){
    long long Mxi = 1;
    long long xi = sqrtl(n / i);
    int sqd = sqrt(xi);
    for(int j = 1; j <= xi / (sqd + 1); j++){
      Mxi -= (xi / j - xi / (j + 1)) * M_list[j];
    }
    for(int j = 2; j <= sqd; j++){
      if(xi / j <= D){
        Mxi -= M_list[xi / j];
      }else{
        Mxi -= Mxi_list[I - j * j * i - 1];
      }
    }
    Mxi_list.push_back(Mxi);
    Mxi_sum += Mxi;
  }
  long long s2 = Mxi_sum - (I - 1) * M_list.back();
  return s1 + s2;
}

#endif