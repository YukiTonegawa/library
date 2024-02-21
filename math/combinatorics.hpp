#ifndef _COMBINATORICS_H_
#define _COMBINATORICS_H_
#include "mod.hpp"
#include "fps_extra.hpp"
#include "../misc/random_number.hpp"

template<typename mint>
mint factorial_mod_large(long long n){
  if(n >= mint::mod()) return 0;
  if(!n) return 1;

  using fps = formal_power_series<mint>;
  int block_size = round(sqrt(mint::mod())); // たくさんクエリがくる場合は /= 4する
  int block_n = mint::mod() / block_size;
  static std::vector<mint> x, fx;

  // init O(k * log^2 k), k = max(block_size, mod / block_size)
  // O(√mod * log^2 mod) if block_size = block_n = √mod
  if(fx.empty()){
    x.resize(block_n);
    for(int i = 0; i < block_n; i++) x[i] = (i + 1) * block_size;
    std::vector<fps> tmp(block_size);
    for(int i = 0; i < block_size; i++) tmp[i] = {-i, 1};
    fps f = multiply_constant_degree<mint>(tmp); // x(x - 1)...(x - block_size + 1)
    fx = multipoint_evaluation<mint>(f, x);
    for(int i = 1; i < block_n; i++) fx[i] *= fx[i - 1];
  }
  // O(block_size)
  auto restore_fact = [&](long long a)->mint{
    int b = a / block_size;
    mint res = (!b ? 1 : fx[b - 1]);
    for(long long i = (long long)b * block_size + 1; i <= a; i++) res *= i;
    return res;
  };
  return restore_fact(n);
}
template<typename mint>
mint combination_mod_large(long long n, long long k){
  if(n < 0 || k < 0 || n < k) return 0;
  if(n == k || k == 0) return 1;
  // [1]
  // comb(n, k) = (n * n - 1 .... n - k + 1) / k !
  // n - k + 1 > k となるように変形(nCk = nC(n - k))したとき, 分子がmod上で0を跨いでいるなら0
  k = std::min(k, n - k);
  long long numerator_min = n - k + 1;
  if(numerator_min % mint::mod() == 0 || (n / mint::mod() != numerator_min / mint::mod())) return 0;

  // [1]より1 <= k < mod かつ nPk ≡ (n%mod)Pk
  assert(1 <= k && k < mint::mod());
  n %= mint::mod();

  return factorial_mod_large<mint>(n) * (factorial_mod_large<mint>(k) * factorial_mod_large<mint>(n - k)).inv();
}

template<typename mint>
mint permutation_mod_large(long long n, long long k){
  if(n < 0 || k < 0 || n < k) return 0;
  if(k == 0) return 1;
  // [1]
  // perm(n, k) = n * n - 1 .... n - k + 1
  // mod上で0を跨いでいるなら0
  long long numerator_min = n - k + 1;
  if(numerator_min % mint::mod() == 0 || (n / mint::mod() != numerator_min / mint::mod())) return 0;

  // [1]より1 <= k < mod かつ nPk ≡ (n%mod)Pk
  assert(1 <= k && k < mint::mod());
  n %= mint::mod();

  return factorial_mod_large<mint>(n) / factorial_mod_large<mint>(n - k);
}
template<typename mint>
formal_power_series<mint> stirling_number_first_kind(int n){
  if(!n) return {1};
  std::vector<formal_power_series<mint>> f(n);
  for(int i = 0; i < n; i++) f[i] = {-i, 1};
  return multiply_constant_degree<mint>(f);
}

// max_length: 文字数制限
// max_lengthに収まるように埋め込む
template<typename mint>
void __umekomi_factorial_mod(int max_length){
  int w = size(std::to_string(mint::mod())) + 1; // 1要素 + 区切り文字のサイズ
  int n = std::max(2, max_length / w), m = mint::mod() / n;
  std::cout << "int umekomi_fact_size = " << n + 1 << ';' << '\n';
  std::cout << "int umekomi_fact_interval = " << m << ';' << '\n';
  mint z = 1 % mint::mod();
  std::cout << "std::vector<int> umekomi_fact{" << z << ',';
  int next = m;
  for(int i = 1; i < mint::mod(); i++){
    z *= i;
    if(i == next){
      std::cout << z << ',';
      next += m;
    }
  }
  std::cout << "};" << '\n';
}

int umekomi_fact_size;
int umekomi_fact_interval;
std::vector<int> umekomi_fact;

template<typename mint>
mint umekomi_factorial_mod(long long n){
  if(n >= mint::mod()) return 0;
  int b = n / umekomi_fact_interval;
  assert(b < umekomi_fact_size);
  mint x = umekomi_fact[b];
  for(int i = b * umekomi_fact_interval + 1; i <= n; i++) x *= i;
  return x;
}

template<typename mint>
mint umekomi_combination_mod(long long n, long long k){
  if(n < 0 || k < 0 || n < k) return 0;
  if(n == k || k == 0) return 1;
  // [1]
  // comb(n, k) = (n * n - 1 .... n - k + 1) / k !
  // n - k + 1 > k となるように変形(nCk = nC(n - k))したとき, 分子がmod上で0を跨いでいるなら0
  k = std::min(k, n - k);
  long long numerator_min = n - k + 1;
  if(numerator_min % mint::mod() == 0 || (n / mint::mod() != numerator_min / mint::mod())) return 0;
  // [1]より1 <= k < mod かつ nPk ≡ (n%mod)Pk
  assert(1 <= k && k < mint::mod());
  n %= mint::mod();
  return umekomi_factorial_mod<mint>(n) * (umekomi_factorial_mod<mint>(k) * umekomi_factorial_mod<mint>(n - k)).inv();
}
template<typename mint>
mint umekomi_permutation_mod(long long n, long long k){
  if(n < 0 || k < 0 || n < k) return 0;
  if(k == 0) return 1;
  // [1]
  // perm(n, k) = n * n - 1 .... n - k + 1
  // mod上で0を跨いでいるなら0
  long long numerator_min = n - k + 1;
  if(numerator_min % mint::mod() == 0 || (n / mint::mod() != numerator_min / mint::mod())) return 0;
  // [1]より1 <= k < mod かつ nPk ≡ (n%mod)Pk
  assert(1 <= k && k < mint::mod());
  n %= mint::mod();
  return umekomi_factorial_mod<mint>(n) / umekomi_factorial_mod<mint>(n - k);
}

struct four_square_theorem{
  int n;
  std::vector<int> A;
  std::vector<int> sz;
  // 0 <= i < Nについて以下のテーブルを作る
  // sum(A) = iを満たしlen(A)が最小(常に4つ以下になる)
  // O(N√N)
  // N = 10^6で1sec程度
  void init(int _n){
    n = _n;
    A.resize(n, -1);
    sz.resize(n, 5);
    sz[0] = 0;
    for(int i = 1; i < n; i++){
      for(int j = 1; j * j <= i; j++){
        if(sz[i - j * j] + 1 < sz[i]){
          A[i] = j;
          sz[i] = sz[i - j * j] + 1;
        }
      }
    }
  }
  four_square_theorem(int _n){
    init(_n);
  }
  // x = ∑A[i]^2 を満たす要素数最小のA(高々4つ)
  std::vector<int> get_square(int x){
    assert(n);
    assert(0 <= x && x < n);
    std::vector<int> ret;
    while(x){
      int y = A[x];
      ret.push_back(y);
      x -= y * y;
    }
    return ret;
  }
  // x = ∑A[i]^2 を満たすA(要素数最小とは限らない, N = 10^5, x < 10^18なら高々6要素になる)
  // @param x <= 10^18
  std::vector<unsigned long long> get_square_large(unsigned long long x){
    std::vector<unsigned long long> ret;
    while(n <= x){
      unsigned long long y = sqrtl(x);
      ret.push_back(y);
      x -= y * y;
    }
    while(x){
      int y = A[x];
      ret.push_back(y);
      x -= y * y;
    }
    return ret;
  }
  // x = ∑ A[i] * (A[i] + 1) / 2 を満たすA(3つ以下)
  std::vector<int> get_trianglar(int x){
    auto v = get_square(8 * x + 3);
    assert(v.size() == 3);
    std::vector<int> ret;
    for(int i = 0; i < 3; i++) if(v[i] != 1) ret.push_back(v[i] / 2);
    return ret;
  }
  // x = ∑ A[i] * (A[i] + 1) / 2 を満たすA(要素数最小とは限らない, N = 10^5, x < 10^18なら高々5要素になる)
  // @param x <= 10^18
  std::vector<unsigned long long> get_trianglar_large(unsigned long long x){
    std::vector<unsigned long long> ret;
    while(n <= 8 * x + 3){
      unsigned long long y = sqrtl(2 * x);
      while(y * (y + 1) > 2 * x) y--;
      ret.push_back(y);
      x -= (y & 1 ? y * ((y >> 1) + 1) : (y >> 1) * (y + 1));
    }
    auto v = get_square(8 * x + 3);
    for(int i = 0; i < 3; i++) if(v[i] != 1) ret.push_back(v[i] / 2);
    return ret;
  }
};

// ピタゴラス数について
// a^2 + b^2 = c^3を満たす自然数の組(a, b, c)をピタゴラス数と呼ぶ
// bが偶数, a, cが奇数でなければならない
// 自然数x, y(x > y)を用いて(a, b, c) = (x^2 - y^2, 2xy, x^2 + y^2)を表すことができる　


// A[i] := c1 * nC1 + c2 * nC2 + .... cn * nCn
// cが与えられた時, 上の式を満たすAを返す
// cが全部1なら二項係数の式変形で解ける

/*
int main(){
  io_init();
  int n;
  std::cin >> n;
  modcomb<mint> mcb(n);
  simple_tree t(n);
  range(i, 0, n - 1){
    int a, b;
    std::cin >> a >> b;
    a--, b--;
    t.add_edge(a, {a, b});
    t.add_edge(b, {b, a});
  }

  vector<mint> c(n);
  range(i, 0, n){
    c[n - 1]++;
    for(auto e : t[i]){
      int sz = (e.t == t.par(i) ? n - t.size(i) : t.size(e.t));
      c[sz - 1]--;
    }
  }
  // T(n) = 2 * T(n / 2) + O(nlog(n))
  // ci * (1 + x) ^ i

  int N = 1 << bit_highest(n);
  if(N != n) N <<= 1;
  vector<fps> p2(19);
  p2[0] = {1, 1};
  range(i, 1, 19) p2[i] = p2[i - 1] * p2[i - 1];

  auto f = [&](auto &&f, int l, int r) -> fps {
    if(r - l == 1){
      if(l >= n) return fps{0};
      return fps{c[l], c[l]};
    }
    int mid = (l + r) / 2;
    auto L = f(f, l, mid);
    auto R = f(f, mid, r);
    int h = bit_highest(mid - l);
    assert((1 << h) == mid - l);
    R = R * p2[h] + L;
    return R;
  };
  fps ans = f(f, 0, N);
  range(i, 1, n + 1) std::cout << (ans.size() <= i ? 0 : ans[i]) << '\n';
}
*/

/*
template<typename mint>
struct dot_sum_perm{
  std::vector<mint> sum;
  modcomb<mint> mcb;

  void push_back(mint x){

  }
  void pop_back(){

  }
  void query(int r){
    //
  }
};
*/
// res[i] := ∑{0 <= j <= i} (i + s)Cj * v[j]
template<typename mint>
std::vector<mint> dot_sum_comb(std::vector<mint> v, int s){
  assert(s >= 0);
  int n = v.size();
  modcomb<mint> mcb(s + n);
  for(int i = 0; i < n; i++) v[i] *= mcb.finv(i);
  std::vector<mint> f(n);
  for(int i = 0; i < n; i++) f[i] = mcb.finv(s + i);
  f = convolution_mod<mint>(v, f);
  for(int i = 0; i < n; i++) f[i] *= mcb.fac(i + s);
  return std::vector<mint>(f.begin(), f.begin() + n);
}

// 
/*
struct online_dot_product_comb{
  
};
*/
#endif