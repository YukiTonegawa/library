#ifndef _MOD_EXTRA_H_
#define _MOD_EXTRA_H_
#include "mod.hpp"
#include "../data_structure/basic/hash_table.hpp"

// 64bit modint


// mod == 2: 定数時間
// modが素数: O(min(n, mod) + log(n))
template<int id>
struct lucas_prime{
  using mint = dynamic_modint<id>;
  modcomb<mint> mcb;
  void set_mod(int mod){
    mint::set_mod(mod);
  }
  int comb(long long n, long long r){
    if(mint::mod() == 1 || n < 0 || r < 0 || n < r) return 0;
    if(mint::mod() == 2) return (n & r) == r;
    mcb.recalc(std::min(n, (long long)mint::mod()));
    mint res = 1;
    while(n){
      int x = n % mint::mod(), y = r % mint::mod();
      res *= mcb.comb(x, y);
      n /= mint::mod(), r /= mint::mod();
    }
    return res.val();
  }
};
#include <random>
// modpow(x,2,mod) == aとなるxを返す
// 存在しないなら-1
// mod は素数
long long modsqrt(long long a, long long mod){
  a %= mod;
	if(a == 0) return 0LL;
	if(mod == 2) return 1LL;
	if(modpow(a, (mod - 1) / 2, mod) != 1) return -1LL;
	if(mod % 4 == 3) return modpow(a, mod / 4 + 1, mod);
	long long q = mod - 1, m = 0;
	while(q % 2 == 0) q >>= 1, m++;
	std::mt19937 mt;
	long long z;
	do{
		z = mt() % mod;
	}while(modpow(z, (mod - 1) / 2, mod) != mod - 1);
	long long c = modpow(z, q, mod);
	long long t = modpow(a, q, mod);
	long long r = modpow(a, (q + 1) >> 1, mod);
	for(; m > 1; --m) {
	  long long tmp = modpow(t, 1LL << (m - 2), mod);
		if(tmp != 1) r = r * c % mod, t = t * (c * c % mod) % mod;
		c = c * c % mod;
	}
	return r;
}

//n次以下の多項式に対し
//f(0) ~ f(n)を与えf(p)を求める
//O(n log(MOD))
template<typename mint>
mint __lagrange(const std::vector<mint> &y, mint p, modcomb<mint> &mcb){
  int sz = y.size();
  mcb.recalc(sz);
  mint M = 1, res = 0;
  std::vector<mint> itable(sz, 1), num(sz);
  if(p.val() < sz) return y[p.val()];
  for(int i = 0; i < sz; i++){
    M *= p - i;
    num[i] = p - i;
  }
  uint32_t cnt = mint::mod() - 2;
  while(cnt){
    if(cnt & 1){
      for(int i = 0; i < sz; i++) itable[i] *= num[i], num[i] *= num[i];
    }else{
      for(int i = 0; i < sz; i++) num[i] *= num[i];
    }
    cnt >>= 1;
  }
  for(int i = 0; i < sz; i++){
    mint iQ = mcb.finv(i) * mcb.finv(sz - 1 - i);
    if((sz - i - 1) & 1) iQ *= -1;
    res += y[i] * iQ * itable[i];
  }
  return res * M;
}
template<typename mint>
mint lagrange(const std::vector<mint> &y, mint p){
  modcomb<mint> mcb;
  return __lagrange(y, p, mcb);
}

template<typename mint>
mint riid(mint r, int d, long long n){
  if(n == 0) return 0;
  if(r.val() == 0){
    return d == 0 ? mint(1) : 0;
  }
  n--;
  std::vector<mint> y(d + 2), ipow(d + 2, 1), tbl(d + 2);
  for(int i = 0; i < d + 2; i++) tbl[i] = i;
  int cnt = d;
  while(cnt){
    if(cnt & 1){
      for(int i = 0; i < d + 2; i++) ipow[i] *= tbl[i], tbl[i] *= tbl[i];
    }else{
      for(int i = 0; i < d + 2; i++) tbl[i] *= tbl[i];
    }
    cnt >>= 1;
  }
  mint tmp = 0, rpow = 1, last = r.pow(n % (mint::mod() - 1));
  n %= mint::mod();
  for(int i = 0; i < d + 2; i++){
    tmp += rpow * ipow[i];
    rpow *= r;
    y[i] = tmp;
  }
  modcomb<mint> mcb(y.size());

  if(r.val() == 1) return __lagrange<mint>(y, n, mcb);
  ipow[0] = 1, ipow[1] = -r;
  mint comb = 1, c = 0;
  for(int i = 2; i < d + 2; i++) ipow[i] = ipow[i - 1] * ipow[1];
  for(int i = 0; i < d + 1; i++){
    comb *= mcb.inv(i + 1) * mint(d + 1 - i);
    c += comb * ipow[d - i] * y[i];
  }
  mint di = mint(1 - r);
  c *= di.pow(d + 1).inv();
  mint powerRinv = 1, rinv = r.inv();
  for(int i = 0; i < d + 1; i++){
    y[i] -= c;
    y[i] *= powerRinv;
    powerRinv *= rinv;
  }
  y.pop_back();
  return c + last * __lagrange<mint>(y, n, mcb);
}

// x^k ≡ yであるような最小のk(存在しない場合-1)
// 0^0 = 1とする
// 素数modである必要はない
template<typename mint>
long long discrete_logarithm_bsgs(mint x, mint y){
  if(y.val() == 1) return 0;
  if(y.val() == 0){
    if(x.val() == 0) return mint::mod() == 1 ? 0 : 1;
    mint tmp = x;
    for(int i = 2; i <= 64; i++){
      tmp *= x;
      if(tmp.val() == 0) return i;
    }
    return -1;
  }
  long long sq = std::sqrt(mint::mod()) + 1;
  uhash_map<unsigned long long, int> table;
  mint xby = y;
  for(int b = 0; b < sq; b++, xby *= x) table.emplace_replace(xby.val(), b);
  mint xH = x.pow(sq);
  mint xaH = xH;
  for(int a = 1; a <= sq; a++, xaH *= xH){
    auto [f, id] = table.at(xaH.val());
    if(!f) continue;
    long long res = (long long)a * sq - id;
    if(x.pow(res) != y) return -1;
    return res;
  }
  return -1;
}
#endif