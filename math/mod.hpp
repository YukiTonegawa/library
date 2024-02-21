#ifndef _MOD_H_
#define _MOD_H_
#include <vector>
#include <cassert>
#include <numeric>
#include <type_traits>
#include <iostream>
#include <ostream>
#include "minior/mod_base.hpp"

template<int m>
long long modpow(long long a, long long b){
  assert(0 <= b);
  assert(0 < m);
  a = safe_mod(a, m);
  long long ret = 1;
  while(b){
    if(b & 1) ret = (ret * a) % m;
    a = (a * a) % m;
    b >>= 1;
  }
  return ret;
}
// @param 0 <= b, 0 < m
long long modpow(long long a, long long b, int m){
  assert(0 <= b);
  assert(0 < m);
  a = safe_mod(a, m);
  long long ret = 1;
  while(b){
    if(b & 1) ret = (ret * a) % m;
    a = (a * a) % m;
    b >>= 1;
  }
  return ret;
}

struct modint_base {};
struct static_modint_base : modint_base {};

template <int m, std::enable_if_t<(1 <= m)>* = nullptr>
struct static_modint : static_modint_base{
  using mint = static_modint;
public:
  static constexpr int mod(){return m;}
  static mint raw(int v) {
    mint x;
    x._v = v;
    return x;
  }
  static_modint(): _v(0){}
  template <class T>
  static_modint(T v){
    long long x = v % (long long)umod();
    if (x < 0) x += umod();
    _v = x;
  }
  unsigned int val()const{return _v;}
  mint& operator++(){
    _v++;
    if (_v == umod()) _v = 0;
    return *this;
  }
  mint& operator--(){
    if (_v == 0) _v = umod();
    _v--;
    return *this;
  }
  mint operator++(int){
    mint result = *this;
    ++*this;
    return result;
  }
  mint operator--(int){
    mint result = *this;
    --*this;
    return result;
  }
  mint& operator+=(const mint& rhs){
    _v += rhs._v;
    if (_v >= umod()) _v -= umod();
    return *this;
  }
  mint& operator-=(const mint& rhs){
    _v -= rhs._v;
    if (_v >= umod()) _v += umod();
    return *this;
  }
  mint& operator*=(const mint& rhs){
    unsigned long long z = _v;
    z *= rhs._v;
    _v = (unsigned int)(z % umod());
    return *this;
  }
  mint& operator/=(const mint& rhs){return *this = *this * rhs.inv();}
  mint operator+()const{return *this;}
  mint operator-()const{return mint() - *this;}
  mint pow(long long n)const{
    assert(0 <= n);
    mint x = *this, r = 1;
    while(n){
      if (n & 1) r *= x;
      x *= x;
      n >>= 1;
    }
    return r;
  }
  mint inv()const{
    if(prime){
      assert(_v);
      return pow(umod() - 2);
    }else{
      auto eg = inv_gcd(_v, m);
      assert(eg.first == 1);
      return eg.second;
    }
  }
  friend mint operator+(const mint& lhs, const mint& rhs){return mint(lhs) += rhs;}
  friend mint operator-(const mint& lhs, const mint& rhs){return mint(lhs) -= rhs;}
  friend mint operator*(const mint& lhs, const mint& rhs){return mint(lhs) *= rhs;}
  friend mint operator/(const mint& lhs, const mint& rhs){return mint(lhs) /= rhs;}
  friend bool operator==(const mint& lhs, const mint& rhs){return lhs._v == rhs._v;}
  friend bool operator!=(const mint& lhs, const mint& rhs){return lhs._v != rhs._v;}
private:
  unsigned int _v;
  static constexpr unsigned int umod(){return m;}
  static constexpr bool prime = is_prime<m>;
};

template<int id> 
struct dynamic_modint : modint_base{
  using mint = dynamic_modint;
public:
  static int mod(){return (int)(bt.umod());}
  static void set_mod(int m){
    assert(1 <= m);
    bt = barrett(m);
  }
  static mint raw(int v){
    mint x;
    x._v = v;
    return x;
  }
  dynamic_modint(): _v(0){}
  template <class T>
  dynamic_modint(T v){
    long long x = v % (long long)(mod());
    if (x < 0) x += mod();
    _v = x;
  }
  unsigned int val()const{return _v;}
  mint& operator++(){
    _v++;
    if(_v == umod()) _v = 0;
    return *this;
  }
  mint& operator--(){
    if (_v == 0) _v = umod();
    _v--;
    return *this;
  }
  mint operator++(int){
    mint result = *this;
    ++*this;
    return result;
  }
  mint operator--(int){
    mint result = *this;
    --*this;
    return result;
  }
  mint& operator+=(const mint& rhs){
    _v += rhs._v;
    if(_v >= umod()) _v -= umod();
    return *this;
  }
  mint& operator-=(const mint& rhs){
    _v += mod() - rhs._v;
    if(_v >= umod()) _v -= umod();
    return *this;
  }
  mint& operator*=(const mint& rhs){
    _v = bt.mul(_v, rhs._v);
    return *this;
  }
  mint& operator/=(const mint& rhs){return *this = *this * rhs.inv();}
  mint operator+()const{return *this;}
  mint operator-()const{return mint() - *this;}
  mint pow(long long n)const{
    assert(0 <= n);
    mint x = *this, r = 1;
    while(n){
      if (n & 1) r *= x;
      x *= x;
      n >>= 1;
    }
    return r;
  }
  mint inv()const{
    auto eg = inv_gcd(_v, mod());
    assert(eg.first == 1);
    return eg.second;
  }
  friend mint operator+(const mint& lhs, const mint& rhs){return mint(lhs) += rhs;}
  friend mint operator-(const mint& lhs, const mint& rhs){return mint(lhs) -= rhs;}
  friend mint operator*(const mint& lhs, const mint& rhs){return mint(lhs) *= rhs;}
  friend mint operator/(const mint& lhs, const mint& rhs){return mint(lhs) /= rhs;}
  friend bool operator==(const mint& lhs, const mint& rhs){return lhs._v == rhs._v;}
  friend bool operator!=(const mint& lhs, const mint& rhs){return lhs._v != rhs._v;}
private:
  unsigned int _v;
  static barrett bt;
  static unsigned int umod(){return bt.umod();}
};
template <int id>
barrett dynamic_modint<id>::bt(998244353);
using modint = dynamic_modint<-1>;

using modint998244353 = static_modint<998244353>;
using modint1000000007 = static_modint<1000000007>;
template <class T>
using is_modint = std::is_base_of<modint_base, T>;
template <class T>
using is_modint_t = std::enable_if_t<is_modint<T>::value>;
template <class T>
using is_static_modint = std::is_base_of<static_modint_base, T>;
template <class T>
using is_static_modint_t = std::enable_if_t<is_static_modint<T>::value>;
template <class> struct is_dynamic_modint : public std::false_type {};
template <int id>
struct is_dynamic_modint<dynamic_modint<id>> : public std::true_type {};
template <class T>
using is_dynamic_modint_t = std::enable_if_t<is_dynamic_modint<T>::value>;
template<int m>
std::ostream &operator<<(std::ostream &dest, const static_modint<m> &a){
  dest << a.val();
  return dest;
}
template<int id>
std::ostream &operator<<(std::ostream &dest, const dynamic_modint<id> &a){
  dest << a.val();
  return dest;
}

// 0 <= n < m <= int_max
// 前処理 O(n + log(m))
// 各種計算 O(1)
// 変数 <= n
template<typename mint, is_modint<mint>* = nullptr>
struct modcomb{
private:
  int n;
  std::vector<mint> f, i, fi;
  void init(int _n){
    assert(0 <= _n && _n < mint::mod());
    if(_n < f.size()) return;
    n = _n;
    f.resize(n + 1), i.resize(n + 1), fi.resize(n + 1);
    f[0] = fi[0] = mint(1);
    if(n) f[1] = fi[1] = i[1] = mint(1);
    for(int j = 2; j <= n; j++) f[j] = f[j - 1] * j;
    fi[n] = f[n].inv();
    for(int j = n; j >= 2; j--){
      fi[j - 1] = fi[j] * j;
      i[j] = f[j - 1] * fi[j];
    }
  }
public:
  modcomb(): n(-1){}
  modcomb(int _n){
    init(_n);
  }
  void recalc(int _n){
    init(std::min(mint::mod() - 1, 1 << ceil_pow2(_n)));
  }
  mint comb(int a, int b){
    if((a < 0) || (b < 0) || (a < b)) return 0;
    return f[a] * fi[a - b] * fi[b];
  }
  mint combinv(int a, int b){
    assert(0 <= b && b <= a);
    return fi[a] * f[a - b] * f[b];
  }
  mint perm(int a, int b){
    if((a < 0) || (b < 0) || (a < b)) return 0;
    return f[a] * fi[a - b];
  }
  mint perminv(int a, int b){
    assert(0 <= b && b <= a);
    return fi[a] * f[a - b];
  }
  mint fac(int x){
    assert(0 <= x && x <= n);
    return f[x];
  }
  mint inv(int x){
    assert(0 < x && x <= n);
    return i[x];
  }
  mint finv(int x){
    assert(0 <= x && x <= n);
    return fi[x];
  }
};
template<typename mint>
mint combination_small_r(int n, int r){
  assert(r < mint::mod());
  if(n < r) return 0;
  mint res = 1, d = 1;
  for(int i = 0; i < r; i++) res *= n - i, d *= i + 1;
  return res / d;
}
template<typename mint, is_modint<mint>* = nullptr>
struct modpow_table{
  std::vector<mint> v;
  // x^maxkまで計算できる
  modpow_table(){}
  void init(int x, int maxk){
    v.resize(maxk + 1);
    v[0] = 1;
    for(int i = 1; i <= maxk; i++) v[i] = v[i - 1] * x;
  }
  mint pow(int k){
    assert(0 <= k && k < v.size());
    return v[k];
  }
};
// p/q
// mod : 素数
// 0 <= a < mod
// p, q < √modであるような既約分数p, qが存在する場合一意に定まる
std::pair<int, int> mod_frac_restore(int mod, int a){
  assert(0 <= a && a < mod);
  using mint = dynamic_modint<0>;
  mint::set_mod(mod);
  int sq = sqrt(mod);
  for(long long q = 1; q <= sq; q++){
    mint p = a * mint(q);
    if(p.val() <= sq) return {p.val(), q};
  }
  // 存在しない場合
  return {-1, -1};
}
#endif
