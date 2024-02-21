#ifndef _FPS_H_
#define _FPS_H_
#include "mod.hpp"
#include "convolution.hpp"
#include <queue>
#include <algorithm>

template<typename mint>
struct formal_power_series;

// 全ての要素の次数が同じだと仮定して総積を求める
// 全ての次数が1の場合 O(Nlog^2N)
template<typename mint>
formal_power_series<mint> multiply_constant_degree(std::vector<formal_power_series<mint>> v, int max_size = -1){
  if(max_size == -1) max_size = std::numeric_limits<int>::max();
  int n = v.size();
  for(int d = 1; d < n; d <<= 1){
    for(int j = 0; j < n; j += 2 * d){
      if(j + d < n){
        v[j] *= v[j + d];
        if(v[j].size() > max_size) v[j].resize(max_size);
      }
    }
  }
  return v[0];
}
// 次数が低い順に計算
// 次数の総和をNとしてO(Nlog^2N)
template<typename mint>
formal_power_series<mint> multiply_lower_degree(std::vector<formal_power_series<mint>> v, int max_size = std::numeric_limits<int>::max()){
  using p = std::pair<int, int>;
  std::priority_queue<p, std::vector<p>, std::greater<p>> que;
  int n = v.size();
  if(n == 0) return {mint(1)};
  for(int i = 0; i < n; i++){
    if(v[i].size() > max_size) v[i].resize(max_size);
    que.push({v[i].size(), i});
  }
  while(que.size() > 1){
    int a = que.top().second;
    que.pop();
    int b = que.top().second;
    que.pop();
    v[a] *= v[b];
    if(v[a].size() > max_size) v[a].resize(max_size);
    que.push({v[a].size(), a});
  }
  return v[que.top().second];
}

template<typename mint>
struct formal_power_series : std::vector<mint>{
  using std::vector<mint>::vector;
  using fps = formal_power_series<mint>;

  formal_power_series(const std::vector<mint> &v){
    this->resize(v.size());
    std::copy(v.begin(), v.end(), this->begin());
  }
  fps operator *= (const fps &vr){
    return *this = convolution_int_mod<mint>(*this, vr);
  }
  fps operator /= (fps &vr){
    return (*this) *= vr.inv();
  }
  fps operator += (const fps &vr){
    int n = this->size(), m = vr.size();
    if(n < m) this->resize(m);
    for(int i = 0; i < m; i++) (*this)[i] += vr[i];
    return *this;
  }
  fps operator -= (const fps &vr){
    int n = this->size(), m = vr.size();
    if(n < m) this->resize(m);
    for(int i = 0; i < m; i++) (*this)[i] -= vr[i];
    return *this;
  }
  fps operator %= (const fps &vr){
    int n = this->size(), m = vr.size();
    if(n < m) return *this;
    n = n - m + 1;
    fps r = ((rev().prefix(n) * vr.rev().inv(n)).prefix(n).rev(n)) * vr;
    return (*this - r).prefix(m - 1);
  }
  // 係数が0でない最大の次数, deg(0) = -1とする
  int deg()const{
    int n = this->size();
    for(int i = n - 1; i >= 0; i--) if((*this)[i].val()) return i;
    return -1;
  }
  // 係数が0でない最大の次数, deg(0) = -1とする, 余分な0を消す
  int deg_fix(){
    int n = this->size();
    for(int i = n - 1; i >= 0; i--){
      if(this->back().val() != 0) return i;
      else this->pop_back();
    }
    return -1;
  }
  // f(x) = q(x) * g(x) + r(x)となるような商q(x)と余りr(x)を求める
  // f(x)がn - 1次, g(x)がm - 1次のとき, q(x)はn - m次, r(x)はm - 2次以下になる
  std::pair<fps, fps> polynomial_div(fps g){
    int n = this->size(), m = g.size();
    if(deg_fix() < g.deg_fix()) return {fps{0}, *this};
    n = n - m + 1;
    fps q = (rev().prefix(n) * g.rev().inv(n)).prefix(n).rev(n);
    return {q, ((*this) - q * g).prefix(m - 1)};
  }
  // どちらかの次数が-1の場合空のfpsを返す, O(nm)
  static fps gcd_naive(const fps &f, const fps &g){
    if(g.deg() == -1) return g;
    if(f.deg() == -1) return f;
    fps st[2] = {g, f};
    bool flip = false;
    while(st[!flip].deg_fix() != -1){
      int n = st[flip].deg_fix(), m = st[!flip].deg_fix(), deg_diff = n - m;
      if(n >= m){
        mint c = st[flip][n] / st[!flip][m];
        for(int i = deg_diff; i <= n; i++) st[flip][i] -= c * st[!flip][i - deg_diff];
        n = st[flip].deg_fix(), m = st[!flip].deg_fix();
      }
      if(n < m) flip ^= 1;
    }
    return st[flip];
  }
  // f(x) * h(x)がmod g(x)でgcd(f(x), g(x))になるような{h(x), gcd(f(x), g(x))}を求める
  // O(nm)
  static std::pair<fps, fps> polynomial_inv_naive(const fps &f, const fps &g){
    if(f.deg() == -1) return std::make_pair(g, fps{0}); // fが0
    if(g.deg() == -1) return std::make_pair(fps{}, fps{}); // gで割れない
    fps st[2] = {g, f};
    fps ms[2] = {fps{0}, fps{1}};
    bool flip = false;
    while(st[!flip].deg_fix() != -1){
      int n = st[flip].deg_fix(), m = st[!flip].deg_fix(), deg_diff = n - m;
      if(n >= m){
        mint c = st[flip][n] / st[!flip][m];
        for(int i = deg_diff; i <= n; i++) st[flip][i] -= c * st[!flip][i - deg_diff];
        int o = ms[!flip].deg_fix() + 1;
        if(o + deg_diff > ms[flip].size()) ms[flip].resize(o + deg_diff, 0);
        for(int i = 0; i < o; i++) ms[flip][i + deg_diff] -= ms[!flip][i] * c;
        n = st[flip].deg_fix(), m = st[!flip].deg_fix();
      }
      if(n < m) flip ^= 1;
    }
    return {st[flip], ms[flip]};
  }
  fps operator += (const mint &vr){
    int n = this->size();
    if(n == 0) this->resize(1, 0);
    (*this)[0] += vr;
    return *this;
  }
  fps operator -= (const mint &vr){
    int n = this->size();
    if(n == 0) this->resize(1, 0);
    (*this)[0] -= vr;
    return *this;
  }
  fps operator *= (const mint &vr){
    int n = this->size();
    for(int i = 0; i < n; i++) (*this)[i] *= vr;
    return *this;
  }
  fps operator /= (const mint &vr){
    mint ir = vr.inv();
    int n = this->size();
    for(int i = 0; i < n; i++) (*this)[i] *= ir;
    return *this;
  }
  fps operator + (const fps& vr)const{return fps(*this) += vr;}
  fps operator - (const fps& vr)const{return fps(*this) -= vr;}
  fps operator * (const fps& vr)const{return fps(*this) *= vr;}
  fps operator / (const fps& vr)const{return fps(*this) /= vr;}
  fps operator % (const fps& vr)const{return fps(*this) %= vr;}
  fps operator + (const mint& vr)const{return fps(*this) += vr;}
  fps operator - (const mint& vr)const{return fps(*this) -= vr;}
  fps operator * (const mint& vr)const{return fps(*this) *= vr;}
  fps operator / (const mint& vr)const{return fps(*this) /= vr;}

  void print(int printsize = -1)const{
    if(printsize == -1) printsize = (int)this->size();
    printsize = std::min(printsize, (int)this->size());
    if(printsize == 0){
      std::cout << '\n';
      return;
    }
    for(int i = 0; i < printsize; i++){
      if(i == printsize - 1) std::cout << (*this)[i].val() << '\n';
      else std::cout << (*this)[i].val() << ' ';
    }
  }
  fps rev(int deg = -1)const{
    fps res(*this);
    if(deg != -1) res.resize(deg, 0);
    std::reverse(res.begin(), res.end());
    return res;
  }
  fps prefix(int deg){
    int n = std::min((int)this->size(), deg);
    return fps(this->begin(), this->begin() + n);
  }
  fps rshift(int s){
    fps res(*this);
    if(res.size() <= s) return {};
    res.erase(res.begin(), res.begin() + s);
    return res;
  }
  fps lshift(int s){
    fps res(*this);
    res.insert(res.begin(), s, mint(0));
    return res;
  }
  fps inv_any_mod(int deg = -1){
    assert((*this)[0].val());
    int n = this->size();
    if(deg == -1) deg = n;
    fps res{(*this)[0].inv()};
    for(int i = 1; i < deg; i <<= 1){
      res = (res + res - (res * res * prefix(i << 1))).prefix(i << 1);
    }
    return res.prefix(deg);
  }
  fps inv(int deg = -1){
    assert((*this)[0].val());
    if(mint::mod() != 998244353) return inv_any_mod(deg);
    int n = this->size();
    if(deg == -1) deg = n;
    fps res{(*this)[0].inv()};
    for(int i = 1, z = i << 2; i < deg; i <<= 1, z <<= 1){
      std::vector<mint> res_old = res, f_prefix = prefix(i << 1);
      res_old.resize(z);
      butterfly(res_old);
      f_prefix.resize(z);
      butterfly(f_prefix);
      for(int j = 0; j < z; j++) res_old[j] *= res_old[j] * f_prefix[j];
      butterfly_inv(res_old);
      mint iz = mint(z).inv();
      res_old.resize(z >> 1);
      for(int j = 0; j < (z >> 1); j++) res_old[j] *= iz;
      res = res + res - (fps)res_old;
    }
    return res.prefix(deg);
  }
  fps diff(){
    int n = (int) this->size();
    fps res(std::max(0, n - 1));
    for(int i = 1; i < n; i++) res[i - 1] = (*this)[i] * i;
    return res;
  }
  fps integral(){
    static modcomb<mint> mcb;
    int n = (int) this->size();
    fps res(n + 1);
    mcb.recalc(n);
    res[0] = 0;
    for(int i = 0; i < n; i++) res[i + 1] = (*this)[i] * mcb.inv(i + 1);
    return res;
  }
  fps log(int deg = -1){
    assert((*this)[0].val() == 1);
    int n = (int)this->size();
    if(deg == -1) deg = n;
    return (this->diff() * this->inv(deg)).prefix(deg - 1).integral();
  }
  fps exp_any_mod(int deg = -1){
    assert((*this)[0].val() == 0);
    int n = this->size();
    if(deg == -1) deg = n;
    fps res{1};
    for(int i = 1; i < deg; i <<= 1){
      res = (res * (prefix(i << 1) + mint(1) - res.log(i << 1))).prefix(i << 1);
    }
    return res.prefix(deg);
  }
  fps exp(int deg = -1){
    if(mint::mod() != 998244353) return exp_any_mod(deg);
    assert((*this)[0].val() == 0);
    int n = this->size();
    if(deg == -1) deg = n;
    fps res{1}, diff_res, diff_res_prev, inv_res, inv_res_old, inv_res_prev, log_res, log_res_prev;
    for(int i = 1; i < deg; i <<= 1){
      diff_res_prev = diff_res;
      diff_res.resize(i - 1);
      for(int j = std::max(1, i >> 1); j < i; j++) diff_res[j - 1] = res[j] * j;
      if(i == 1){
        inv_res = inv_res_old = {1};
      }else{
        fps i1 = inv_res_old;
        i1.resize(i << 1);
        butterfly(i1);
        fps i2 = res.prefix(i);
        i2.resize(i << 1);
        butterfly(i2);
        for(int j = 0; j < (i << 1); j++) i1[j] *= i1[j] * i2[j];
        butterfly_inv(i1);
        i1.resize(i);
        mint iz = mint(i << 1).inv();
        for(int j = 0; j < i; j++) i1[j] *= iz;
        inv_res_old = inv_res_old + inv_res_old - i1;

        i1 = inv_res_old;
        i1.resize(i << 2);
        butterfly(i1);
        i2 = res.prefix(i << 1);
        i2.resize(i << 2);
        butterfly(i2);
        for(int j = 0; j < (i << 2); j++) i1[j] *= i1[j] * i2[j];
        butterfly_inv(i1);
        i1.resize(i << 1);
        iz = mint(i << 2).inv();
        for(int j = 0; j < (i << 1); j++) i1[j] *= iz;
        inv_res = inv_res_old + inv_res_old - i1;
      }
      log_res = (diff_res * inv_res).prefix((i << 1) - 1).integral();
      res = (res * (prefix(i << 1) + mint(1) - log_res)).prefix(i << 1);
    }
    return res.prefix(deg);
  }
  fps pow(long long k, int deg = -1){
    int n = (int) this->size();
    if(deg == -1) deg = n;
    if(!k){
      fps res(deg, 0);
      res[0] = 1;
      return res;
    }
    for(long long i = 0; i < n; i++){
      if((*this)[i] != 0) {
        fps C = (*this) * (*this)[i].inv();
        fps D(deg);
        for(int j = i; j < n; j++) D[j - i] = C[j];
        D = (D.log(deg) * mint(k)).exp(deg) * (*this)[i].pow(k);
        fps E(deg);
        if(i * k > deg) return E;
        for(long long j = 0; j + __int128_t(i) * k < deg && j < D.size(); j++) E[j + i * k] = D[j];
        return E;
      }
    }
    return *this;
  }
  //初めて0以外が現れるのが奇数次数 or 初めに現れる数が平方根でない -> 解無し(空のfpsを返す)
  fps sqrt(int deg = -1){
    int n = (int)this->size();
    if(deg == -1) deg = n;
    if((*this)[0].val()==0){
      for(int i = 1; i < n;i++){
        if((*this)[i].val() != 0){
          if(i & 1) return {};
          if(deg - i / 2 <= 0) return fps(deg, 0);
          fps res = rshift(i).sqrt(deg - i / 2);
          if(res.empty()) return {};
          res = res.lshift(i / 2);
          if(res.size() < deg) res.resize(deg, 0);
          return res;
        }
      }
      return fps(deg, 0);
    }
    int x = modsqrt((*this)[0].val(), mint::mod());
    if(x == -1) return {};
    fps res({mint(x)});
    mint i2 = mint(2).inv();
    for(int i = 1; i < deg; i <<= 1) res = (res + prefix(i << 1) * res.inv(i << 1)) * i2;
    return res;
  }
  // reference: http://www.eecs.harvard.edu/~htk/publication/1978-jacm-brent-kung.pdf
  // N次多項式 f(x)にM次多項式 g(x)を合成した結果を求める
  // deg(f) = deg(g) = deg(ans) = Nとして、
  // f(x)をk := √N ブロックに平方分割すると　f(x) = f_0(x) + f_1(x) x^k ... と約k項になる
  // f((g(x))) = f_0(g(x)) + f_1(g(x))g(x)^k + ...
  // 1. f_i(g(x))はk項とM項の合成なので、g(x)^i (0<=i<=k)を前処理しておくとO(Nk)
  // 2. g(x)^kiを求める、かけるのは共に一回あたりO(NlogN)
  // (1, 2)をkブロック分行うのでO(k × (Nk + NlogN)) = O(N^2 + N^1.5 * logN)
  static fps composition(const fps &f, const fps &g, int deg = -1){
    int n = f.size();
    //if(deg == -1) deg = (n - 1) * (m - 1) + 1;
    if(deg == -1) deg = n;
    int k = (int)std::sqrt(n) + 1; // ブロック数
    int d = (n + k - 1) / k; // 1ブロックあたりの次数
    std::vector<fps> gpow(k + 1);
    gpow[0] = {mint(1)};
    for(int i = 1; i <= k; i++){
      gpow[i] = gpow[i - 1] * g;
      if(gpow[i].size() > deg) gpow[i].resize(deg);
    }
    std::vector<fps> fi(k, fps(deg, 0));
    for(int i = 0; i < k; i++){
      for(int j = 0; j < d; j++){
        int idx = i * d + j;
        if(idx >= n) break;
        for(int t = 0; t < gpow[j].size(); t++) fi[i][t] += gpow[j][t] * f[idx];
      }
    }
    fps res(deg, 0), gd = {1};
    for(int i = 0; i < k; i++){
      fi[i] *= gd;
      int sz = std::min(deg, (int)fi[i].size());
      for(int j = 0; j < sz; j++) res[j] += fi[i][j];
      gd *= gpow[d];
      if(gd.size() > deg) gd.resize(deg);
    }
    return res;
  }
};
#endif
