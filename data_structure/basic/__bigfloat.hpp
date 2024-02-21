#ifndef _BIGFLOAT_H_
#define _BIGFLOAT_H_
#include "../../math/mod.hpp"
#include "../../math/convolution.hpp"
#include <string>
#include <cassert>
#include <iostream>
#include <algorithm>

struct bigfloat : std::vector<unsigned long long>{
  static int prec;
  // 使う前に初期化する　
  static void set_prec(int p){
    prec = p;
  }
  bool s;
  int exp_10 = 0;
private:
  using ull = unsigned long long;
  using ull128 = __uint128_t;
  using std::vector<ull>::vector;
  static constexpr ull base = 1000000000000000000ULL; // 10^18
  static constexpr int base_log_10 = 18;
  void set_zero(){
    this->clear();
    this->push_back(0);
    exp_10 = 0;
    s = false;
  }
  static std::vector<unsigned long long> convolution_bigfloat(const std::vector<unsigned long long> &a, const std::vector<unsigned long long> &b){
    static constexpr ull128 base_sq = (ull128)base * base;
    int n = int(a.size()), m = int(b.size());
    if(!n || !m) return {};
    if(std::min(n, m) <= 60){
      std::vector<ull128> c2(n + m, 0);
      for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
          c2[i + j] += (ull128)a[i] * b[j];
          if(c2[i + j] >= base_sq) c2[i + j] -= base_sq, c2[i + j + 1] += base;
        }
      }
      ull up = 0;
      std::vector<ull> c(n + m, 0);
      for(int i = 0; i < n + m; i++){
        ull128 tmp = c2[i] + up;
        up = tmp / base;
        c[i] = tmp - (ull128)up * base;
      }
      assert(up == 0);
      if(c.back() == 0) c.pop_back();
      return c;
    }
    static constexpr unsigned long long MOD1 = 754974721;  // 2^24
    static constexpr unsigned long long MOD2 = 998244353;  // 2^23
    static constexpr unsigned long long MOD3 = 167772161;  // 2^25
    static constexpr unsigned long long MOD4 = 469762049;  // 2^26
    static constexpr unsigned long long MOD5 = 1224736769; // 2^24
    static constexpr unsigned long long M1M2 = MOD1 * MOD2;
    static constexpr __uint128_t M1M2M3 = (__uint128_t)M1M2 * MOD3;
    static constexpr __uint128_t M1M2M3M4 = M1M2M3 * MOD4;
    static constexpr __uint128_t Q = M1M2M3M4 / base;
    static constexpr __uint128_t R = M1M2M3M4 % base;
    static constexpr unsigned long long IM1_M2 = inv_gcd(MOD1, MOD2).second;
    static constexpr unsigned long long IM1M2_M3 = inv_gcd(MOD1 * MOD2, MOD3).second;
    static constexpr unsigned long long IM1M2M3_M4 = inv_gcd((long long)(M1M2M3 % MOD4), MOD4).second;
    static constexpr unsigned long long IM1M2M3M4_M5 = inv_gcd((long long)(M1M2M3M4 % MOD5), MOD5).second;
    auto c1 = convolution_mod<MOD1>(a, b);
    auto c2 = convolution_mod<MOD2>(a, b);
    auto c3 = convolution_mod<MOD3>(a, b);
    auto c4 = convolution_mod<MOD4>(a, b);
    auto c5 = convolution_mod<MOD5>(a, b);
    std::vector<unsigned long long> c(n + m - 1);
    __uint128_t up = 0;
    for(int i = 0; i < n + m - 1; i++){
      __uint128_t x = c1[i], y, t;
      t = ((c2[i] < x ? MOD2 + c2[i] - x : c2[i] - x) * IM1_M2) % MOD2;
      x += t * MOD1; // x < 2^60
      y = x % MOD3;
      t = ((c3[i] < y ? MOD3 + c3[i] - y : c3[i] - y) * IM1M2_M3) % MOD3;
      x += t * M1M2; // x < 2^90
      y = x % MOD4;
      t = ((c4[i] < y ? MOD4 + c4[i] - y : c4[i] - y) * IM1M2M3_M4) % MOD4;
      x += t * M1M2M3; // x < 2^120
      y = x % MOD5;
      t = ((c5[i] < y ? MOD5 + c5[i] - y : c5[i] - y) * IM1M2M3M4_M5) % MOD5;
      x += up + t * R;
      up = Q * t;
      __uint128_t z = x / base;
      x -= base * z;
      up += z;
      c[i] = x;
    }
    assert(up < base);
    if(up) c.push_back(up);
    return c;
  }
  static std::vector<unsigned long long> square_bigfloat(const std::vector<unsigned long long> &a){
    static constexpr ull128 base_sq2 = (ull128)base * base * 2;
    int n = int(a.size());
    if(!n) return {};
    if(n <= 60){
      std::vector<ull128> c2(2 * n, 0);
      for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
          c2[i + j] += 2 * ull128(a[i]) * a[j];
          if(c2[i + j] >= base_sq2) c2[i + j] -= base_sq2, c2[i + j + 1] += (base << 1);
        }
        c2[2 * i] += ull128(a[i]) * a[i];
        if(c2[2 * i] >= base_sq2) c2[2 * i] -= base_sq2, c2[2 * i + 1] += (base << 1);
      }
      ull up = 0;
      std::vector<ull> c(2 * n, 0);
      for(int i = 0; i < 2 * n; i++){
        ull128 tmp = c2[i] + up;
        up = tmp / base;
        c[i] = tmp - (ull128)up * base;
      }
      assert(up == 0);
      if(c.back() == 0) c.pop_back();
      return c;
    }
    static constexpr unsigned long long MOD1 = 754974721;  // 2^24
    static constexpr unsigned long long MOD2 = 998244353;  // 2^23
    static constexpr unsigned long long MOD3 = 167772161;  // 2^25
    static constexpr unsigned long long MOD4 = 469762049;  // 2^26
    static constexpr unsigned long long MOD5 = 1224736769; // 2^24
    static constexpr unsigned long long M1M2 = MOD1 * MOD2;
    static constexpr __uint128_t M1M2M3 = (__uint128_t)M1M2 * MOD3;
    static constexpr __uint128_t M1M2M3M4 = M1M2M3 * MOD4;
    static constexpr __uint128_t Q = M1M2M3M4 / base;
    static constexpr __uint128_t R = M1M2M3M4 % base;
    static constexpr unsigned long long IM1_M2 = inv_gcd(MOD1, MOD2).second;
    static constexpr unsigned long long IM1M2_M3 = inv_gcd(MOD1 * MOD2, MOD3).second;
    static constexpr unsigned long long IM1M2M3_M4 = inv_gcd((long long)(M1M2M3 % MOD4), MOD4).second;
    static constexpr unsigned long long IM1M2M3M4_M5 = inv_gcd((long long)(M1M2M3M4 % MOD5), MOD5).second;
    auto c1 = square_mod<MOD1>(a);
    auto c2 = square_mod<MOD2>(a);
    auto c3 = square_mod<MOD3>(a);
    auto c4 = square_mod<MOD4>(a);
    auto c5 = square_mod<MOD5>(a);
    std::vector<unsigned long long> c(2 * n - 1, 0);
    __uint128_t up = 0;
    for(int i = 0; i < 2 * n - 1; i++){
      __uint128_t x = c1[i], y, t;
      t = ((c2[i] < x ? MOD2 + c2[i] - x : c2[i] - x) * IM1_M2) % MOD2;
      x += t * MOD1; // x < 2^60
      y = x % MOD3;
      t = ((c3[i] < y ? MOD3 + c3[i] - y : c3[i] - y) * IM1M2_M3) % MOD3;
      x += t * M1M2; // x < 2^90
      y = x % MOD4;
      t = ((c4[i] < y ? MOD4 + c4[i] - y : c4[i] - y) * IM1M2M3_M4) % MOD4;
      x += t * M1M2M3; // x < 2^120
      y = x % MOD5;
      t = ((c5[i] < y ? MOD5 + c5[i] - y : c5[i] - y) * IM1M2M3M4_M5) % MOD5;
      x += up + t * R;
      up = Q * t;
      __uint128_t z = x / base;
      x -= base * z;
      up += z;
      c[i] = x;
    }
    assert(up < base);
    if(up) c.push_back(up);
    return c;
  }
  // 

  bigfloat square()const{
    auto v = square_bigfloat(*this);
    bigfloat ret;
    ret.exp_10 = 2 * exp_10; // 符号は正, exp_10は2倍になる
    
    int m = v.size();
    ret.resize(m);
    for(int i = 0; i < m; i++) ret[i] = v[i];
    return ret;
  }
  static void add(bigfloat &a, bigfloat b){
    int e = std::min(a.exp_10, b.exp_10);
    if(e < a.exp_10){
      a.lshift_inplace(a.exp_10 - e);
      a.exp_10 = e;
    }else if(e < b.exp_10){
      b.lshift_inplace(b.exp_10 - e);
    }
    int n = a.size(), m = b.size();
    if(n < m) a.resize(m, 0);
    bool up = 0;
    for(int i = 0; i < m; i++){
      if(a[i] + up < base - b[i]){
        a[i] += b[i] + up;
        up = 0;
      }else{
        a[i] = (a[i] + up) - (base - b[i]);
        up = 1;
      }
    }
    for(int i = m; i < n; i++){
      if(!up) break;
      if(a[i] == base - 1) a[i] = 0;
      else{
        a[i]++;
        up = 0;
        break;
      }
    }
    // 桁が増える
    if(up) a.push_back(1);
  }
  static void sub(bigfloat &a, bigfloat b){
    int e = std::min(a.exp_10, b.exp_10);
    if(e < a.exp_10){
      a.lshift_inplace(a.exp_10 - e);
      a.exp_10 = e;
    }else if(e < b.exp_10){
      b.lshift_inplace(b.exp_10 - e);
    }
    int n = a.size(), m = b.size(), k = -1;
    bool is_this_bigger = false;
    for(int i = std::max(n, m) - 1; i >= 0; i--){
      ull _l = (i >= n ? 0 : a[i]), _r = (i >= m ? 0 : b[i]);
      if(_l != _r){
        k = i;
        is_this_bigger = _l > _r;
        break;
      }
    }
    // this == r
    if(k == -1){
      a.set_zero();
      return;
    }
    a.resize(k + 1, 0);
    if(is_this_bigger){ // this > r, n >= k
      bool up = false;
      for(int i = 0; i <= k; i++){
        ull sub = i >= m ? 0 : b[i];
        if(a[i] < sub + up){
          a[i] += base - sub - up;
          up = 1;
        }else{
          a[i] -= sub + up;
          up = 0;
        }
      }
    }else{// this < r, m >= k
      a.s ^= 1;
      bool up = false;
      for(int i = 0; i <= k; i++){
        ull sub = i >= n ? 0 : a[i];
        if(b[i] < sub + up){
          a[i] = base - sub - up + b[i];
          up = 1;
        }else{
          a[i] = b[i] - sub - up;
          up = 0;
        }
      }
    }
  }
  bigfloat inv_newton(int p)const{
    int d = __deg();
    if(d == -1){
      std::cerr << "Divide by zero" << '\n';
      exit(1);
    }
    bigfloat x(1ULL);
    x.exp_10 = -(d + 1);
    bigfloat pre = __prefix(base_log_10);
    for(int i = 0; i < 5; i++){
      bigfloat y = pre * x.square();
      x = (x + x - y);
      x = x.__prefix(base_log_10);
    }
    x = x.__prefix(base_log_10);
    for(int i = 1; i < p; i <<= 1){
      bigfloat y = __prefix(std::max(base_log_10, i << 1)) * x.square();
      x = (x + x - y);
      x = x.__prefix(std::max(base_log_10, i << 1));
    }
    return x.__prefix(p);
  }
  int __deg()const{
    int n = this->size();
    for(int i = n - 1; i >= 0; i--){
      if((*this)[i]){
        ull x = base / 10;
        for(int j = base_log_10 - 1; j >= 0; j--){
          if((*this)[i] >= x) return (i * base_log_10) + j;
          x /= 10;
        }
      }
    }
    return -1; // 0
  }
  // 10進数でk桁ずらす
  void lshift_inplace(int k){
    static std::vector<ull> ten_pow;
    if(ten_pow.empty()){
      ten_pow.resize(base_log_10 + 1);
      ten_pow[0] = 1;
      for(int i = 1; i <= base_log_10; i++) ten_pow[i] = ten_pow[i - 1] * 10;
    }
    assert(k >= 0);
    int n = this->size();
    this->resize(n + (k + base_log_10 - 1) / base_log_10, 0);
    int b = k / base_log_10;
    k %= base_log_10;
    n = this->size();
    if(!k){
      for(int i = n - 1; i >= b; i--) (*this)[i] = (*this)[i - b];
      for(int i = b - 1; i >= 0; i--) (*this)[i] = 0;
    }else{
      int upper = base_log_10 - k;
      for(int i = n - 1; i >= 0; i--){
        if(i < b) (*this)[i] = 0;
        else if(i == b){
          (*this)[i] = ((__uint128_t)(*this)[i - b] % ten_pow[upper]) * ten_pow[k];
        }else{
          (*this)[i] = ((__uint128_t)(*this)[i - b] % ten_pow[upper]) * ten_pow[k] + ((*this)[i - b - 1] / ten_pow[upper]);
        }
      }
    }
  }
  // 10進数でk桁ずらす
  bigfloat rshift(int k)const{
    static std::vector<ull> ten_pow;
    if(ten_pow.empty()){
      ten_pow.resize(base_log_10 + 1);
      ten_pow[0] = 1;
      for(int i = 1; i <= base_log_10; i++) ten_pow[i] = ten_pow[i - 1] * 10;
    }
    assert(k >= 0);
    int b = k / base_log_10;
    k %= base_log_10;
    int n = this->size(), m = n - b;
    bigfloat res;
    if(n <= b){
      res.set_zero();
      return res;
    }
    res.s = s;
    res.exp_10 = exp_10;
    res.resize(m, 0);
    if(!k){
      for(int i = 0; i < m; i++) res[i] = (*this)[i + b];
    }else{
      int upper = base_log_10 - k;
      for(int i = 0; i < m; i++){
        if(i + b == n - 1){
          res[i] = (*this)[i + b] / ten_pow[k];
        }else{
          res[i] = (*this)[i + b] / ten_pow[k] + ((__uint128_t)(*this)[i + b  + 1] % ten_pow[k]) * ten_pow[upper];
        }
      }
    }
    return res;
  }
  // 先頭d桁
  bigfloat __prefix(int d)const{
    int deg = __deg();
    int n = std::min(deg, d);
    bigfloat res = rshift(deg - n);
    res.exp_10 += deg - n;
    return res;
  }
  std::string __str()const{
    std::string res = s ? "-" : "";
    int n = this->size();
    bool f = false;
    for(int i = n - 1; i >= 0; i--){
      if((*this)[i] && !f){ // omit leading zeros
        f = true;
        res += std::to_string((*this)[i]);
      }else if(f){ // use leading zeros
        std::string t = std::to_string((*this)[i]);
        res += std::string(base_log_10 - t.size(), '0') + t;
      }
    }
    if(!f) res = '0';
    return res;
  }
  bigfloat div1(const bigfloat &vr){
    if(vr > (*this)){
      set_zero();
      return *this;
    }
    int d = __deg();
    if(d < base_log_10) return (*this)[0] / vr[0];
    return (*this) *= vr.inv_newton(d + 1);
  }
public:
  // return {this / vr, this % vr}
  std::pair<bigfloat, bigfloat> div2(const bigfloat &vr)const{
    bigfloat q = bigfloat(*this).div1(vr);
    q = q.rshift(-q.exp_10);
    q.exp_10 = 0;
    bigfloat r = *this - (vr * q);
    if(vr <= r){
      r -= vr;
      q = q + 1LL;
    }
    return {q, r};
  }
  bigfloat operator += (const bigfloat &r){
    if(s != r.s){
      s ^= 1;
      (*this) -= r;
      s ^= 1;
      return *this;
    }
    add(*this, r);
    return *this;
  }
  bigfloat operator -= (const bigfloat &r){
    if(s != r.s){
      s ^= 1;
      (*this) += r;
      s ^= 1;
      return *this;
    }
    sub(*this, r);
    return *this;
  }
  bigfloat operator *= (const bigfloat &r){
    if(this->empty() || r.empty()){
      set_zero();
      return *this;
    }
    s ^= r.s;
    exp_10 += r.exp_10;
    auto v = convolution_bigfloat(*this, r);
    int n = v.size();
    this->resize(n);
    for(int i = 0; i < n; i++) (*this)[i] = v[i];
    return *this;
  }
  bigfloat operator /= (const bigfloat &vr){return *this = div2(vr).first;}
  bigfloat operator %= (const bigfloat &vr){return *this = div2(vr).second;}
  bigfloat operator + (const bigfloat& vr)const{return bigfloat(*this) += vr;}
  bigfloat operator - (const bigfloat& vr)const{return bigfloat(*this) -= vr;}
  bigfloat operator * (const bigfloat& vr)const{return bigfloat(*this) *= vr;}
  bigfloat operator / (const bigfloat &vr)const{return div2(vr).first;}
  bigfloat operator % (const bigfloat &vr)const{return div2(vr).second;}
  bigfloat operator += (long long r){return *this += bigfloat(r);}
  bigfloat operator -= (long long r){return *this -= bigfloat(r);}
  bigfloat operator *= (long long r){return *this *= bigfloat(r);}
  bigfloat operator /= (long long r){return *this /= bigfloat(r);}
  bigfloat operator %= (long long r){return *this %= bigfloat(r);}
  bigfloat operator + (long long r){return bigfloat(*this) += bigfloat(r);}
  bigfloat operator - (long long r){return bigfloat(*this) -= bigfloat(r);}
  bigfloat operator * (long long r){return bigfloat(*this) *= bigfloat(r);}
  bigfloat operator / (long long r){return bigfloat(*this) /= bigfloat(r);}
  bigfloat operator % (long long r){return bigfloat(*this) %= bigfloat(r);}
  bool operator > (const bigfloat &vr)const{
    int n = this->size(), m = vr.size();
    for(int i = std::max(n, m) - 1; i >= 0; i--){
      ull x = i >= n ? 0 : (*this)[i];
      ull y = i >= m ? 0 : vr[i];
      if(x > y) return true;
      if(x < y) return false;
    }
    return false;
  }
  bool operator <= (const bigfloat &vr)const{
    return !((*this) > vr);
  }
  bigfloat(): s(0){}
  bigfloat(ull r) : s(0){
    this->push_back(r);
  }
  bigfloat(long long r){
    if(r < 0){
      s = 1;
      r *= -1;
    }else s = 0;
    this->push_back(r);
  }
  bigfloat(std::string &t) : s(0){
    int n = t.size();
    int _st = 0;
    if(t[0] == '-'){
      s = 1;
      _st = 1;
    }
    std::string tmp = "";
    for(int i = n - 1; i >= _st; i--){
      tmp += t[i];
      if(tmp.size() == base_log_10){
        std::reverse(tmp.begin(), tmp.end());
        this->push_back(std::stoull(tmp));
        tmp = "";
      }
    }
    if(tmp != ""){
      std::reverse(tmp.begin(), tmp.end());
      this->push_back(std::stoull(tmp));
    }
  }
  std::string str()const{
    return __str();
  }
};
#endif