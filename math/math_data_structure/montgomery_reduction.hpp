#ifndef _MONTGOMERY_REDUCTION_H_
#define _MONTGOMERY_REDUCTION_H_
#include <vector>
#include <algorithm>
#include <cstdint>
#include <algorithm>

struct barrett_reduction{
  unsigned int mod;
  unsigned long long m;
  barrett_reduction(unsigned int _mod) : mod(_mod){
    m = ((__uint128_t)1 << 64) / mod;
  }
  unsigned int reduce(unsigned int a){
    unsigned long long q = ((__uint128_t)a * m) >> 64;
    a -= q * mod; // 0 <= a < 2 * mod
    // return a;
    return a >= mod ? a - mod : a;
  }
  unsigned int mul(unsigned int a, unsigned int b){
    return reduce((unsigned long long)a * b);
  }
  // {gcd(a, mod), x}, such that a * x ≡ gcd(a, mod)
  std::pair<unsigned int, unsigned int> inv(unsigned int a){
    if(a >= mod) a = reduce(a);
    if(a == 0) return {mod, 0};
    unsigned int s = mod, t = a;
    long long m0 = 0, m1 = 1;
    while(t){
      int u = s / t;
      s -= t * u;
      m0 -= m1 * u;
      std::swap(m0, m1);
      std::swap(s, t);
    }
    if(m0 < 0) m0 += mod / s;
    return {s, m0};
  }
};
// 64bit mod対応
// R = 2^64
// 偶数modだと壊れる
struct montgomery_reduction_64bit{
private:
  // [0, 2 * MOD)
  inline uint64_t reduce_unsafe(__uint128_t x) const{
    x = (x + ((uint64_t)x * (uint64_t)NEG_INV) * MOD) >> 64;
    return x;
  }
  void _set_mod(uint64_t mod){
    assert((mod < (1ULL << 63)) && (mod & 1));
    MOD = mod;
    NEG_INV = 0;
    __uint128_t s = 1, t = 0;
    for(int i = 0; i < 64; i++){
      if (~t & 1) {
        t += MOD;
        NEG_INV += s;
      }
      t >>= 1;
      s <<= 1;
    }
    R2 = ((__uint128_t)1 << 64) % MOD;
    R2 = R2 * R2 % MOD;
    ONE = generate(1);
  }
  __uint128_t MOD, NEG_INV, R2;
  uint64_t ONE;
public:
  montgomery_reduction_64bit(){}
  montgomery_reduction_64bit(uint64_t mod){_set_mod(mod);}
  void set_mod(uint64_t mod){
    _set_mod(mod);
  }
  uint64_t mod()const{
    return MOD;
  }
  inline uint64_t one()const{
    return ONE;
  }
  inline uint64_t generate(uint64_t x)const{
    return reduce((__uint128_t)x * R2);
  }
  inline uint64_t reduce(__uint128_t x)const{
    x = (x + ((uint64_t)x * (uint64_t)NEG_INV) * MOD) >> 64;
    return x < MOD ? x : x - MOD;
  }
  inline uint64_t fix(uint64_t x)const{
    return x < MOD ? x : x - MOD;
  }
  // [0, 2MOD)
  inline uint64_t mul(uint64_t mx, uint64_t my)const{
    return reduce_unsafe((__uint128_t)mx * my);
  }
  inline uint64_t mul_safe(uint64_t mx, uint64_t my)const{
    return reduce((__uint128_t)mx * my);
  }
  // [0, 2MOD)
  inline uint64_t add(uint64_t mx, uint64_t my)const{
    return (mx >= MOD ? mx - MOD : mx) + (my >= MOD ? my - MOD : my);
  }
  inline uint64_t add_safe(uint64_t mx, uint64_t my)const{
    uint64_t res = (mx >= MOD ? mx - MOD : mx) + (my >= MOD ? my - MOD : my);
    return res >= MOD ? res - MOD : res;
  }
  // [0, 2MOD)
  inline uint64_t sub(uint64_t mx, uint64_t my)const{
    if(my >= MOD) my -= MOD;
    return mx >= my ? mx - my : mx + MOD - my;
  }
  inline uint64_t sub_safe(uint64_t mx, uint64_t my)const{
    if(my >= MOD) my -= MOD;
    uint64_t res = mx >= my ? mx - my : mx + MOD - my;
    return res >= MOD ? res - MOD : res;
  }
  inline uint64_t pow(uint64_t ma, uint64_t b)const{
    uint64_t ret = one();
    while(b){
      if(b & 1) ret = mul(ret, ma);
      ma = mul(ma, ma);
      b >>= 1;
    }
    return ret;
  }
  inline uint64_t pow_safe(uint64_t ma, uint64_t b)const{
    return fix(pow(ma, b));
  }
};
unsigned long long mod_pow_mr(unsigned long long a, unsigned long long b, unsigned long long m){
  montgomery_reduction_64bit mr(m);
  return mr.reduce(mr.pow(mr.generate(a), b));
}
#endif