#ifndef _MODINT64BIT_H_
#define _MODINT64BIT_H_
#include "prime.hpp"
#include "math_data_structure/montgomery_reduction.hpp"

// @param mod < 2^63, modは素数
template<int id> 
struct dynamic_modint64bit{
  using mint = dynamic_modint64bit;
public:
  static long long mod(){return (long long)mr.mod();}
  static void set_mod(long long m){
    assert(1 <= m && m < (1ULL << 63));
    mr = montgomery_reduction_64bit(m);
    assert(_miller_rabin_mr(m, mr));
  }
  static mint raw(unsigned long long v){
    mint x;
    x._v = v;
    return x;
  }
  dynamic_modint64bit(): _v(0){}
  template <class T>
  dynamic_modint64bit(T v){
    long long x = v % (long long)(mod());
    if(x < 0) x += mod();
    _v = mr.generate(x);
  }
  unsigned long long val()const{return mr.reduce(_v);}
  mint& operator++(){
    _v++;
    if(_v == (umod() << 1)) _v = 0;
    return *this;
  }
  mint& operator--(){
    if(_v == 0) _v = umod();
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
    _v = mr.add(_v, rhs._v);
    return *this;
  }
  mint& operator-=(const mint& rhs){
    _v = mr.sub(_v, rhs._v);
    return *this;
  }
  mint& operator*=(const mint& rhs){
    _v = mr.mul(_v, rhs._v);
    return *this;
  }
  mint& operator/=(const mint& rhs){return *this = *this * rhs.inv();}
  mint operator+()const{return *this;}
  mint operator-()const{return mint() - *this;}
  mint pow(long long n)const{
    assert(0 <= n);
    return raw(mr.pow(_v, n));
  }
  mint inv()const{
    return raw(mr.pow(_v, mod() - 2));
  }
  friend mint operator+(const mint& lhs, const mint& rhs){return mint(lhs) += rhs;}
  friend mint operator-(const mint& lhs, const mint& rhs){return mint(lhs) -= rhs;}
  friend mint operator*(const mint& lhs, const mint& rhs){return mint(lhs) *= rhs;}
  friend mint operator/(const mint& lhs, const mint& rhs){return mint(lhs) /= rhs;}
  friend bool operator==(const mint& lhs, const mint& rhs){return mr.fix(lhs._v) == mr.fix(rhs._v);}
  friend bool operator!=(const mint& lhs, const mint& rhs){return mr.fix(lhs._v) != mr.fix(rhs._v);}
private:
  unsigned long long _v;
  static montgomery_reduction_64bit mr;
  static unsigned long long umod(){return mr.mod();}
};

template <int id>
montgomery_reduction_64bit dynamic_modint64bit<id>::mr(998244353);
#endif