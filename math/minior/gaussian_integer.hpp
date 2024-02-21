#ifndef _GAUSSIAN_INTEGER_H_
#define _GAUSSIAN_INTEGER_H_
#include <algorithm>
#include <cassert>
template<typename Int>
struct gaussian_integer{
  using gint = gaussian_integer<Int>;
  Int a, b;
  gaussian_integer(Int _a, Int _b): a(_a), b(_b){}
  Int norm(){return a * a + b * b;}
  gint operator *= (gint r){
    Int x = a * r.a - b * r.b;
    Int y = a * r.b + b * r.a;
    a = x;
    b = y;
    return *this;
  }
  gint operator /= (const gint r){
    (*this) *= gint(r.a, -r.b);
    Int q = r.a * r.a + r.b * r.b;
    a = 2 * a + q, b = 2 * b + q;
    q *= 2;
    a = (a / q) - (a % q < 0);
    b = (b / q) - (b % q < 0);
    return *this;
  }
  gint operator += (const gint r){
    a += r.a;
    b += r.b;
    return *this;
  }
  gint operator -= (const gint r){
    a -= r.a;
    b -= r.b;
    return *this;
  }
  gint operator %= (const gint r){
    gint q = (*this) / r;
    (*this) -= q * r;
    return *this;
  }
  gint operator * (const gint r)const{
    gint res(*this);
    return res *= r;
  }
  gint operator / (const gint r)const{
    gint res(*this);
    return res /= r;
  }
  gint operator % (const gint r)const{
    gint res(*this);
    return res %= r;
  }
  gint operator + (const gint r)const{
    gint res(*this);
    return res += r;
  }
  gint operator - (const gint r)const{
    gint res(*this);
    return res -= r;
  }
  static gint gcd(gint a, gint b){
    if(a.norm() == 0) return b;
    if(b.norm() == 0) return a;
    if(a.norm() < b.norm()) std::swap(a, b);
    while(b.norm()){
      a %= b;
      std::swap(a, b);
    }
    return a;
  }
};
#endif