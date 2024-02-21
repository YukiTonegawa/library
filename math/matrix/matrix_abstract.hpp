#ifndef _MATRIX_ABSTRACT_H_
#define _MATRIX_ABSTRACT_H_
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>
#include "../../algebraic_structure/semi_ring.hpp"

template<typename semi_ring>
struct matrix_abstract{
  using Val = typename semi_ring::Val;
  static constexpr auto id_add = semi_ring::id_add;
  static constexpr auto id_mul = semi_ring::id_mul;
  static constexpr auto add = semi_ring::add;
  static constexpr auto mul = semi_ring::mul;

  int n, m;
  using vec = std::vector<Val>;
  using mat = std::vector<vec>;
  using mat_ab = matrix_abstract<semi_ring>;
  mat val;
  static mat_ab _mul(const mat_ab &vl, const mat_ab &vr){
    assert(vl.n && vr.n);
    int N = vl.n, K = vl.m, M = vr.m;
    assert(K == vr.n);
    matrix_abstract ret(N, M);
    auto vr_t = vr.t();
    for(int i = 0; i < N; i++){
      for(int j = 0; j < M; j++){
        for(int k = 0; k < K; k++){
          ret[i][j] = add(ret[i][j], mul(vl.val[i][k], vr_t[j][k]));
        }
      }
    }
    return ret;
  }
  matrix_abstract(): n(0), m(0){}
  matrix_abstract(int _n, int _m) : n(_n), m(_m), val(_n, vec(_m, id_add())){}
  matrix_abstract(int _n, int _m, Val x) : n(_n), m(_m), val(_n, vec(_m, x)){}
  matrix_abstract(const mat_ab &v) : n(v.n), m(v.m), val(v.val){}
  matrix_abstract(const vec &v): n(1), m(v.size()), val(1, vec(v.size())){val[0] = v;}
  matrix_abstract(const mat &v): n(v.size()), m(v[0].size()), val(v){}

  mat_ab operator *  (const mat_ab  &vr){return           _mul(*this, vr);}
  mat_ab operator *= (const mat_ab  &vr){return (*this) = _mul(*this, vr);}
  vec&    operator [] (const int i){return val[i];}

  // n次の単位行列
  static mat_ab eye(int n){
    mat_ab ret(n, n, id_add());
    for(int i = 0; i < n; i++) ret[i][i] = id_mul();
    return ret;
  }
  void print(){
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        std::cout << val[i][j] << (j == m - 1 ? '\n' : ' ');
      }
    }
  }
  mat_ab pow(long long k){
    assert(n && m && n == m); // 正方行列でなければならない
    mat_ab ret = eye(n); // k == 0の場合単位行列を返す
    mat_ab mul_mat(*this);
    while(k){
      if(k & 1) ret *= mul_mat;
      k >>= 1;
      mul_mat *= mul_mat;
    }
    return ret;
  }
  // 転置
  mat_ab t()const{
    mat_ab ret(m, n);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        ret[j][i] = val[i][j];
      }
    }
    return ret;
  }
};

// 静的サイズのN * N正方行列
template<int N, typename semi_ring>
struct square_matrix_abstract{
  using Val = typename semi_ring::Val;
  static constexpr auto id_add = semi_ring::id_add;
  static constexpr auto id_mul = semi_ring::id_mul;
  static constexpr auto add = semi_ring::add;
  static constexpr auto mul = semi_ring::mul;

  using vec = std::array<Val, N>;
  using mat = std::array<vec, N>;
  using mat_ab = square_matrix_abstract<N, semi_ring>;
  mat val;
  static mat_ab _mul(const mat_ab &vl, const mat_ab &vr){
    mat_ab ret;
    auto vr_t = vr.t();
    for(int i = 0; i < N; i++){
      for(int j = 0; j < N; j++){
        for(int k = 0; k < N; k++){
          ret[i][j] = add(ret[i][j], mul(vl.val[i][k], vr_t[j][k]));
        }
      }
    }
    return ret;
  }
  square_matrix_abstract(){for(int i = 0; i < N; i++) val[i].fill(id_add());}
  square_matrix_abstract(const mat_ab &v): val(v.val){}

  bool operator == (const mat_ab &vr)const{
    for(int i = 0; i < N; i++) if(val[i] != vr.val[i]) return false;
    return true;
  }
  bool operator != (const mat_ab &vr)const{return !((*this) == vr);}
  mat_ab operator *  (const mat_ab  &vr){return           _mul(*this, vr);}
  mat_ab operator *= (const mat_ab  &vr){return (*this) = _mul(*this, vr);}
  vec&    operator [] (const int i){return val[i];}

  // n次の単位行列
  static mat_ab eye(){
    mat_ab ret;
    for(int i = 0; i < N; i++) ret[i][i] = id_mul();
    return ret;
  }
  void print(){
    for(int i = 0; i < N; i++){
      for(int j = 0; j < N; j++){
        std::cout << val[i][j] << (j == N - 1 ? '\n' : ' ');
      }
    }
  }
  mat_ab pow(long long k){
    mat_ab ret = eye(); // k == 0の場合単位行列を返す
    mat_ab mul_mat(*this);
    while(k){
      if(k & 1) ret *= mul_mat;
      k >>= 1;
      mul_mat *= mul_mat;
    }
    return ret;
  }
  // 転置
  mat_ab t()const{
    mat_ab ret;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < N; j++){
        ret[j][i] = val[i][j];
      }
    }
    return ret;
  }
};
#endif