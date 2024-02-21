#ifndef _MATRIX_MOD2_H_
#define _MATRIX_MOD2_H_
#include <vector>
#include <iostream>
#include "../../data_structure/bit_sequence/dynamic_bitset.hpp"

struct matrix_mod2{
  int n, m;
  using ull = unsigned long long;
  static constexpr int bitlen = 64;
  static constexpr int bitlen_shift = 6;
  static constexpr int bitlen_mod = 63;

private:
  using vec = std::vector<ull>;
  using matrix = matrix_mod2;

  // n × k 行列と k × m 行列の積(n × m行列)
  // K == 0だと壊れる
  // O(NKM / w)
  static matrix __mul_mat(const matrix &vl, const matrix &vr){
    int N = vl.n, K = vl.m, M = vr.m;
    assert(K == vr.n);
    assert(K);
    if(N == 0) return matrix(0, M, 0);
    if(M == 0) return matrix(N, 0, 0);
    auto vr_t = vr.t();
    matrix ret(N, M, 0);
    for(int i = 0; i < N; i++){
      for(int j = 0; j < M; j++){
        ull S = 0;
        for(int k = 0; k < ((K + bitlen_mod) >> bitlen_shift); k++){
          S ^= vl.val[i][k] & vr_t.val[j][k];
        }
        ret.val[i][j >> bitlen_shift] |= (ull)(__builtin_popcountll(S) & 1) << (j & bitlen_mod);
      }
    }
    return ret;
  }
  // n × m 行列と n × m 行列の和(n × m行列)
  static void __add_mat_inplace(matrix &vl, const matrix &vr){
    assert(vl.n == vr.n && vl.m == vr.m);
    int N = vl.n, M = vl.m;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < ((M + bitlen_mod) >> bitlen_shift) ; j++){
        vl.val[i][j] ^= vr.val[i][j];
      }
    }
  }
  std::vector<vec> val;
public:
  matrix_mod2(): n(0), m(0){}
  matrix_mod2(int _n, int _m, bool f = 0) : n(_n), m(_m), val(_n, vec((_m + bitlen_mod) >> bitlen_shift, 0)){
    if(f){
      int __m = (m + bitlen_mod) >> bitlen_shift;
      int last = m - (__m << bitlen_shift) + bitlen;
      ull last_bit = (last == bitlen ? ~(ull)0 : ((ull)1 << last) - 1);
      for(int i = 0; i < n; i++){
        for(int j = 0; j < __m; j++){
          if(j != __m - 1) val[i][j] = ~(ull)0;
          else val[i][j] = last_bit;
        }
      }
    }
  }
  matrix_mod2(int m, const dynamic_bitset &v): n(1), m(m), val(1, vec((m + bitlen_mod) >> bitlen_shift, 0)){
    assert(m == v.n);
    for(int i = 0; i < v.im; i++) val[0][i] = v.v[i];
  }
  matrix_mod2(const matrix_mod2 &v) : n(v.n), m(v.m), val(v.val){}
  matrix_mod2 operator +  (const matrix_mod2 &vr){matrix_mod2 tmp(*this); return tmp += vr;}
  matrix_mod2 operator *  (const matrix_mod2 &vr){return __mul_mat(*this, vr);}
  matrix_mod2 operator ^  (const long long vr){return pow(vr);}
  matrix_mod2 operator += (const matrix_mod2 &vr){__add_mat_inplace(*this, vr); return *this;}
  matrix_mod2 operator *= (const matrix_mod2 &vr){return (*this) = __mul_mat(*this, vr);}
  matrix_mod2 operator ^= (const long long vr){return (*this) = pow(vr);}

  void set(int i, int j, bool f){
    bool g = (val[i][j >> bitlen_shift] >> (j & bitlen_mod)) & 1;
    if(f != g) val[i][j >> bitlen_shift] ^= (ull)1 << (j & bitlen_mod);
  }
  bool get(int i, int j){
    return (val[i][j >> bitlen_shift] >> (j & bitlen_mod)) & 1;
  }
  // n次の単位行列
  static matrix_mod2 eye(int n){
    matrix_mod2 ret(n, n, 0);
    for(int i = 0; i < n; i++) ret.val[i][i >> bitlen_shift] = (ull)1 << (i & bitlen_mod);
    return ret;
  }
  void print()const{
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        int f = (val[i][j >> bitlen_shift] >> (j & bitlen_mod)) & 1;
        std::cout << f << (j == m - 1 ? '\n' : ' ');
      }
    }
  }
  // O(n^3 log k / w)
  matrix_mod2 pow(long long k){
    assert(n && m && n == m); // 正方行列でなければならない
    matrix_mod2 ret = eye(n); // k == 0の場合単位行列を返す
    matrix_mod2 m(*this);
    while(k){
      if(k & 1) ret *= m;
      m *= m;
      k >>= 1;
    }
    return ret;
  }
  // 転置, O(nm)
  matrix_mod2 t()const{
    matrix_mod2 ret(m, n, 0);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        bool f = (val[i][j >> bitlen_shift] >> (j & bitlen_mod)) & 1; // (i, j)のビット
        ret.val[j][i >> bitlen_shift] |= (ull)f << (i & bitlen_mod); // (j, i)にセット　
      }
    }
    return ret;
  }
  //掃き出し法で上三角行列を作る, {変形後の行列、行のスワップ回数}を返す O(NM^2 / w)
  std::pair<matrix_mod2, int> gaussian_elimination(){
    matrix_mod2 v(*this);
    int row = 0, swp = 0;
    for(int i = 0; i < m && row < n; i++){
      //i列目が0でない行を探す
      int r = -1;
      for(int j = row; j < n; j++){
        if((v.val[j][i >> bitlen_shift] >> (i & bitlen_mod)) & 1){
          r = j;
          break;
        }
      }
      if(r == -1) continue;
      if(r != row){
        swp++;
        std::swap(v.val[r], v.val[row]);
      }
      //i列目が0でない行の処理
      for(int j = row + 1; j < n; j++){
        if(((v.val[j][i >> bitlen_shift] >> (i & bitlen_mod)) & 1) == 0) continue;
        for(int k = (i >> bitlen_shift); k < ((m + bitlen_mod) >> bitlen_shift); k++){
          v.val[j][k] ^= v.val[row][k];
        }
      }
      row++;
    }
    return {v, swp};
  }
  
  //すでに上三角行列になっていることが前提
  int rank(){
    int cnt = 0;
    for(int i = 0; i < n; i++, cnt++){
      bool f = false;
      for(int j = (i >> bitlen_shift); j < ((m + bitlen_mod) >> bitlen_shift); j++){
        if(val[i][j]){
          f = true;
          break;
        }
      }
      if(!f) break;
    }
    return cnt;
  }
  // 行列式 O(N^3 / w)
  bool det(){
    assert(n == m); // 正方行列のみ
    auto [tmp, swp] = gaussian_elimination();
    for(int i = 0; i < n; i++) if(((tmp.val[i][i >> bitlen_shift] >> (i & bitlen_mod)) & 1) == 0) return false;
    return true;
  }
  // (n, m) + (l, m) -> (n + l, m) 縦に結合
  matrix_mod2 concat_vertical(matrix_mod2 vr){
    assert(m == vr.m);
    matrix_mod2 ret(*this);
    for(int i = 0; i < vr.n; i++) ret.val.push_back(vr.val[i]);
    ret.n += vr.n;
    return ret;
  }
  // (n, m) + (n, l) -> (n, m + l)　横に結合
  matrix_mod2 concat_horizontal(matrix_mod2 vr){
    assert(n == vr.n);
    matrix_mod2 ret(*this);

    if((m & bitlen_mod) == 0){
      for(int i = 0; i < n; i++){
        ret.val[i].insert(ret.val[i].end(), vr.val[i].begin(), vr.val[i].end());
      }
    }else{
      int s1 = (m & bitlen_mod), s2 = bitlen - s1;
      int M = (m + bitlen_mod) >> bitlen_shift;
      int M2 = (vr.m + bitlen_mod) >> bitlen_shift;
      int K = (m + vr.m + bitlen_mod) >> bitlen_shift;
      for(int i = 0; i < n; i++){
        ret.val[i].resize(K, 0);
        for(int j = M - 1, k = 0; k <= M2; j++, k++){
          if(k < M2) ret.val[i][j] |= vr.val[i][k] << s1;
          if(k) ret.val[i][j] |= vr.val[i][k - 1] >> s2;
        }
      }
    }
    ret.m += vr.m;
    return ret;
  }
  // (n, m) -> (k, m), (n - k, m)
  std::pair<matrix_mod2, matrix_mod2> split_vertical(int k){
    assert(0 <= k && k <= n);
    matrix_mod2 a(k, m, 0), b(n - k, m, 0);
    for(int i = 0; i < k; i++){
      for(int j = 0; j < ((m + bitlen_mod) >> bitlen_shift); j++){
        a.val[i][j] = val[i][j];
      }
    }
    for(int i = 0; i < n - k; i++){
      for(int j = 0; j < ((m + bitlen_mod) >> bitlen_shift); j++){
        b.val[i][j] = val[k + i][j];
      }
    }
    return {a, b};
  }
  // (n, m) -> (n, k), (n, m - k)
  std::pair<matrix_mod2, matrix_mod2> split_horizontal(int k){
    assert(0 <= k && k <= m);
    matrix_mod2 a(n, k, 0), b(n, m - k, 0);
    int K = (k + bitlen_mod) >> bitlen_shift;
    for(int i = 0; i < n; i++){
      std::copy(val[i].begin(), val[i].begin() + K, a.val[i].begin());
    }
    if((k & bitlen_mod) == 0){
      for(int i = 0; i < n; i++){
        std::copy(val[i].begin() + K, val[i].end(), b.val[i].begin());
      }
      return {a, b};
    }
    int s1 = k & bitlen_mod, s2 = bitlen - s1;
    for(int i = 0; i < n; i++) a.val[i].back() &= ((ull)1 << s1) - 1;
    int M = (m - k + bitlen_mod) >> bitlen_shift;
    for(int i = 0; i < n; i++){
      for(int j = (k >> bitlen_shift), t = 0; t < M; j++, t++){
        b.val[i][t] = val[i][j] >> s1;
        if(j + 1 < val[i].size()) b.val[i][t] |= val[i][j + 1] << s2;
      }
    }
    return {a, b};
  }
  // O(n^3 / w)
  matrix_mod2 inv(){
    assert(n == m);
    auto [tmp, swp] = concat_horizontal(eye(n)).gaussian_elimination();
    for(int i = 0; i < n; i++){
      if(tmp.get(i, i) == 0){
        return matrix_mod2{};
      }
    }
    for(int i = n - 1; i >= 0; i--){
      for(int j = i + 1; j < n; j++){
        if(tmp.get(i, j) == 0) continue;
        for(int k = (j >> bitlen_shift); k < tmp.val[0].size(); k++){
          tmp.val[i][k] ^= tmp.val[j][k];
        }
      }
    }
    return tmp.split_horizontal(n).second;
  }
  /*
  // https://ja.wikipedia.org/wiki/LU%E5%88%86%E8%A7%A3
  std::pair<matrix_mod2, matrix_mod2> lu_decomposition(){
    matrix_mod2 l = eye(n), u(n, n, 0);
    for(int i = 0; i < n; i++){
      // u[i][i]を決定
      u[i][i] = val[i][i];
      for(int j = 0; j < i; j++) u[i][i] -= l[i][j] * u[i][j];
      if(u[i][i].val() == 0) return {matrix_mod{}, matrix_mod{}}; // 不可能
      mint iuii = u[i][i].inv();

      // l[0, n)[i]を決定
      for(int j = i + 1; j < n; j++){
        l[j][i] = val[j][i];
        for(int k = 0; k < i; k++) l[j][i] -= l[j][k] * u[i][k];
        l[j][i] *= iuii;
      }
      // u[i][0, n)を決定
      for(int j = i + 1; j < n; j++){
        u[j][i] = val[i][j];
        for(int k = 0; k < i; k++) u[j][i] -= l[i][k] * u[j][k];
      }
    }
    u = u.t();
    return {l, u};
  }
  */
  // Ax = b
  // (n, m) * (m, 1) -> (n, 1)
  // を満たす連立方程式を解く
  // {解空間の次元, 解の1つ}
  // ない場合は長さ0のbitset(ある場合は長さmのbitset)
  std::pair<int, dynamic_bitset> system_of_linear_equations(const dynamic_bitset &vr){
    assert(vr.n == n);
    matrix_mod2 tmp = concat_horizontal(matrix_mod2(n, vr).t()).gaussian_elimination().first;
    //解空間の次元 = 変数の数 - 階数
    int r = tmp.rank();
    std::vector<int> fc(r, -1);//各行に初めて非零要素が現れる列
    for(int i = 0; i < r; i++){
      bool f = false;
      for(int j = i; j < tmp.m; j++){
        if(tmp.get(i, j) == 0) continue;
        if(j == tmp.m - 1 && !f){
          return {-1, dynamic_bitset(0)}; // 解なし
        }
        if(!f){
          fc[i] = j;
          f = true;
        }
      }
    }
    int d = tmp.m - 1 - r, v = tmp.m - 1;
    dynamic_bitset plus(v, 0);
    for(int i = r - 1; i >= 0; i--){
      int idx = fc[i];
      assert(idx != -1);
      plus.set(idx, tmp.get(i, v));
      for(int j = idx + 1; j < v; j++){
        bool f = plus.get(j) && tmp.get(i, j);
        bool g = plus.get(idx);
        plus.set(idx, f ^ g);
      }
    }
    return {d, plus};
  }
};
#endif
