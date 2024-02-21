#ifndef _BINARY_BASIS_H_
#define _BINARY_BASIS_H_

#include "../../data_structure/bit_sequence/dynamic_bitset.hpp"
#include <algorithm>
#include <vector>

template<typename T, int BITLEN>
struct binary_basis{
  int r, zero;
  T val[BITLEN];
  binary_basis(): r(0), zero(0){}
  binary_basis(T x): r(0), zero(0){
    if(x) val[r++] = x;
    else zero++;
  }
  // 要素数
  int size(){
    return r + zero;
  }
  // ランク
  int rank(){
    return r;
  }
  // マージ O(BITLEN ^ 2)
  void merge(const binary_basis<T, BITLEN> &v){
    for(int i = 0; i < v.r; i++) push(v.val[i]);
  }
  // xを追加, O(BITLEN)
  void push(T x){
    for(int j = 0; j < r; j++) x = std::min(x, x ^ val[j]);
    if(x){
      val[r++] = x;
      for(int j = r - 1; j && val[j] > val[j - 1]; j--) std::swap(val[j], val[j - 1]);
    }else zero++;
    assert(r <= BITLEN);
  }
  // xを追加, O(BITLEN)
  // msbがiの要素がある -> 他の全ての要素のiのビットが0になる　　ように変形
  void push2(T x){
    for(int j = 0; j < r; j++) x = std::min(x, x ^ val[j]);
    if(x){
      val[r++] = x;
      int j = r - 1;
      for(; j && val[j] > val[j - 1]; j--) std::swap(val[j], val[j - 1]);
      for(int i = j + 1; i < r; i++) val[j] = std::min(val[j], val[j] ^ val[i]);
      for(int i = j - 1; i >= 0; i--) val[i] = std::min(val[i], val[i] ^ val[j]);
    }else zero++;
    assert(r <= BITLEN);
  }
  // 0個以上の要素のxorで表せる最大値, O(BITLEN)
  T max_xor(){
    T ans = 0;
    for(int i = 0; i < r; i++) ans = std::max(ans, ans ^ val[i]);
    return ans;
  }
  // 1個以上の要素のxorで表せる最小値, O(1)
  T min_xor(){
    assert(r);
    return val[r - 1];
  }
  // 0個以上のxorでxを作れるか
  bool equation(T x){
    for(int i = 0; i < r; i++) x = std::min(x, x ^ val[i]);
    return x == 0;
  }
};

template<typename T, int BITLEN>
struct binary_basis_constructive{
  int r, zero;
  T val[BITLEN], t[BITLEN], use[BITLEN];
  binary_basis_constructive(): r(0), zero(0){}
  binary_basis_constructive(T x): r(0), zero(0){
    if(x) val[r] = t[r] = x, use[r++] = 1;
    else zero++;
  }
  // 要素数
  int size(){
    return r + zero;
  }
  // ランク
  int rank(){
    return r;
  }
  // マージ O(BITLEN ^ 2)
  void merge(const binary_basis_constructive<T, BITLEN> &v){
    for(int i = 0; i < v.r; i++) push(v.val[i]);
  }
  // xを追加, O(BITLEN)
  void push(T x){
    T y = 0, z = x;
    for(int j = 0; j < r; j++){
      if(x > (x ^ val[j])){
        x ^= val[j];
        y ^= use[j];
      }
    }
    if(x){
      val[r] = x;
      use[r] = y ^ ((T)1 << r);
      t[r++] = z;
      for(int j = r - 1; j && val[j] > val[j - 1]; j--){
        std::swap(val[j], val[j - 1]);
        std::swap(use[j], use[j - 1]);
      }
    }else zero++;
    assert(r <= BITLEN);
  }
  // 0個以上の要素のxorで表せる最大値, O(BITLEN)
  T max_xor(){
    T ans = 0;
    for(int i = 0; i < r; i++) ans = std::max(ans, ans ^ val[i]);
    return ans;
  }
  // 1個以上の要素のxorで表せる最小値, O(1)
  T min_xor(){
    assert(r);
    return val[r - 1];
  }
  // 0個以上のxorでxを作れるか
  bool equation(T x){
    for(int i = 0; i < r; i++) x = std::min(x, x ^ val[i]);
    return x == 0;
  }
  // 0個以上の要素のxorで表せる最大値, その時に使う数(k-bit目が1の時, t[k]を使う), O(BITLEN)
  std::pair<T, T> max_xor_construct(){
    T ans = 0, argans = 0;
    for(int i = 0; i < r; i++){
      if(ans < (ans ^ val[i])){
        ans ^= val[i];
        argans ^= use[i];
      }
    }
    return {ans, argans};
  }
  // 1個以上の要素のxorで表せる最小値, その時に使う数(k-bit目が1の時, t[k]を使う), O(1)
  std::pair<T, T> min_xor_construct(){
    assert(r);
    return {val[r - 1], use[r - 1]};
  }
  // 0個以上の要素のxorでxを作れるか, その時に使う数(k-bit目が1の時, t[k]を使う), O(BITLEN)
  std::pair<bool, T> equation_construct(T x){
    T y = 0;
    for(int i = 0; i < r; i++){
      if(x > (x ^ val[i])){
        x ^= val[i];
        y ^= use[i];
      }
    }
    if(x == 0) return {true, y};
    return {false, 0};
  }
};

// ビット長が長い場合
struct binary_basis_long{
  int r, zero, m;
  // 0 : aのmsbの方が大きい
  // 1 : msbが等しい
  // 2 : bのmsbの方が大きい
  // O(m / w)
  int compare_msb(dynamic_bitset &a, dynamic_bitset &b){
    int i = a.im - 1, j = b.im - 1;
    while(i >= 0 || j >= 0){
      if(i == j){
        if(a.v[i] && b.v[i]){
          int A = find_prev_64bit(a.v[i], 63);
          int B = find_prev_64bit(b.v[i], 63);
          return (A > B ? 0 : A == B ? 1 : 2);
        }else if(a.v[i]){
          // aだけ正
          return 0;
        }else if(b.v[i]){
          // bだけ正
          return 2;
        }else{
          // 両方0
          i--, j--;
          continue;
        }
      }else if(i > j){
        if(a.v[i--] == 0) continue;
        return 0;
      }else{
        if(b.v[j--] == 0) continue;
        return 2;
      }
    }
    return 1; // a == 0 && b == 0
  }
  // O(m / w)
  int get_msb(dynamic_bitset &a){
    for(int i = a.im - 1; i >= 0; i--){
      if(a.v[i] == 0) continue;
      return (i * a.BITLEN) + find_prev_64bit(a.v[i], 63);
    }
    return -1; // a == 0
  }
  std::vector<dynamic_bitset> val;
  binary_basis_long(int m): r(0), zero(0), m(m){}
  binary_basis_long(int m, dynamic_bitset &x): r(0), zero(0), m(m){
    assert(m == x.n);
    if(x.any()) r++, val.push_back(x);
    else zero++;
  }
  // 要素数
  int size(){
    return r + zero;
  }
  // ランク
  int rank(){
    return r;
  }
  // マージ O(BITLEN ^ 2)
  void merge(const binary_basis_long &v){
    assert(m == v.m);
    for(int i = 0; i < v.r; i++) push(v.val[i]);
  }
  // xを追加, O(m ^ 2 / w)
  void push(dynamic_bitset x){
    for(int i = 0; i < r; i++){
      int f = compare_msb(x, val[i]);
      if(f == 0){
        r++;
        assert(r <= m);
        val.insert(val.begin() + i, x);
        return;
      }else if(f == 1){
        x ^= val[i];
      }
    }
    if(x.any()){
      r++;
      val.push_back(x);
      assert(r <= m);
    }else zero++;
  }
  // 0個以上の要素のxorで表せる最大値, O(m ^ 2 / w)
  dynamic_bitset max_xor(){
    dynamic_bitset ans(m, 0);
    for(int i = 0; i < r; i++){
      int msb = get_msb(val[i]);
      if(ans.get(msb)) continue;
      ans ^= val[i];
    }
    return ans;
  }
  // 1個以上の要素のxorで表せる最小値, O(m / w)
  dynamic_bitset min_xor(){
    assert(r);
    return val[r - 1];
  }
  // 0個以上のxorでxを作れるか, O(m ^ 2 / w)
  bool equation(dynamic_bitset x){
    for(int i = 0; i < r; i++){
      int f = compare_msb(x, val[i]);
      if(f == 0){
        return false;
      }else if(f == 1){
        x ^= val[i];
      }
    }
    return x.none();
  }
};
#endif