#ifndef _BITSET2D_H_
#define _BITSET2D_H_
#include <vector>
#include <iostream>
#include <cassert>
struct bitset2d{
  int n, m;
private:
  using mat = std::vector<std::vector<bool>>;
  mat v;
public:
  bitset2d(): n(0), m(0){}
  bitset2d(int _n, int _m, bool _x = 0): n(_n), m(_m), v(n, std::vector<bool>(m, _x)){}
  
  std::pair<int, int> size(){
    return {n, m};
  }
  int popcount(){
    int res = 0;
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        res += v[i][j];
      }
    }
    return res;
  }
  // res[i][j] := (元の行列)[j][i]
  bitset2d t() const {
    bitset2d res(m, n);
    for(int i = 0; i < m; i++){
      for(int j = 0; j < n; j++){
        res[i][j] = v[j][i];
      }
    }
    return res;
  }
  // 時計回りに90度回転
  // res[i][j] = (元の行列)[n - 1 - j][i]
  bitset2d rot() const {
    bitset2d res(m, n);
    for(int i = 0; i < m; i++){
      for(int j = 0; j < n; j++){
        res[i][j] = v[n - 1 - j][i];
      }
    }
    return res;
  }
  // 右にdx, 下にdyずらす. (0, 0)が(dx, dy)になる
  // ずれた分は0が初期値
  // @param dx, dy >= 0
  bitset2d shift(int dx, int dy) const {
    assert(0 <= dx && 0 <= dy);
    bitset2d res(n + dx, m + dy);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        res[i + dx][j + dy] = v[i][j];
      }
    }
    return res;
  }
  // 元の行列の[0, rx) × [0, ry)の部分
  // rx < nの場合, はみ出た部分が消える
  // n < rxの場合, 足りない部分を0埋めする
  // yについても同様
  // 必ずrx × ryの行列を返す
  // @param rx, ry >= 0
  bitset2d submat(int rx, int ry) const {
    assert(0 <= rx && 0 <= ry);
    bitset2d res(rx, ry);
    for(int i = 0; i < std::min(n, rx); i++){
      for(int j = 0; j < std::min(m, ry); j++){
        res[i][j] = v[i][j];
      }
    }
    return res;
  }
  // 元の行列の[lx, rx) × [ly, ry)の部分
  // rx < nの場合, はみ出た部分が消える
  // n < rxの場合, 足りない部分を0埋めする
  // yについても同様
  // 必ず(rx - lx) × (ry - ly)の行列を返す
  // @param 0 <= lx <= rx, 0 <= ly <= ry
  bitset2d submat(int lx, int rx, int ly, int ry) const {
    assert(0 <= lx && lx <= rx);
    assert(0 <= ly && ly <= ry);
    bitset2d res(rx - lx, ry - ly);
    for(int i = lx; i < std::min(n, rx); i++){
      for(int j = ly; j < std::min(m, ry); j++){
        res[i - lx][j - ly] = v[i][j];
      }
    }
    return res;
  }
  // max(a.n, b.n) × max(a.m, b.m)
  // bitwise_or
  static bitset2d or_mat(const bitset2d &a, const bitset2d &b){
    bitset2d res(std::max(a.n, b.n), std::max(a.m, b.m));
    for(int i = 0; i < a.n; i++){
      for(int j = 0; j < a.m; j++){
        res[i][j] = a.v[i][j];
      }
    }
    for(int i = 0; i < b.n; i++){
      for(int j = 0; j < b.m; j++){
        res[i][j] = res[i][j] || b.v[i][j];
      }
    }
    return res;
  }
  // min(a.n, b.n) × min(a.m, b.m)
  // bitwise_and
  static bitset2d and_mat(const bitset2d &a, const bitset2d &b){
    bitset2d res(std::min(a.n, b.n), std::min(a.m, b.m));
    for(int i = 0; i < std::min(a.n, b.n); i++){
      for(int j = 0; j < std::min(a.m, b.m); j++){
        res[i][j] = a.v[i][j] && b.v[i][j];
      }
    }
    return res;
  }
  // max(a.n, b.n) × max(a.m, b.m)
  // bitwise_xor
  static bitset2d xor_mat(const bitset2d &a, const bitset2d &b){
    bitset2d res(std::max(a.n, b.n), std::max(a.m, b.m));
    for(int i = 0; i < a.n; i++){
      for(int j = 0; j < a.m; j++){
        res[i][j] = a.v[i][j];
      }
    }
    for(int i = 0; i < b.n; i++){
      for(int j = 0; j < b.m; j++){
        res[i][j] = res[i][j] ^ b.v[i][j];
      }
    }
    return res;
  }
  bitset2d operator | (const bitset2d &b) const {return or_mat(*this, b);}
  bitset2d operator & (const bitset2d &b) const {return and_mat(*this, b);}
  bitset2d operator ^ (const bitset2d &b) const {return xor_mat(*this, b);}
  std::vector<bool> &operator [](int i){return v[i];}
  void print(){
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++) std::cout << v[i][j];
      std::cout << '\n';
    }
  }
};
#endif