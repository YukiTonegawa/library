#ifndef _SPARSE_TABLE2D_H_
#define _SPARSE_TABLE2D_H_
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <limits>
#include "sparse_table.hpp"

template<typename T, T (*merge)(T, T), T (*id)()>
struct sparse_table2d{
  std::vector<std::vector<std::vector<T>>> table;
  int n, m, m2;
  sparse_table2d(): n(0), m(0){}
  sparse_table2d(const std::vector<std::vector<T>> &v): n(v.size()), m(0){
    if(n == 0) return;
    int n2 = 0;
    m = v[0].size(), m2 = 0;
    while((1 << n2) <= n) n2++;
    while((1 << m2) <= m) m2++;
    table.resize(n2 * m2, std::vector<std::vector<T>>(n, std::vector<T>(m)));
    table[0] = v;
    // y
    for(int i = 0; i < n2; i++){
      for(int j = 1; j < m2; j++){
        int yshift = 1 << (j - 1);
        int zidx = i * m2 + j;
        for(int x = 0; x < n; x++){
          for(int y = 0; y + yshift < m; y++){
            table[zidx][x][y] = merge(table[zidx - 1][x][y], table[zidx - 1][x][y + yshift]);
          }
        }
      }
    }
    // x
    for(int j = 0; j < m2; j++){
      for(int i = 1; i < n2; i++){
        int xshift = 1 << (i - 1);
        int zidx = i * m2 + j;
        for(int x = 0; x + xshift < n; x++){
          for(int y = 0; y < m; y++){
            table[zidx][x][y] = merge(table[zidx - m2][x][y], table[zidx - m2][x + xshift][y]);
          }
        }
      }
    }
  }
  // 最左bit
  int msb(int x){
    return (x == 1 ? 0 : 31 - __builtin_clz(x - 1));
  }
  T query(int lx, int rx, int ly, int ry){
    lx = std::max(lx, 0);
    rx = std::min(rx, n);
    ly = std::max(ly, 0);
    ry = std::min(ry, m);
    if(lx >= rx || ly >= ry) return id();
    int len_x = rx - lx, len_y = ry - ly;
    int bx = msb(len_x), by = msb(len_y), bz = bx * m2 + by;
    T lower = merge(table[bz][lx][ly], table[bz][lx][ry - (1 << by)]);
    T upper = merge(table[bz][rx - (1 << bx)][ly], table[bz][rx - (1 << bx)][ry - (1 << by)]);
    return merge(lower, upper);
  }
};

template<typename T>
struct rmq2d{
private:
  static constexpr T min_func(T a, T b){
    return std::min(a, b);
  }
  static constexpr T id(){
    return std::numeric_limits<T>::max();
  }
  static constexpr int block_size = 4;
  static constexpr int block_size_log = 2;
  sparse_table2d<T, min_func, id> st;
  std::vector<std::vector<rmq<T>>> row, col; // 2 * n, 2 * m 個のrmq
public:
  int n, m;
  rmq2d(): n(0), m(0){}
  rmq2d(const std::vector<std::vector<T>> &v): n(v.size()), m(0){
    if(v.empty()) return;
    m = v[0].size();
    int n2 = (n + block_size - 1) / block_size, m2 = (m + block_size - 1) / block_size;
    std::vector<std::vector<T>> v2(n2, std::vector<T>(m2));
    for(int i = 0; i < n2; i++){
      for(int j = 0; j < m2; j++){
        v2[i][j] = v[i * block_size][j * block_size];
        for(int x = i * block_size, k = 0; x < n && k < block_size; x++, k++){
          for(int y = j * block_size, t = 0; y < m && t < block_size; y++, t++){
            v2[i][j] = min_func(v2[i][j], v[x][y]);
          }
        }
      }
    }
    st = sparse_table2d<T, min_func, id>(v2);
    v2 = v;
    row.resize(block_size_log, std::vector<rmq<T>>(n));
    for(int k = 0; k < block_size_log; k++){
      for(int i = 0; i < n; i++) row[k][i] = rmq<T>(v[i]);
      if(k + 1 != block_size_log){
        for(int i = 0; i + (1 << k) < n; i++){
          for(int j = 0; j < m; j++){
            v2[i][j] = std::min(v2[i][j], v2[i + (1 << k)][j]);
          }
        }
      }
    }
    col.resize(block_size_log, std::vector<rmq<T>>(m));
    v2.resize(m, std::vector<T>(n));
    for(int i = 0; i < n; i++) for(int j = 0; j < m; j++) v2[j][i] = v[i][j];
    for(int k = 0; k < block_size_log; k++){
      for(int j = 0; j < m; j++) col[k][j] = rmq<T>(v2[j]);
      if(k + 1 != block_size_log){
        for(int j = 0; j + (1 << k) < m; j++){
          for(int i = 0; i < m; i++){
            v2[j][i] = std::min(v2[j][i], v2[j + (1 << k)][i]);
          }
        }
      }
    }
  }
  T query(int lx, int rx, int ly, int ry){
    lx = std::max(lx, 0);
    rx = std::min(rx, n);
    ly = std::max(ly, 0);
    ry = std::min(ry, m);
    if(lx >= rx || ly >= ry) return id();
    int lxb = (lx + block_size - 1) / block_size, rxb = rx / block_size;
    int lyb = (ly + block_size - 1) / block_size, ryb = ry / block_size;
    int lxceil = lxb * block_size, rxfloor = rxb * block_size;
    int lyceil = lyb * block_size, ryfloor = ryb * block_size;
    T ret = st.query(lxb, rxb, lyb, ryb); 
    if(lxb > rxb){
      if(rx - lx >= 2) ret = min_func(ret, row[1][lx].query(ly, ry)), lx += 2;
      if(rx > lx)      ret = min_func(ret, row[0][lx].query(ly, ry));
      return ret;
    }
    if(lyb > ryb){
      if(ry - ly >= 2) ret = min_func(ret, col[1][ly].query(lx, rx)), ly += 2;
      if(ry > ly)      ret = min_func(ret, col[0][ly].query(lx, rx));
      return ret;
    }
    if(lxceil - lx >= 2) ret = min_func(ret, row[1][lx].query(ly, ry)), lx += 2;
    if(lxceil - lx >= 1) ret = min_func(ret, row[0][lx].query(ly, ry));
    if(rx - rxfloor >= 2) ret = min_func(ret, row[1][rxfloor].query(ly, ry)), rxfloor += 2;
    if(rx - rxfloor >= 1) ret = min_func(ret, row[0][rxfloor].query(ly, ry));
    if(lyceil - ly >= 2) ret = min_func(ret, col[1][ly].query(lx, rx)), ly += 2;
    if(lyceil - ly >= 1) ret = min_func(ret, col[0][ly].query(lx, rx));
    if(ry - ryfloor >= 2) ret = min_func(ret, col[1][ryfloor].query(lx, rx)), ryfloor += 2;
    if(ry - ryfloor >= 1) ret = min_func(ret, col[0][ryfloor].query(lx, rx));
    return ret;
  }
};
#endif