#ifndef _BINARY_INDEXED_TREE2D_H_
#define _BINARY_INDEXED_TREE2D_H_
#include <vector>
#include <limits>
#include <cassert>

template<typename Val>
struct binary_indexed_tree2d{
  int h, w;
  std::vector<std::vector<Val>> sum;
  binary_indexed_tree2d(){}
  binary_indexed_tree2d(int _h, int _w): h(_h), w(_w), sum(h + 1, std::vector<Val>(w + 1, 0)){}
  binary_indexed_tree2d(const std::vector<std::vector<Val>> &v): h(v.size()){
    if(v.empty()) return;
    w = v[0].size();
    sum.resize(h + 1, std::vector<Val>(w + 1, 0));
    for(int i = 1; i <= h; i++){
      int next_i = i + (i & (-i));
      for(int j = 1; j <= w; j++){
        sum[i][j] += v[i - 1][j - 1];
        int next_j = j + (j & (-j));
        if(next_i <= h) sum[next_i][j] += sum[i][j];
        if(next_j <= w) sum[i][next_j] += sum[i][j];
      }
    }
  }
  // 点(x, y)にzを足す
  void update(int x, int y, Val z){
    for(int i = x + 1; i <= h; i += (i & (-i))){
      for(int j = y + 1; j <= w; j += (j & (-j))){
        sum[i][j] += z;
      }
    }
  }
  // [0, rx) × [0, ry)の和
  Val query(int rx, int ry){
    Val ret = 0;
    for(int i = rx; i > 0; i -= (i & (-i))){
      for(int j = ry; j > 0; j -= (j & (-j))){
        ret += sum[i][j];
      }
    }
    return ret;
  }
  // [lx, rx)　× [ly, ry)の和
  Val query(int lx, int rx, int ly, int ry){
    if(lx >= rx || ly >= ry) return 0;
    return query(rx, ry) - query(rx, ly) - query(lx, ry) + query(lx, ly);
  }
};
#endif