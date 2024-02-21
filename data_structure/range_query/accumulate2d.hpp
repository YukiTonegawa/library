#ifndef _ACCUMULATE2d_H_
#define _ACCUMULATE2d_H_
#include <vector>
#include <iostream>
#include <functional>
#include <cassert>

template<typename Val>
struct accumulate2d{
  int n, m;
  std::vector<std::vector<Val>> sum;
  accumulate2d(){}
  accumulate2d(const std::vector<std::vector<Val>> &v){
    n = v.size(), m = (!n ? 0 : v[0].size());
    sum.resize(n + 1, std::vector<Val>(m + 1, 0));
    for(int i = 1; i <= n; i++){
      for(int j = 1; j <= m; j++){
        sum[i][j] = sum[i][j - 1] + v[i - 1][j - 1];
      }
    }
    for(int j = 1; j <= m; j++){
      for(int i = 1; i <= n; i++){
        sum[i][j] += sum[i - 1][j];
      }
    }
  }
  // [lx, rx) × [ly, ry)の和, 範囲外は全て単位元とする
  Val query(int lx, int rx, int ly, int ry){
    lx = std::max(lx, 0);
    ly = std::max(ly, 0);
    rx = std::min(rx, n);
    ry = std::min(ry, m);
    if(lx >= rx || ly >= ry) return 0;
    Val upper_left = sum[lx][ly];
    Val lower_left = sum[rx][ly];
    Val upper_right = sum[lx][ry];
    Val lower_right = sum[rx][ry];
    return lower_right - lower_left - upper_right + upper_left;
  }
};
#endif