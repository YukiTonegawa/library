#ifndef _RANGE_FREQ_H_
#define _RANGE_FREQ_H_
#include <vector>
#include <cassert>
#include <algorithm>

// 全ての値が[0, v.size()))
struct static_range_freq{
  int n;
  std::vector<std::vector<int>> pos; // pos[x] := xのある位置の集合
  static_range_freq(const std::vector<int> &v): n(v.size()), pos(n){
    for(int i = 0; i < n; i++){
      assert(0 <= v[i] && v[i] < n);
      pos[v[i]].push_back(i);
    }
  }
  // [0, r)のxの個数 O(logN)
  int rank(int r, int x){
    assert(0 <= x && x < n);
    return std::lower_bound(pos[x].begin(), pos[x].end(), r) - pos[x].begin();
  }
  int freq(int l, int r, int x){
    return rank(r, x) - rank(l, x);
  }
  // k番目のxの場所, ない場合は-1 O(1)
  int select(int k, int x){
    assert(0 <= x && x < n);
    return pos[x].size() <= k ? -1 : pos[x][k];
  }
};
#endif