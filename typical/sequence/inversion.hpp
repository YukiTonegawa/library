#ifndef _INVERSION_H_
#define _INVERSION_H_
#include <vector>
#include <algorithm>
#include "../../minior/binary_indexed_tree_set.hpp"
// 転倒数
template<typename T>
long long inversion(const std::vector<T> &v){
  int n = v.size();
  std::vector<std::pair<T, int>> V(n);
  for(int i = 0; i < n; i++) V[i] = {v[i], i};
  std::sort(V.begin(), V.end());
  binary_indexed_tree_set bit(n);
  long long ans = 0;
  for(auto [_, i] : V){
    ans += bit.rank0(i);
    bit.insert(i);
  }
  return ans;
}
#endif