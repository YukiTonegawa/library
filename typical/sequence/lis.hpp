#ifndef _LIS_H_
#define _LIS_H_

#include <vector>
#include <algorithm>
// LIS(元の数列のインデックス)
template<typename T>
std::vector<int> lis(const std::vector<T> &a){
  int n = a.size();
  if(!n) return {};
  std::vector<std::pair<T, int>> v{{a[0], 0}};
  std::vector<int> rev(n, -1);
  for(int i = 1; i < n; i++){
    int j = std::lower_bound(v.begin(), v.end(), std::pair<T, int>{a[i], 0}) - v.begin();
    if(j) rev[i] = v[j - 1].second;
    if(j == v.size()) v.push_back({a[i], i});
    else v[j] = {a[i], i};
  }
  std::vector<int> ans(v.size());
  int s = v.back().second, i = v.size() - 1;
  while(s != -1){
    ans[i--] = s;
    s = rev[s];
  }
  return ans;
}
#endif