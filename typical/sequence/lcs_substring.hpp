#ifndef _LCSUBSTRING_H_
#define _LCSUBSTRING_H_
#include "../../string/string.hpp"
// {長さ, sの始点, tの始点}
template<typename T>
std::tuple<int, int, int> longest_common_substring(std::vector<T> s, const std::vector<T> &t){
  static constexpr T minf = std::numeric_limits<T>::min();
  if(s.empty() || t.empty()) return {0, 0, 0};
  int n = s.size(), m = t.size();
  s.push_back(minf);
  s.insert(s.begin() + n + 1, t.begin(), t.end());
  auto sa = suffix_array<T>(s);
  auto lcp = lcp_array<T>(s, sa);
  int lenmx = 0, lenmxi = -1;
  for(int i = 1; i + 1 < sa.size(); i++){
    if((sa[i] > n && sa[i + 1] < n) || (sa[i] < n && sa[i + 1] > n)){
      int len = min({n, m, lcp[i]});
      if(lenmx < len) lenmx = len, lenmxi = i;
    }
  }
  if(lenmxi == -1) return {0, 0, 0};
  int l1 = sa[lenmxi], l2 = sa[lenmxi + 1];
  if(l1 > n) std::swap(l1, l2);
  l2 -= n + 1;
  return {lenmx, l1, l2};
}
#endif