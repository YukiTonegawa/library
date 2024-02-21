#ifndef _CONSTRUCTIVE_GRAPH_H_
#define _CONSTRUCTIVE_GRAPH_H_
#include <vector>
#include <cassert>
#include <numeric>
#include <algorithm>
#include <unordered_map>

// |E|がO(N^2)の最大マッチングの特殊な場合

// x, y(共に長さn)で同じ値に辺を張ってはいけない2部グラフの完全マッチング
template<typename T>
std::vector<int> contrast(const std::vector<T> &x, const std::vector<T> &y, bool is_sorted = false){
  assert(x.size() == y.size());
  int n = x.size();
  std::vector<int> x_idx(n), y_idx, ret(n);
  std::iota(x_idx.begin(), x_idx.end(), 0);
  y_idx = x_idx;
  if(is_sorted){
    std::reverse(y_idx.begin(), y_idx.end());
  }else{
    std::sort(x_idx.begin(), x_idx.end(), [&x](int A, int B){return x[A] < x[B];});
    std::sort(y_idx.begin(), y_idx.end(), [&y](int A, int B){return y[A] > y[B];});
  }
  T dup_element;
  for(int i=0;i<n;i++){
    if(x[x_idx[i]] == y[y_idx[i]]){
      dup_element = x[x_idx[i]];
      break;
    }
  }
  std::vector<int> safe, dup;
  for(int i=0;i<n;i++){
    if(x[x_idx[i]] == y[y_idx[i]]){
      dup.push_back(i);
    }
    if(x[x_idx[i]] != dup_element && y[y_idx[i]] != dup_element){
      safe.push_back(i);
    }
  }
  if(!dup.empty()){
    if(dup.size() > safe.size()) return std::vector<int>();
    while(!dup.empty()){
      std::swap(y_idx[dup.back()], y_idx[safe.back()]);
      dup.pop_back();
      safe.pop_back();
    }
  }
  for(int i=0;i<n;i++) ret[x_idx[i]] = y_idx[i];
  return ret;
}

// x(長さn), y(長さm)で同じ値に辺を張ってはいけない2部グラフの最大マッチング
template<typename T>
std::pair<int, std::vector<int>> contrast_maximum_matching(const std::vector<T> &x, const std::vector<T> &y, bool is_sorted = false){
  int n = x.size(), m = y.size();
  std::vector<int> ret(n, -1), x_idx, y_idx;;
  std::unordered_map<T, int> cnt;
  if(n < m){
    x_idx.resize(n);
    std::iota(x_idx.begin(), x_idx.end(), 0);
    for(int i=0;i<n;i++) cnt[x[i]]++;
    for(int i=0;i<m;i++){
      if(cnt[y[i]] < n || n - y_idx.size() == m - i){
        y_idx.push_back(i);
        cnt[y[i]]++;
      }
      if(y_idx.size() == n) break;
    }
  }else if(n > m){
    y_idx.resize(m);
    std::iota(y_idx.begin(), y_idx.end(), 0);
    for(int i=0;i<m;i++) cnt[y[i]]++;
    for(int i=0;i<n;i++){
      if(cnt[x[i]] < m || m - x_idx.size() == n - i){
        x_idx.push_back(i);
        cnt[x[i]]++;
      }
      if(x_idx.size() == m) break;
    }
  }
  n = std::min(n, m);
  if(is_sorted){
    std::reverse(y_idx.begin(), y_idx.end());
  }else{
    std::sort(x_idx.begin(), x_idx.end(), [&x](int A, int B){return x[A] < x[B];});
    std::sort(y_idx.begin(), y_idx.end(), [&y](int A, int B){return y[A] > y[B];});
  }
  T dup_element;
  for(int i=0;i<n;i++){
    if(x[x_idx[i]] == y[y_idx[i]]){
      dup_element = x[x_idx[i]];
      break;
    }
  }
  std::vector<int> safe, dup;
  for(int i=0;i<n;i++){
    if(x[x_idx[i]] == y[y_idx[i]]){
      dup.push_back(i);
    }
    if(x[x_idx[i]] != dup_element && y[y_idx[i]] != dup_element){
      safe.push_back(i);
    }
  }
  if(!dup.empty()){
    while(!dup.empty() && !safe.empty()){
      std::swap(y_idx[dup.back()], y_idx[safe.back()]);
      dup.pop_back();
      safe.pop_back();
    }
  }
  for(int i=0;i<n;i++) ret[x_idx[i]] = y_idx[i];
  for(int d:dup) ret[x_idx[d]] = -1;
  return {n - dup.size(), ret};
}
#endif
