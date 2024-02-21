#ifndef _SORT_H_
#define _SORT_H_

#include <vector>
#include <cassert>
// {begin, end, f: Itr::value_type -> uint}
template<typename Itr, typename F>
void bucket_sort(Itr first, Itr last, F f){
  std::vector<std::vector<typename Itr::value_type>> tmp(1);
  int sz = 1;
  for(auto itr = first; itr != last; itr++){
    int key = f(*itr);
    while(sz <= key) tmp.resize(sz *= 2);
    tmp[key].push_back(*itr);
  }
  int row = 0, col = 0;
  for(auto itr = first; itr != last; itr++, col++){
    while(tmp[row].size() <= col) row++, col = 0;
    *itr = tmp[row][col];
  }
}
// keyが[0, max_elem)
template<typename Itr, typename F>
void bucket_sort_fast(Itr first, Itr last, F f, int max_elem){
  static std::vector<std::vector<typename Itr::value_type>> tmp;
  if(tmp.size() < max_elem) tmp.resize(max_elem);
  std::vector<int> idx(max_elem, 0);
  for(auto itr = first; itr != last; itr++){
    int key = f(*itr);
    if(idx[key] == tmp[key].size()){
      tmp[key].push_back(*itr);
      idx[key]++;
    }else tmp[key][idx[key]++] = *itr;
  }
  auto itr = first;
  for(int row = 0; row < max_elem && itr != last; row++){
    for(int col = 0; col < idx[row]; col++, itr++){
      *itr = tmp[row][col];
    }
  }
}
// ソート済の[l, m)と[m, r)をマージ(inplace)
template<typename T, typename Compare>
void __merge_sorted(std::vector<T> &a, int l, int m, int r, Compare f){
  static std::vector<T> tmp;
  if(tmp.size() < r - l) tmp.resize(r - l);
  int i = l, j = m, k = 0;
  while(i < m || j < r){
    if(i == m){
      tmp[k++] = a[j++];
    }else if(j == r || f(a[i], a[j])){
      tmp[k++] = a[i++];
    }else{
      tmp[k++] = a[j++];
    }
  }
  for(i = 0; i < r - l; i++) a[l + i] = tmp[i];
}
template<typename T, typename Compare>
void merge_sort(std::vector<T> &a, Compare f){
  auto calc = [&](auto &&calc, int l, int r) -> void {
    if(r - l <= 1) return;
    int m = (l + r) / 2;
    calc(calc, l, m);
    calc(calc, m, r);
    __merge_sorted(a, l, m, r, f);
  };
  calc(calc, 0, a.size());
}
// ソート済の[l1, r1), [l2, r2)を比較関数fでマージ
template<typename Itr1, typename Itr2, typename Compare>
std::vector<typename Itr1::value_type> __merge_sorted(Itr1 l1, Itr1 r1, Itr2 l2, Itr2 r2, Compare f){
  std::vector<typename Itr1::value_type> res;
  std::merge(l1, r1, l2, r2, std::back_inserter(res), f);
  return res;
}
template<typename T>
std::vector<T> __merge_unique(const std::vector<T> &L, const std::vector<T> &R){
  int n = L.size(), m = R.size();
  int lidx = 0, ridx = 0;
  std::vector<T> res;
  while(lidx < n || ridx < m){
    if(lidx == n){
      if(res.empty() || R[ridx] != res.back()) res.push_back(R[ridx]);
      ridx++;
    }else if(ridx == m || L[lidx] < R[ridx]){
      if(res.empty() || L[lidx] != res.back()) res.push_back(L[lidx]);
      lidx++;
    }else{
      if(res.empty() || R[ridx] != res.back()) res.push_back(R[ridx]);
      ridx++;
    }
  }
  return res;
}
#endif