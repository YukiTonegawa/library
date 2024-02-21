#ifndef _RANGE_INVERSION_H_
#define _RANGE_INVERSION_H_
#include "../mo_algorithm.hpp"
#include "../../../minior/binary_indexed_tree_set.hpp"

template<typename T>
std::vector<long long> offline_static_range_inversion(const std::vector<T> &V, const std::vector<std::pair<int, int>> &Q){
  struct range_inversion_st{
    binary_indexed_tree_set b;
    std::vector<int> a;
    int n, all = 0;
    long long sum = 0;
    range_inversion_st(int n): n(n){}
    void add_left(int i){
      all++;
      sum += b.rank1(a[i]);
      b.insert(a[i]);
    }
    void del_left(int i){
      all--;
      sum -= b.rank1(a[i]);
      b.erase(a[i]);
    }
    void add_right(int i){
      sum += all - b.rank1(a[i] + 1);
      all++;
      b.insert(a[i]);
    }
    void del_right(int i){
      sum -= all - b.rank1(a[i] + 1);
      all--;
      b.erase(a[i]);
    }
  };
  int n = V.size(), q = Q.size();
  range_inversion_st ri(n);
  std::vector<int> tmp(n);
  std::vector<std::pair<T, int>> z;
  for(int i = 0; i < n; i++) z.push_back({V[i], i});
  std::sort(z.begin(), z.end());
  for(int i = 0; i < n; i++) tmp[z[i].second] = i;
  ri.a = tmp;
  ri.b = binary_indexed_tree_set(n);
  mo_algorithm mo(n, ri);
  for(auto [l, r] : Q) mo.insert(l, r);
  std::vector<long long> ans(q);
  while(true){
    auto [i, id] = mo.process();
    if(i == -1) break;
    ans[i] = ri.sum;
  }
  return ans;
}
#endif