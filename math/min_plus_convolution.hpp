#ifndef _MIN_PLUS_CONVOLUTION_H_
#define _MIN_PLUS_CONVOLUTION_H_
#include <vector>
#include <limits>
#include <queue>

//aが下に凸な関数として、 c_k := min(a_i + b_j, i + j = k)
template<typename T>
std::vector<T> min_plus_convolution(const std::vector<T> &a, const std::vector<T> &b){
  static constexpr T inf = std::numeric_limits<T>::max();
  int n = a.size(), m = b.size();
  std::vector<T> c(n + m - 1);
  //kを固定した時のiのargmin, 探索範囲は[l, r)
  auto calculate_argmin = [&](int k, int l, int r)->std::pair<int, T>{
    T ret = inf;
    int idx = l;
    for(int i = l; i < r; i++){
      T val = (i >= m || k < i || k - i >= n ? inf : a[k - i] + b[i]);
      if(val < ret) ret = val, idx = i;
    }
    return {idx, ret};
  };
  std::queue<std::tuple<int, int, int, int>> q; //l_k, r_k, l_i, r_i
  q.push({0, n + m - 1, 0, std::max(n, m)});
  while(!q.empty()){
    auto [l_k, r_k, l_i, r_i] = q.front();
    q.pop();
    int mid_k = (r_k + l_k) / 2;
    auto calc = calculate_argmin(mid_k, l_i, r_i);
    c[mid_k] = calc.second;
    if(r_k > mid_k + 1) q.push({mid_k + 1, r_k, calc.first, r_i});
    if(mid_k > l_k) q.push({l_k, mid_k, l_i, calc.first + 1});
  }
  return c;
}

//aが上に凸な関数として、 c_k := max(a_i + b_j, i + j = k)
template<typename T>
std::vector<T> max_plus_convolution(const std::vector<T> &a, const std::vector<T> &b){
  static constexpr T minf = std::numeric_limits<T>::min();
  int n = a.size(), m = b.size();
  std::vector<T> c(n + m - 1);
  //kを固定した時のiのargmax, 探索範囲は[l, r)
  auto calculate_argmax = [&](int k, int l, int r)->std::pair<int, T>{
    T ret = minf;
    int idx = l;
    for(int i = l; i < r; i++){
      T val = (i >= m || k < i || k - i >= n ? minf : a[k - i] + b[i]);
      if(val > ret) ret = val, idx = i;
    }
    return {idx, ret};
  };
  std::queue<std::tuple<int, int, int, int>> q; //l_k, r_k, l_i, r_i
  q.push({0, n + m - 1, 0, std::max(n, m)});
  while(!q.empty()){
    auto [l_k, r_k, l_i, r_i] = q.front();
    q.pop();
    int mid_k = (r_k + l_k) / 2;
    auto calc = calculate_argmax(mid_k, l_i, r_i);
    c[mid_k] = calc.second;
    if(r_k > mid_k + 1) q.push({mid_k + 1, r_k, calc.first, r_i});
    if(mid_k > l_k) q.push({l_k, mid_k, l_i, calc.first + 1});
  }
  return c;
}
#endif