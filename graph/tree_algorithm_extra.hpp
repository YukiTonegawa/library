#ifndef _TREE_ALGORITHM_EXTRA_H_
#define _TREE_ALGORITHM_EXTRA_H_
#include "tree_algorithm.hpp"
#include "../math/convolution.hpp"

namespace tree_algorithm{
  template<typename edge>
  std::vector<long long> frequency_table_of_tree_distance(std::vector<std::vector<edge>> &G){
    int n = G.size();
    auto [Gc, root, size_i, dep_i, par_i] = centroid_decomposition<edge>(G);
    auto make_deplist = [&](auto &&make_deplist, int v, int p, int d, std::vector<long long> &deplist, const std::vector<int> &i, const int root_depth)->void {
      if(deplist.size() == d) deplist.resize(2 * d);
      deplist[d]++;
      for(auto e: G[v]){
        if(e.t == p || i[e.t] <= root_depth) continue;
        make_deplist(make_deplist, e.t, v, d + 1, deplist, i, root_depth);
      }
    };
    auto convo = [&](std::vector<long long> &a, std::vector<long long> &b)->void {
      static constexpr int lim = 998244353;
      long long asum = std::accumulate(a.begin(), a.end(), 0);
      long long bsum = (&a == &b ? asum : std::accumulate(b.begin(), b.end(), 0));
      if(asum * bsum >= lim) a = convolution_ll(a, b);
      else a = convolution_mod<lim>(a, b);
    };
    auto square = [&](std::vector<long long> &a)->void {
      static constexpr int lim = 998244353;
      long long asum = std::accumulate(a.begin(), a.end(), 0);
      if(asum * asum >= lim) a = square_ll(a);
      else a = square_mod<lim>(a);
    };
    std::vector<long long> ans(n, 0);
    for(int i = 0; i < n; i++){
      if(size_i[i] == 1) continue;
      std::vector<std::vector<long long>> A{{1}};
      std::vector<long long> B{1};
      for(auto e: G[i]){
        if(dep_i[e.t] <= dep_i[i]) continue;
        A.push_back({0});
        make_deplist(make_deplist, e.t, -1, 1, A.back(), dep_i, dep_i[i]);
        for(int k = 0; k < A.back().size(); k++){
          if(B.size() == k) B.resize(2 * k);
          B[k] += A.back()[k];
        }
      }
      // 次数が2の場合, そのまま xy を計算した方が
      // (x + y) ^ 2 - x^2 - y^2 を計算するよりも早い
      if(A.size() == 3){
        std::vector<long long> &x = A[1], &y = A[2];
        for(int j = 0; j < x.size(); j++) ans[j] += x[j] << 1;
        for(int j = 0; j < y.size(); j++) ans[j] += y[j] << 1;
        convo(x, y);
        for(int j = 0; j < std::min(n, (int)x.size()); j++) ans[j] += x[j] << 1;
      }else{
        square(B);
        for(std::vector<long long> &a: A){
          square(a);
          for(int k = 0; k < std::min(n, (int)a.size()); k++) B[k] -= a[k];
        }
        for(int j = 0; j < std::min(n, (int)B.size()); j++) ans[j] += B[j];
      }
    }
    ans[0] = n;
    for(int i = 1; i < n; i++) ans[i] = ans[i] / 2;
    return ans;
  }
  /*
template<typename edge, typename mint>
std::vector<mint> dist_sum_mod(std::vector<std::vector<edge>> &G, vector<mint> val){
  int n = G.size();
  auto [Gc, root, size_i, dep_i, par_i] = tree_algorithm::centroid_decomposition<edge>(G);
  auto make_deplist = [&](auto &&make_deplist, int v, int p, int d, std::vector<mint> &deplist, const std::vector<int> &i, const int root_depth)->void{
    if(deplist.size() == d) deplist.resize(2 * d);
    deplist[d] += val[v];
    for(auto e: G[v]){
      if(e.t == p || i[e.t] <= root_depth) continue;
      make_deplist(make_deplist, e.t, v, d + 1, deplist, i, root_depth);
    }
  };
  std::vector<mint> ans(n, 0);
  for(int i = 0; i < n; i++){
    if(size_i[i] == 1) continue;
    std::vector<std::vector<mint>> A{{val[i]}};
    std::vector<mint> B{val[i]};
    for(auto e: G[i]){
      if(dep_i[e.t] <= dep_i[i]) continue;
      A.push_back({0});
      make_deplist(make_deplist, e.t, -1, 1, A.back(), dep_i, dep_i[i]);
      for(int k = 0; k < A.back().size(); k++){
        if(B.size() == k) B.resize(2 * k);
        B[k] += A.back()[k];
      }
    }
    if(false){
      std::vector<mint> &x = A[1], &y = A[2];
      for(int j = 0; j < x.size(); j++) ans[j] += x[j] * 2;
      for(int j = 0; j < y.size(); j++) ans[j] += y[j] * 2;
      x = convolution_mod<mint>(x, y);
      for(int j = 0; j < std::min(n, (int)x.size()); j++) ans[j] += x[j] * 2;
    }else{
      B = square_mod<mint>(B);
      for(std::vector<mint> &a: A){
        a = square_mod<mint>(a);
        for(int k = 0; k < std::min(n, (int)a.size()); k++) B[k] -= a[k];
      }
      for(int j = 0; j < std::min(n, (int)B.size()); j++) ans[j] += B[j];
    }
  }
  for(int i = 0; i < n; i++) ans[0] += val[i];
  mint i2 = mint(2).inv();
  for(int i = 1; i < n; i++) ans[i] = ans[i] * i2;
  return ans;
}
*/
};

// centroid decomposition bfs


// square root decomposition

#endif