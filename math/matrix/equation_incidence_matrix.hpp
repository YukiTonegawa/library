#ifndef _EQUATION_INCIDENCE_MATRIX_H_
#define _EQUATION_INCIDENCE_MATRIX_H_
#include <vector>
#include <cassert>
// 有向グラフにおいて
// 各辺(1 <= i <= |E|)に重みwiを割り当てる
// 全ての辺について S_始点 -= wi, S_終点 += wi
// としたときSが一致するような辺への重みの割り当てを求める
// 存在しない場合空の配列を返す
template<typename T>
std::vector<T> linear_equation_on_incidence_matrix(std::vector<T> S, const std::vector<std::pair<int, int>> &E){
  int n = S.size(), m = E.size();
  std::vector<std::vector<std::pair<int, int>>> E2(n); // {終点, 辺番号}
  for(int i = 0; i < m; i++){
    auto [s, t] = E[i];
    assert(0 <= s && s < n);
    assert(0 <= t && t < n);
    E2[t].push_back({s, i + m}); // m以上の場合逆辺
    E2[s].push_back({t, i});
  }
  std::vector<T> ans(m, 0);
  std::vector<bool> used(n, 0);
  auto dfs = [&](auto &&dfs, int v) -> void {
    assert(!used[v]);
    used[v] = 1;
    for(auto [to, id] : E2[v]){
      if(used[to]) continue;
      dfs(dfs, to);
      if(id < m) ans[id] = S[to];
      else ans[id - m] = -S[to];
      S[v] -= S[to];
    }
  };
  for(int i = 0; i < n; i++){
    if(!used[i]){
      dfs(dfs, i);
      if(S[i] != 0) return {}; // 解なし
    }
  }
  return ans;
}
#endif