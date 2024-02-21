#ifndef _POINT_SET_LEX_SORT_H_
#define _POINT_SET_LEX_SORT_H_
#include <vector>
#include <algorithm>
// input: {initial array, changes(index, value)}
// output: ans[i] := priority of Vi(after ith change)
template<typename T>
std::vector<int> point_set_lex_sort(const std::vector<T> &v, const std::vector<std::pair<int, T>> &c){
  int n = v.size();
  std::vector<std::vector<std::pair<int, T>>> C0(n); // C0[i] := index iが変わった{時刻, 値}
  for(int i = 0; i < n; i++) C0[i].push_back({0, v[i]});
  for(int i = 0; i < c.size(); i++){
    auto [pos, val] = c[i];
    C0[pos].push_back({i + 1, val});
  }
  std::vector<std::vector<std::pair<int, int>>> C(n); // C[i] := index iが変わった{時刻, 値}
  for(int i = 0; i < n; i++){
    if(C0[i].empty()) continue;
    std::sort(C0[i].begin(), C0[i].end(), [](const std::pair<int, T> &a, const std::pair<int, T> &b){ return a.second < b.second;});
    C[i].resize(C0[i].size());
    T pre = C0[i][0].second;
    for(int j = 0, k = 0; j < C0[i].size(); j++){
      T tmp = C0[i][j].second;
      if(pre != tmp){
        pre = tmp;
        k++;
      }
      C[i][j] = {C0[i][j].first, k}; 
    }
    std::sort(C[i].begin(), C[i].end());
  }
  C0.clear();
  int max_val = n + c.size(), max_time = c.size() + 1;
  std::vector<std::vector<std::pair<int, int>>> S2(max_val), S1(max_val); // S1[i] := 
  auto solve_f = [&](auto &&solve_f, int l, int r) -> std::vector<std::pair<int, int>> {
    if(r - l == 1){
      auto res = C[l];
      C[l].clear();
      return res;
    }
    int mid = (l + r) / 2;
    auto L = solve_f(solve_f, l, mid);
    auto R = solve_f(solve_f, mid, r);
    int N = L.size(), M = R.size();
    L.push_back({max_time, -1});
    R.push_back({max_time, -1});
    std::vector<std::pair<int, int>> res;
    int prelx = -1, prerx = -1, lidx = 0, ridx = 0;
    while(lidx < N || ridx < M){
      auto [lt, lx] = L[lidx];
      auto [rt, rx] = R[ridx];
      if(lt == rt){
        assert(lt == 0);
        S2[rx].push_back({0, lx});
        res.push_back({0, -1});
        prelx = lx, prerx = rx;
        lidx++, ridx++;
      }else if(lt < rt){
        S2[prerx].push_back({res.size(), lx});
        res.push_back({lt, -1});
        prelx = lx;
        lidx++;
      }else{
        S2[rx].push_back({res.size(), prelx});
        res.push_back({rt, -1});
        prerx = rx;
        ridx++;
      }
    }
    L.clear(), R.clear();
    for(int i = 0; i < M; i++){
      for(auto [t, lx] : S2[i]) S1[lx].push_back({t, i});
      S2[i].clear();
    }
    for(int i = 0, k = -1; i < N; i++){
      prerx = -1;
      for(auto [t, rx] : S1[i]){
        if(prerx != rx){
          k++;
          prerx = rx;
        }
        res[t].second = k;
      }
      S1[i].clear();
    }
    return res;
  };
  auto tmp = solve_f(solve_f, 0, n);
  assert(tmp.size() == c.size() + 1);
  std::vector<int> ans(tmp.size());
  for(int i = 0; i < tmp.size(); i++) ans[i] = tmp[i].second;
  return ans;
}
#endif