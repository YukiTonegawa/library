#ifndef _BACKTRACK_H_
#define _BACKTRACK_H_
#include <vector>
#include <cassert>
#include <queue>
#include <array>
struct backtrack{
  int n;
  std::vector<std::vector<int>> g, grev;
  std::vector<int> outdeg;
  backtrack(const std::vector<std::vector<int>> &g): n(g.size()), g(g), grev(n), outdeg(n, 0){
    for(int i = 0; i < n; i++){
      outdeg[i] = g[i].size();
      for(int j : g[i]) grev[j].push_back(i);
    }
  }
  // 0: 先手の勝ち
  // 1: 後手の勝ち
  // -1: 未決定(互いに負けを押し付け合う状態, 奇数手なら先手の勝ち, 偶数手なら後手の勝ち)

  // val[i] := {頂点iの先手から, 頂点iの後手から}
  std::vector<std::array<int, 2>> solve(const std::vector<std::array<int, 2>> &val){
    auto ans = val;
    std::vector<std::array<int, 2>> deg(n, {0, 0});
    // 移動可能な頂点が無く, 初期値が勝ちでない頂点を負けにする
    for(int i = 0; i < n; i++){
      if(outdeg[i]) continue;
      if(ans[i][0] == -1) ans[i][0] = 1;
      if(ans[i][1] == -1) ans[i][1] = 0;
    }
    std::queue<std::pair<int, int>> q; // {頂点, 手番}
    for(int i = 0; i < n; i++){
      if(ans[i][0] != -1) q.push({i, 0});
      if(ans[i][1] != -1) q.push({i, 1});
    }
    while(!q.empty()){
      auto [v, fs] = q.front();
      q.pop();
      assert(ans[v][fs] == 0 || ans[v][fs] == 1);
      if(fs != ans[v][fs]){
        for(int u : grev[v]){
          assert(ans[u][!fs] != fs);
          if(ans[u][!fs] == -1) ans[u][!fs] = !fs, q.push({u, !fs});
        }
      }else{
        for(int u : grev[v]){
          if(ans[u][!fs] == -1 && ++deg[u][!fs] == outdeg[u]){
            ans[u][!fs] = fs;
            q.push({u, !fs});
          }
        }
      }
    }
    return ans;
  }
};
#endif
