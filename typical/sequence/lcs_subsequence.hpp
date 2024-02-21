#ifndef _LCS_H_
#define _LCS_H_
#include <vector>
#include <algorithm>
#include <cassert>
template<typename T>
int lcs(const std::vector<T> &a, const std::vector<T> &b){
  int n = a.size(), m = b.size();
  std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
  for(int i = 0; i <= n; i++){
    for(int j = 0; j <= m; j++){
      if(i) dp[i][j] = dp[i - 1][j];
      if(j) dp[i][j] = std::max(dp[i][j], dp[i][j - 1]);
      if(i && j && a[i - 1] == b[j - 1]) dp[i][j] = std::max(dp[i][j], dp[i - 1][j - 1] + 1);
    }
  }
  return dp[n][m];
}

template<typename T>
std::vector<T> lcs_restore(const std::vector<T> &a, const std::vector<T> &b){
  int n = a.size(), m = b.size();
  std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
  std::vector<std::vector<std::pair<int, int>>> prev(n + 1, std::vector<std::pair<int, int>>(m + 1, {-1, -1}));
  for(int i = 0; i <= n; i++){
    for(int j = 0; j <= m; j++){
      if(i && dp[i][j] < dp[i - 1][j]){
        dp[i][j] = dp[i - 1][j];
        prev[i][j] = {i - 1, j};
      }
      if(j && dp[i][j] < dp[i][j - 1]){
        dp[i][j] = dp[i][j - 1];
        prev[i][j] = {i, j - 1};
      }
      if(i && j && a[i - 1] == b[j - 1] && dp[i][j] <= dp[i - 1][j - 1]){
        dp[i][j] = dp[i - 1][j - 1] + 1;
        prev[i][j] = {i - 1, j - 1};
      }
    }
  }
  std::vector<T> res;
  int x = dp[n][m], y = n, z = m;
  while(x){
    std::tie(y, z) = prev[y][z];
    if(x != dp[y][z]){
      assert(a[y] == b[z]);
      res.push_back(a[y]);
    }
    x = dp[y][z];
  }
  std::reverse(res.begin(), res.end());
  return res;
}
#endif