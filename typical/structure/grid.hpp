#ifndef _GRID_H_
#define _GRID_H_
#include <vector>
#include <cassert>

struct grid{
  /*
  0 1 2 3 4
  5 6 7 8 9
  ...
  */
  static std::vector<std::pair<int, int>> order_yoko(int n, int m){
    std::vector<std::pair<int, int>> res;
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        res.push_back({i, j});
      }
    }
    return res;
  }
  /*
  0 4 ...
  1 5 ...
  2 6 ...
  3 7 ...
  */
  static std::vector<std::pair<int, int>> order_tate(int n, int m){
    std::vector<std::pair<int, int>> res;
    for(int j = 0; j < m; j++){
      for(int i = 0; i < n; i++){
        res.push_back({i, j});
      }
    }
    return res;
  }
  /*
  0 1 3 6
  2 4 7..
  5 8....
  9......
  */
  static std::vector<std::pair<int, int>> order_naname(int n, int m){
    std::vector<std::pair<int, int>> res;
    for(int j =  0; j < m; j++){
      for(int i = 0; i < std::min(n, j + 1); i++){
        int x = i, y = j - i;
        res.push_back({x, y});
      }
    }
    for(int j = n - 2; j >= 0; j--){
      for(int i = 0; i < std::min(m, j + 1); i++){
        int x = (n - 1 - j) + i, y = m - 1 - i;
        res.push_back({x, y});
      }
    }
    return res;
  }
  /*
  0  1  2  3  4
  15 16 17 18 5
  14 23 24 19 6
  13 22 21 20 7
  12 11 10 9  8
  */
  static std::vector<std::pair<int, int>> order_circle(int n, int m){
    assert(n && m);
    std::vector<bool> used(n * m, false);
    used[0] = true;
    std::vector<std::pair<int, int>> res{{0, 0}};
    int dx[4] = {0, 1, 0, -1}, dy[4] = {1, 0, -1, 0};
    int x = 0, y = 0, dir = 0;
    while((int)res.size() < n * m){
      int a = x + dx[dir], b = y + dy[dir];
      if(a < 0 || a >= n || b < 0 || b >= m || used[a * m + b]){
        dir = (dir + 1) & 3;
        continue;
      }
      used[a * m + b] = true;
      res.push_back({a, b});
      x = a, y = b;
    }
    return res;
  }
  // 右に90度回転
  template<typename T>
  static std::vector<std::vector<T>> rotate(const std::vector<std::vector<T>> &v){
    int n = v.size();
    if(!n) return {{}};
    int m = v[0].size();
    std::vector<std::vector<T>> res(m, std::vector<T>(n));
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        res[j][n - 1 - i] = v[i][j];
      }
    }
    return res;
  }
  // 転置
  template<typename T>
  static std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>> &v){
    int n = v.size();
    if(!n) return {{}};
    int m = v[0].size();
    std::vector<std::vector<T>> res(m, std::vector<T>(n));
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        res[j][i] = v[i][j];
      }
    }
    return res;
  }
};
#endif