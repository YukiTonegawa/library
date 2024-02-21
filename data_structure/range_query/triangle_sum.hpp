#ifndef _TRIANGLE_SUM_H_
#define _TRIANGLE_SUM_H_
#include <vector>
#include <cassert>

template<typename Val>
struct triangle_sum{ 
private:
  using vec = std::vector<Val>;
  using mat = std::vector<vec>;
  int n, m;
  mat rect, lltri, ultri;
  // 範囲外でないことが保証されている
  Val __lower_left_tri(int lx, int rx, int ly){
    int ry = ly + (rx - lx);
    assert(0 <= lx && lx <= rx && rx <= n);
    assert(0 <= ly && ly <= ry && ry <= m);
    return lltri[rx][ry] - lltri[lx][ly] - query_rect(lx, rx, 0, ly);
  }
  Val __upper_left_tri(int lx, int rx, int ly){
    int ry = ly + (rx - lx);
    assert(0 <= lx && lx <= rx && rx <= n);
    assert(0 <= ly && ly <= ry && ry <= m);
    return ultri[rx][ly + 1] - (ry == m ? 0 : ultri[lx][ry + 1]) - query_rect(lx, rx, 0, ly);
  }
public:
  triangle_sum(const mat &v): n(v.size()), m(n ? v[0].size() : 0), rect(n + 1, vec(m + 1)), lltri(n + 1, vec(m + 1)), ultri(n + 1, vec(m + 1)){
    std::fill(rect[0].begin(), rect[0].end(), 0);
    std::fill(lltri[0].begin(), lltri[0].end(), 0);
    std::fill(ultri[0].begin(), ultri[0].end(), 0);
    for(int i = 1; i <= n; i++){
      rect[i][0] = 0;
      assert(v[i - 1].size() == m);
      for(int j = 1; j <= m; j++) rect[i][j] = rect[i][j - 1] + v[i - 1][j - 1];
    }
    for(int i = 1; i <= n; i++){
      lltri[i][0] = 0;
      ultri[i][m] = rect[i][m];
      for(int j = 0; j < m; j++) ultri[i][j] += ultri[i - 1][j + 1] + rect[i][j];
      for(int j = 1; j <= m; j++){
        lltri[i][j] = lltri[i - 1][j - 1] + rect[i][j];
        rect[i][j] += rect[i - 1][j];
      }
    }
  }
  // sum [lx, rx) × [ly, rx)
  // 範囲外には0が書かれているとする
  Val query_rect(int lx, int rx, int ly, int ry){
    lx = std::max(lx, 0);
    rx = std::min(rx, n);
    ly = std::max(ly, 0);
    ry = std::min(ry, m);
    if(lx >= rx || ly >= ry) return 0;
    return rect[rx][ry] - rect[lx][ry] - rect[rx][ly] + rect[lx][ly];
  }
  /*
  範囲外には0が書かれているとする
  x x x
  . x x
  . . x
  */
  Val query_upper_right_tri(int lx, int rx, int ly){
    int ry = ly + (rx - lx);
    return query_rect(lx, rx, ly, ry) - query_lower_left_tri(lx + 1, rx, ly);
  }
  /*
  範囲外には0が書かれているとする
  x . .  
  x x .  
  x x x 
  */
  Val query_lower_left_tri(int lx, int rx, int ly){
    if(lx >= n || ly >= m) return 0;
    rx = std::min(rx, n);
    Val res = 0;
    if(lx < 0){
      res += query_rect(0, rx, ly, ly - lx);
      ly -= lx;
      lx = 0;
    }
    if(ly < 0){
      lx -= ly;
      ly = 0;
    }
    if(lx >= rx || ly >= m) return res;
    // この時点で0 <= lx < rx <= n かつ 0 <= ly < m
    int ry = ly + (rx - lx);
    if(ry > m){
      int diff = ry - m;
      res += query_rect(rx - diff, rx, ly, m);
      rx -= diff;
    }
    return res + __lower_left_tri(lx, rx, ly);
  }
  /*
  範囲外には0が書かれているとする
  x x x
  x x .  
  x . . 
  */
  Val query_upper_left_tri(int lx, int rx, int ly){
    if(lx >= n || ly >= m) return 0;
    Val res = 0;
    lx = std::max(lx, 0);
    if(rx > n){
      int diff = rx - n;
      res += query_rect(lx, n, ly, ly + diff);
      rx = n;
      ly += diff;
    }
    if(ly < 0){
      rx += ly;
      ly = 0;
    }
    if(lx >= rx || ly >= m) return res;
    // この時点で0 <= lx < rx <= n かつ 0 <= ly < m
    int ry = ly + (rx - lx);
    if(ry > m){
      int diff = ry - m;
      res += query_rect(lx, lx + diff, ly, m);
      lx += diff;
    }
    return res + (lx < rx ? __upper_left_tri(lx, rx, ly) : 0);
  }
  /*
  範囲外には0が書かれているとする
  . . x
  . x x
  x x x
  */
  Val query_lower_right_tri(int lx, int rx, int ly){
    int ry = ly + (rx - lx);
    return query_rect(lx, rx, ly, ry) - query_upper_left_tri(lx, rx - 1, ly);
  }
};
#endif