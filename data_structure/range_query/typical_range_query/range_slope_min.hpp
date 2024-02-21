#ifndef _RANGE_SLOPE_MIN_H_
#define _RANGE_SLOPE_MIN_H_
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <deque>

template<typename Val>
struct monotone_cht{
private:
  static constexpr Val inf = std::numeric_limits<Val>::max();
  std::deque<std::tuple<Val, Val, Val>> v; // 最小値を取る可能性がある{前の値との傾きの境界, a, b}
  // ceil(a / b)
  Val ceildiv(Val a, Val b){
    assert(b > 0);
    Val c = a / b;
    if(a > 0 && a - b * c) c++;
    return c;
  }
public:
  monotone_cht(){}
  // ax+bを追加, ならしO(1)
  // 過去に追加した任意の直線より大きい傾きをもつ
  void push_back(Val a, Val b){
    if(v.empty()){
      v.push_back({inf, a, b});
      return;
    }
    assert(std::get<1>(v.back()) < a);
    while(v.size() > 1){
      auto [__, e, f] = v[(int)v.size() - 2];
      auto [_, c, d] = v.back();
      // (f - b) / (a - e) < (f - d) / (c - e)
      if(((__int128_t)f - b) * (c - e) < ((__int128_t)f - d) * (a - e)) break;
      v.pop_back();
    }
    // これ未満の場合 ax+b < cx+d
    auto [_, c, d] = v.back();
    Val t = ceildiv(d - b, a - c);
    v.push_back({t, a, b});
  }
  // ax+bを追加, ならしO(1)
  // 過去に追加した任意の直線より小さい傾きをもつ
  void push_front(Val a, Val b){
    if(v.empty()){
      v.push_back({inf, a, b});
      return;
    }
    assert(a < std::get<1>(v[0]));
    while(v.size() > 1){
      auto [__, e, f] = v[1];
      auto [_, c, d] = v[0];
      // (b - d) / (c - a) > (b - f) / (e - a)
      if(((__int128_t)b - d) * (e - a) > ((__int128_t)b - f) * (c - a)) break;
      v.pop_front();
    }
    // これ未満の場合 ax+b < cx+d
    auto [_, c, d] = v[0];
    Val t = ceildiv(b - d, c - a);
    v[0] = {t, c, d};
    v.push_front({inf, a, b});
  }
  // xでの最小値, O(logN)
  Val query(Val x){
    int i = std::partition_point(v.begin(), v.end(), [x](auto t){return std::get<0>(t) > x;}) - v.begin();
    assert(0 < i);
    auto [_, a, b] = v[i - 1];
    return a * x + b;
  }
};

template<typename Val>
struct offline_range_slope_min{
private:
  std::vector<Val> v;
  struct Query{
    int id, l, r;
    Val a, b;
  };
  std::vector<Query> Q;
  void solve_inner(int l, int r, int ql, int qr, std::vector<Val> &ans){
    if(ql == qr) return;
    if(r - l == 1){
      for(int i = ql; i < qr; i++) ans[Q[i].id] = std::min(ans[Q[i].id], v[l] + l * Q[i].a + Q[i].b);
      return;
    }
    int m = (l + r) / 2;
    // mより左の区間, mより右の区間, mを含む区間
    std::vector<Query> tmpl, tmpr, tmpm;
    for(int i = ql; i < qr; i++){
      if(Q[i].r < m) tmpl.push_back(Q[i]);
      else if(Q[i].l > m) tmpr.push_back(Q[i]);
      else tmpm.push_back(Q[i]);
    }
    int lsz = tmpl.size();
    std::copy(tmpl.begin(), tmpl.end(), Q.begin() + ql);
    std::copy(tmpr.begin(), tmpr.end(), Q.begin() + ql + lsz);
    solve_inner(l, m, ql, ql + lsz, ans);
    solve_inner(m, r, ql + lsz, ql + lsz + tmpr.size(), ans);
    // 右側
    std::sort(tmpm.begin(), tmpm.end(), [](Query &a, Query &b){return a.r < b.r;});
    monotone_cht<Val> hr;
    int j = 0;
    while(j < tmpm.size() && tmpm[j].r == m) j++;
    for(int i = m; i < r; i++){
      hr.push_back(i, v[i]);
      while(j < tmpm.size() && i + 1 == tmpm[j].r){
        Val y = hr.query(tmpm[j].a) + tmpm[j].b;
        ans[tmpm[j].id] = std::min(ans[tmpm[j].id], y);
        j++;
      }
    }
    // 左側
    std::sort(tmpm.begin(), tmpm.end(), [](Query &a, Query &b){return a.l > b.l;});
    monotone_cht<Val> hl;
    j = 0;
    while(j < tmpm.size() && tmpm[j].l == m) j++;
    for(int i = m - 1; i >= l; i--){
      hl.push_front(i, v[i]);
      while(j < tmpm.size() && i == tmpm[j].l){
        Val y = hl.query(tmpm[j].a) + tmpm[j].b;
        ans[tmpm[j].id] = std::min(ans[tmpm[j].id], y);
        j++;
      }
    }
  }
public:
  offline_range_slope_min(const std::vector<Val> &v): v(v){}
  // [l, r)にax + bを足したと仮定したときの[l, r)のmin
  void query(int l, int r, Val a, Val b){
    int id = Q.size();
    Q.push_back({id, l, r, a, b});
  }
  std::vector<Val> solve(){
    Val inf = std::numeric_limits<Val>::max();
    std::vector<Val> ans(Q.size(), inf);
    solve_inner(0, v.size(), 0, Q.size(), ans);
    return ans;
  }
};
#endif