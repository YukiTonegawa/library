#ifndef _RECTANGLE_INTERSECTION_
#define _RECTANGLE_INTERSECTION_
#include <array>
#include <tuple>
#include <algorithm>
#include "rectangle_sum.hpp"
#include "../segment_tree/lazy_segment_tree.hpp"

//       ┃         ┃
//   1   ┃    8    ┃    7
//       ┃         ┃  
// ━━━━━━━━━━━━━━━━━━━━━━━━━━
//       ┃         ┃
//   2   ┃         ┃    6
//       ┃         ┃
// ━━━━━━━━━━━━━━━━━━━━━━━━━━
//       ┃         ┃
//   3   ┃    4    ┃    5
//       ┃         ┃

// 質問される長方形に対して交差する(共通部分の面積が0より大きい)長方形の重みの和は
//                                  (1)                     (2)                     (3)                      (4)
// (追加されている長方形の和) - (右上が1, 2にある長方形) - (左上が3, 4にある長方形) - (左下が5, 6にある長方形) - (右下が7, 8にある長方形)

// これらは点追加長方形和で求められる クエリあたり O(logN)(静的) O(log^2 N)(動的)

template<typename Idx, typename Val>
struct offline_static_rectangle_intersect_sum{
private:
  offline_static_rectangle_sum<Idx, Val> rect_ur, rect_ul, rect_lr, rect_ll;
  static constexpr Idx inf = std::numeric_limits<Idx>::max() / 2;
  static constexpr Idx minf = std::numeric_limits<Idx>::min() / 2;
  bool built = false;
  int q = 0;
  Val sum = 0;
public:
  // [lx, rx) × [ly, ry), 重みzを追加
  void update(Idx lx, Idx rx, Idx ly, Idx ry, Val z){
    assert(!built);
    if(lx == rx || ly == ry) return;
    sum += z;
    // 右上
    rect_ur.update(lx, ry, z);
    // 左上
    rect_ul.update(lx, ly, z);
    // 左下
    rect_ll.update(rx, ly, z);
    // 右下
    rect_lr.update(rx, ry, z);
  }
  // [lx, rx) × [ly, ry)との共通部分の面積が0より大きい長方形の和
  void query(Idx lx, Idx rx, Idx ly, Idx ry){
    assert(!built);
    q++;
    // (1)
    // 1, 2
    rect_ur.query(minf, rx, minf, ly + 1);

    // (2)
    // 3, 4
    rect_ul.query(rx, inf, minf, ry);

    // (3)
    // 5, 6
    rect_ll.query(lx + 1, inf, ry, inf);

    // (4)
    // 7 8
    rect_lr.query(minf, lx + 1, ly + 1, inf);
  }
  std::vector<Val> solve(){
    assert(!built);
    built = true;
    std::vector<Val> ans(q, sum);
    auto qur = rect_ur.solve();
    auto qul = rect_ul.solve();
    auto qll = rect_ll.solve();
    auto qlr = rect_lr.solve();
    for(int i = 0; i < q; i++){
      ans[i] -= qur[i];
      ans[i] -= qul[i];
      ans[i] -= qll[i];
      ans[i] -= qlr[i];
    }
    return ans;
  }
};
template<typename Idx, typename Val>
struct offline_dynamic_rectangle_intersect_sum{
private:
  offline_point_add_rectangle_sum<Idx, Val> rect_ur, rect_ul, rect_lr, rect_ll;
  static constexpr Idx inf = std::numeric_limits<Idx>::max() / 2;
  static constexpr Idx minf = std::numeric_limits<Idx>::min() / 2;
  int q = 0;
  Val sum = 0;
  bool built = false;
  std::vector<Val> ans;
public:
  // [lx, rx) × [ly, ry)を追加
  void add_rect(Idx lx, Idx rx, Idx ly, Idx ry, Val z){
    assert(!built);
    if(lx == rx || ly == ry) return;
    sum += z;
    // 右上
    rect_ur.update(lx, ry, z);
    // 左上
    rect_ul.update(lx, ly, z);
    // 左下
    rect_ll.update(rx, ly, z);
    // 右下
    rect_lr.update(rx, ry, z);
  }
  // [lx, rx) × [ly, ry)との共通部分の面積が0より大きい長方形の数
  void query(Idx lx, Idx rx, Idx ly, Idx ry){
    assert(!built);
    ans.push_back(sum); // クエリ時点での長方形の和
    q++;
    // (1)
    // 1, 2
    rect_ur.query(minf, rx, minf, ly + 1);

    // (2)
    // 3, 4
    rect_ul.query(rx, inf, minf, ry);

    // (3)
    // 5, 6
    rect_ll.query(lx + 1, inf, ry, inf);

    // (4)
    // 7 8
    rect_lr.query(minf, lx + 1, ly + 1, inf);
  }
  std::vector<Val> solve(){
    assert(!built);
    built = true;
    auto qur = rect_ur.solve();
    auto qul = rect_ul.solve();
    auto qll = rect_ll.solve();
    auto qlr = rect_lr.solve();
    for(int i = 0; i < q; i++){
      ans[i] -= qur[i];
      ans[i] -= qul[i];
      ans[i] -= qll[i];
      ans[i] -= qlr[i];
    }
    return ans;
  }
};


// 長方形の和集合
// v[i] = [lx, rx) × [ly, ry)
// @param Idx: インデックス, Idx2: 答えの最大値(インデックスの2乗)
template<typename Idx = int, typename Idx2 = long long>
long long union_of_rectangles(std::vector<std::array<Idx, 4>> v){
  int n = v.size();
  if(!n) return 0;
  struct event{
    Idx x, ly, ry;
    bool is_add;
  };
  std::vector<event> e;
  std::vector<std::pair<Idx, int>> y;
  for(int i = 0; i < n; i++){
    auto [l, r, d, u] = v[i];
    y.push_back({d, i});
    y.push_back({u, i + n});
  }
  std::sort(y.begin(), y.end());
  // y_dist[i] := y[i + 1] - y[i]
  std::vector<std::pair<int, Idx>> y_dist;
  Idx y_unique = 0, prev_y = y[0].first;
  for(int i = 0; i < 2 * n; i++){
    if(y[i].first != prev_y){
      y_dist.push_back({0, y[i].first - prev_y});
      prev_y = y[i].first;
      y_unique++;
    }
    if(y[i].second < n) v[y[i].second][2] = y_unique;
    else v[y[i].second - n][3] = y_unique;
  }
  for(auto [l, r, d, u] : v){
    e.push_back(event{l, d, u, true});
    e.push_back(event{r, d, u, false});
  }
  std::sort(e.begin(), e.end(), [](event &a, event &b){return a.x < b.x;});
  Idx y_all = y.back().first - y[0].first;
  lazy_segment_tree<range_add_range_min_count, std::pair<int, Idx>, int> seg(y_dist);
  Idx2 ans = 0;
  Idx cur_x = 0;
  for(int i = 0; i < 2 * n; i++){
    auto [x, ly, ry, is_add] = e[i];
    auto [m, mcnt] = seg.query_all();
    Idx dx = x - cur_x;
    if(m == 0) ans += (Idx2)(y_all - mcnt) * dx;
    else ans += (Idx2)y_all * dx;
    cur_x = x;
    if(is_add) seg.update(ly, ry, 1);
    else seg.update(ly, ry, -1);
  }
  return ans;
}

#endif