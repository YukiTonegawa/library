#ifndef _MANHATTAN_H_
#define _MANHATTAN_H_
#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>
#include <tuple>
#include "../../data_structure/union_find/union_find.hpp"
#include "../../data_structure/segment_tree/binary_indexed_tree_min.hpp"

template<typename Idx>
struct manhattan_mst{
  struct point{
    Idx x, y;
    int label;
    int cx, cy; // 座標圧縮後のx, y
    point(Idx x, Idx y, int label): x(x), y(y), label(label){}
    Idx dist(const point &r){
      return abs(x - r.x) + abs(y - r.y);
    }
  };
  std::vector<point> p;
private:
  void __compress(){
    int n = p.size();
    std::vector<std::pair<Idx, int>> x(n), y(n);
    for(int i = 0; i < n; i++){
      x[i] = {p[i].x, i};
      y[i] = {p[i].y, i};
    }
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    for(int i = 0; i < n; i++){
      p[x[i].second].cx = i;
      p[y[i].second].cy = i;
    }
  }
  static constexpr std::pair<Idx, int> __id(){
    return {std::numeric_limits<Idx>::max(), -1};
  }
  std::vector<std::pair<int, int>> __nearest(){
    std::vector<std::pair<int, int>> E; // {s, t}
    __compress();
    int n = p.size();
    binary_indexed_tree_min<std::pair<Idx, int>, __id> seg_x(n), seg_y(n);

    std::sort(p.begin(), p.end(), [](const point &a, const point &b){
      return a.x + a.y < b.x + b.y;
    });
    // x + yの昇順
    for(int i = 0; i < n; i++){
      const point &pcur = p[i];
      // 3, x + yがpcurより小さく, xがpcur以上
      auto [dist, id] = seg_x.query(n - 1 - pcur.cx);
      if(id != -1) E.push_back({pcur.label, id});

      // 6, x + yがpcurより小さく, yがpcur以上
      std::tie(dist, id) = seg_y.query(n - 1 - pcur.cy);
      if(id != -1) E.push_back({pcur.label, id});

      seg_x.update(n - 1 - pcur.cx, {pcur.x - pcur.y, pcur.label});
      seg_y.update(n - 1 - pcur.cy, {-pcur.x + pcur.y, pcur.label});
    }
    seg_x.reset();
    seg_y.reset();
    // x + yの降順
    for(int i = n - 1; i >= 0; i--){
      const point &pcur = p[i];
      // 7, x + yがpcurより大きく, xがpcur以下
      auto [dist, id] = seg_x.query(pcur.cx);
      if(id != -1) E.push_back({pcur.label, id});

      // 2, x + yがpcurより大きく, yがpcur以下
      std::tie(dist, id) = seg_y.query(pcur.cy);
      if(id != -1) E.push_back({pcur.label, id});

      seg_x.update(pcur.cx, {-pcur.x + pcur.y, pcur.label});
      seg_y.update(pcur.cy, {pcur.x - pcur.y, pcur.label});
    }
    std::sort(p.begin(), p.end(), [](const point &a, const point &b){
      return a.x - a.y < b.x - b.y;
    });
    seg_x.reset();
    seg_y.reset();
    // x - yの昇順
    for(int i = 0; i < n; i++){
      const point &pcur = p[i];
      // 0, x - yがpcurより小さく, xがpcur以上
      auto [dist, id] = seg_x.query(n - 1 - pcur.cx);
      if(id != -1) E.push_back({pcur.label, id});

      // 5, x - yがpcurより小さく, yがpcur以下
      std::tie(dist, id) = seg_y.query(pcur.cy);
      if(id != -1) E.push_back({pcur.label, id});

      seg_x.update(n - 1 - pcur.cx, {pcur.x + pcur.y, pcur.label});
      seg_y.update(pcur.cy, {-pcur.x - pcur.y, pcur.label});
    }
    seg_x.reset();
    seg_y.reset();
    // x - yの降順
    for(int i = n - 1; i >= 0; i--){
      const point &pcur = p[i];
      // 4, x - yがpcurより大きく, xがpcur以下
      auto [dist, id] = seg_x.query(pcur.cx);
      if(id != -1) E.push_back({pcur.label, id});

      // 1, x + yがpcurより大きく, yがpcur以上
      std::tie(dist, id) = seg_y.query(n - 1 - pcur.cy);
      if(id != -1) E.push_back({pcur.label, id});

      seg_x.update(pcur.cx, {-pcur.x - pcur.y, pcur.label});
      seg_y.update(n - 1 - pcur.cy, {pcur.x + pcur.y, pcur.label});
    }
    std::sort(p.begin(), p.end(), [&](const point &a, const point &b){return a.label < b.label;});
    return E;
  }
  std::vector<std::pair<int, int>> __mst(){
    int n = p.size();
    auto E = __nearest();
    union_find uf(n);
    std::sort(E.begin(), E.end(), [&](std::pair<int, int> a, std::pair<int, int> b){
      return p[a.first].dist(p[a.second]) < p[b.first].dist(p[b.second]);
    });
    std::vector<std::pair<int, int>> res;
    for(auto [a, b] : E){
      if(!uf.same(a, b)){
        res.push_back({a, b});
        uf.unite(a, b);
      }
    }
    return res;
  }
public:
  manhattan_mst(){}
  void add_point(Idx x, Idx y, int label){
    p.push_back(point(x, y, label));
  }
  /*      y
          ┃
        7 ┃ 0
      6   ┃   1
  ━━━━━━━━━━━━━━━━━━  x
      5   ┃   2
       4  ┃ 3
          ┃
  */
  // 各点に対して8方向の最近点をそれぞれ求め, 合計8n個以下の点を返す. O(nlog(n))
  std::vector<std::pair<int, int>> nearest(){
    return __nearest();
  }
  // 最小全域木を構成する辺 {a, b}, a, bは点のラベル
  std::vector<std::pair<int, int>> mst(){
    return __mst();
  }
  // pのi番目の点(labelがiとは限らない)
  point& operator [](int i){return p[i];}
};

template<typename Idx>
std::pair<Idx, Idx> manhattan_rotate(Idx x, Idx y){
  return {x + y, x - y};
}


template<typename Idx, int Dim>
struct chebyshev{
  static constexpr int manhattan_dim = Dim;
  static constexpr int chebyshev_dim = 1 << (Dim - 1);
  using manhattan_point = std::array<Idx, manhattan_dim>;
  using chebyshev_point = std::array<Idx, chebyshev_dim>;

  chebyshev_point manhattan2chebyshev(const manhattan_point &p){
    if(Dim == 2) return {p[0] + p[1], p[0] - p[1]};
    if(Dim == 3) return {p[0] + p[1] + p[2], p[0] + p[1] - p[2], p[0] - p[1] + p[2], p[0] - p[1] - p[2]};
    chebyshev_point res;
    res.fill(p[0]);
    for(int i = 0; i < (1 << chebyshev_dim); i++){
      for(int j = 0; j < manhattan_dim - 1; j++){
        if((i >> j) & 1) res[i] -= p[manhattan_dim - 1 - j];
        else res[i] += p[manhattan_dim - 1 - j];
      }
    }
    return res;
  }
  std::array<std::pair<Idx, Idx>, chebyshev_dim> max_min;
  chebyshev(){
    for(int i = 0; i < chebyshev_dim; i++){
      max_min[i] = {std::numeric_limits<Idx>::min(), std::numeric_limits<Idx>::max()};
    }
  }
  void add_point(const manhattan_point &p){
    chebyshev_point q = manhattan2chebyshev(p);
    for(int i = 0; i < chebyshev_dim; i++){
      max_min[i].first = std::max(max_min[i].first, q[i]);
      max_min[i].second = std::min(max_min[i].second, q[i]);
    }
  }
  Idx max_dist(const manhattan_point &p){
    assert(!max_min.empty());
    chebyshev_point q = manhattan2chebyshev(p);
    Idx res = std::numeric_limits<Idx>::min();
    for(int i = 0; i < chebyshev_dim; i++){
      res = std::max(res, abs(max_min[i].first - q[i]));
      res = std::max(res, abs(max_min[i].second - q[i]));
    }
    return res;
  }
};
#endif