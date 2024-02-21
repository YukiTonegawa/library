#ifndef _RANGE_UNION_FIND_H_
#define _RANGE_UNION_FIND_H_
#include <vector>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <array>
#include "../range_query/set_segments.hpp"

template<typename Idx>
struct range_union_find{
private:
  set_segments<Idx, false> segment;
  std::unordered_map<Idx, std::pair<Idx, 2>> info; // {区間の左端の座標, {親, サイズ}}
  Idx __find(Idx u){
    std::pair<Idx, Idx> &ui = info.at(u);
    if(ui.first == u) return u;
    return ui.first = __find(ui.first);
  }
  void isolated_check(Idx u){
    auto [f, l, r] = segment.find(u);
    if(!f){
      segment.merge(u, u + 1);
      info.emplace(u, {u, 1});
    }
  }
  void __unite(Idx u, Idx v){
    u = __find(u);
    v = __find(v);
    if(u == v) return;
    std::pair<Idx, Idx> &ui = info.at(u), &vi = info.at(v);
    if(vi.second > ui.second) vi.second += ui.second, ui.first = v;
    else ui.second += vi[0], vi.first = u;
  }
public:
  range_union_find(){}
  // uの親
  Idx find(Idx u){
    auto [f, l, r] = segment.find(u);
    if(!f) return u;
    return __find(l);
  }
  // uを含む連結成分の要素数
  Idx size(Idx u){
    auto [f, l, r] = segment.find(u);
    if(!f) return 1;
    return info.at(__find(l)).second;
  }
  // u, vが同じ連結成分か
  bool same(Idx u, Idx v){
    return find(u) == find(v);
  }
  // u, vを同じ連結成分にする
  void unite(Idx u, Idx v){
    isolated_check(u);
    isolated_check(v);
    u = find(u);
    v = find(v);
    if(u == v) return;
    std::pair<Idx, Idx> &ui = info.at(u), &vi = info.at(v);
    if(vi.second > ui.second) vi.second += ui.second, ui.first = v;
    else ui.second += vi[0], vi.first = u;
  }
  // [l, r)を同じ連結成分にする
  void range_unite(Idx l, Idx r){
    isolated_check(l);
    Idx lp = find(l);
    Idx S = 0; // すでに存在し, [l, r)と被っている要素のサイズ
    for(auto [L, R] : segment.enumerate_intersect(l, r)){
      S += std::min(r, R) - std::max(l, L);
      __unite(lp, L);
    }
    lp = __find(lp);
    std::pair<Idx, Idx> &li = info.at(lp);
    li.second += (r - l) - S;
    segment.merge(l, r);
  }
};
#endif