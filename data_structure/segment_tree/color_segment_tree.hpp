#ifndef _COLOR_SEGMENT_TREE_H_
#define _COLOR_SEGMENT_TREE_H_
#include <vector>
#include <cassert>
#include "../range_query/compress_bbst.hpp"
template<typename monoid>
struct color_segment_tree{
  using Color = int;
  using Val = typename monoid::Val;
  using Cmp = compress_bbst<Color>;
  using St = splay_tree_monoid<int, monoid>;
  typename Cmp::node *root = nullptr;
  std::unordered_map<Color, typename St::node*> mp;
  color_segment_tree(const std::vector<std::pair<Color, Val>> &v){
    for(int i = 0; i < v.size(); i++){
      auto [c, x] = v[i];
      auto a = Cmp::make_node(c, c);
      auto b = new typename St::node(i, x);
      root = Cmp::merge(root, a);
      auto itr = mp.find(c);
      if(itr == mp.end()) mp.emplace(c, b);
      else itr->second = St::merge(itr->second, b);
    }
  }
  // k番目の要素の色をcに, 値をvalにする
  void set(int k, Color c, Val val){
    auto [x, y, z] = Cmp::split3(root, k, k + 1);
    auto itr = mp.find(y->cmp);
    y->cmp = c;
    root = Cmp::merge3(x, y, z);
    itr->second = St::erase(itr->second, k);
    itr = mp.find(c);
    if(itr == mp.end()){
      mp.emplace(new typename St::node(k, val));
    }else{
      itr->second = St::insert(itr->second, k, val, 0);
    }
  }
  // [l, r)の要素の色をcにする
  void set_color(int l, int r, Color c){
    if(l >= r) return;
    auto [x, y, z] = Cmp::split3(root, l, r);
    std::vector<std::tuple<int, int, Color>> cmp;
    y = Cmp::compress2(y, cmp);
    y->cmp = c;
    root = Cmp::merge3(x, y, z);
    typename St::node *v = nullptr;
    for(auto [L, R, C] : cmp){
      auto itr = mp.find(C);
      auto [X, Y, Z] = St::split3(itr->second, L, R + 1);
      itr->second = St::merge(X, Z);
      v = St::merge(v, Y);
    }
    auto itr = mp.find(c);
    if(itr == mp.end()){
      mp.emplace(c, v);
    }else{
      auto [X, Y] = St::split(itr->second, l);
      itr->second = St::merge(St::merge(X, v), Y);
    }
  }
  // [l, r)の色cの要素(を順序を保って並べた時)のモノイド積
  Val query(int l, int r, Color c){
    auto itr = mp.find(c);
    if(itr == mp.end() || l >= r) return monoid::id();
    auto [x, y, z] = St::split3(itr->second, l, r);
    Val ans = (y ? y->sum : monoid::id());
    itr->second = St::merge3(x, y, z);
    return ans;
  }
};
#endif