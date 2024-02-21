#ifndef _SPARSE_BINARY_H_
#define _SPARSE_BINARY_H_
#include <vector>
#include <string>
#include <cassert>
#include "../segment_tree/sparse_lazy_segment_tree.hpp"

template<typename Idx = int>
struct sparse_binary{
  sparse_lazy_segment_tree<range_add_range_sum<Idx>> seg;
  sparse_binary(){}
  sparse_binary(std::string s){
    if(s.empty()) return;
    std::reverse(s.begin(), s.end());
    int n = s.size();
    std::vector<std::pair<Idx, bool>> v;
    for(int i = 0; i < n; i++){
      bool f = (s[i] == '1');
      if(v.empty() || v.back().second != f) v.push_back({i, f});
    }
    v.push_back({n, 0});
    for(int i = 0; i < v.size() - 1; i++){
      auto [l, f] = v[i];
      if(f) seg.update(l, v[i + 1].first, f);
    }
  }
  Idx msb(){
    auto [s, v] = seg.bisect_from_right(seg.inf, [&](Idx x){return x > 0;});
    if(!v) return -1;
    return v->rx - 1;
  }
  Idx lsb(){
    select1(0);
  }
  Idx find_next1(Idx k){
    auto [s, v] = seg.bisect_from_left(k, [&](Idx x){return x > 0;});
    if(!v) return -1;
    return std::max(k, v->lx);
  }
  Idx find_next0(Idx k){
    auto [s, v] = seg.bisect_from_left2(k, [&](Idx x, Idx lx, Idx rx){return (rx - lx) - x > 0;});
    if(!v) return -1;
    return std::max(k, v->lx);
  }
  Idx find_prev1(Idx k){
    auto [s, v] = seg.bisect_from_right(k, [&](Idx x){return x > 0;});
    if(!v) return -1;
    return std::min(k, v->rx - 1);
  }
  Idx find_prev0(Idx k){
    auto [s, v] = seg.bisect_from_right2(k, [&](Idx x, Idx lx, Idx rx){return (rx - lx) - x > 0;});
    if(!v) return -1;
    return std::min(k, v->rx - 1);
  }
  Idx select1(Idx k){
    auto [s, v] = seg.bisect_from_left(0, [&](Idx x){return x > k;});
    if(!v) return -1;
    return std::max((Idx)0, v->lx) + (k - s);
  }
  Idx select0(Idx k){
    auto [s, v] = seg.bisect_from_left2(0, [&](Idx x, Idx lx, Idx rx){return (rx - lx) - x > k;});
    if(!v) return -1;
    return std::max((Idx)0, v->lx) + (k - v->lx + s);
  }
  Idx rank1(Idx r){
    return seg.query(0, r);
  }
  Idx rank0(Idx r){
    return r - rank1(r);
  }
  Idx popcount(){
    return seg.query_all();
  }
  bool is_zero(){
    return popcount() == 0;
  }
  void set(Idx k, bool f){
    seg.set(k, f);
  }
  void set_range(Idx l, Idx r, bool f){
    seg.reset(l, r);
    if(f) seg.update(l, r, 1);
  }
  bool get(Idx k){
    return seg.get(k);
  }
  // add 2^k
  void add(Idx k){
    Idx r = find_next0(k);
    seg.set(r, 1);
    seg.update(k, r, -1);
  }
  // sub 2^k
  void sub(Idx k){
    Idx r = find_next1(k);
    seg.set(r, 0);
    seg.update(k, r, 1);
  }
  std::string to_string(Idx size_limit = 1000000){
    Idx m = std::min(size_limit, msb() + 1);
    std::string ret(m, '.');
    for(auto [l, r, f] : seg.enumerate()){
      if(l >= m) continue;
      l = std::max(0, l);
      r = std::min(r, m);
      for(int j = l; j < r; j++){
        ret[j] = (f ? '1' : '0');
      }
    }
    std::reverse(ret.begin(), ret.end());
    return ret;
  }
};
#endif
