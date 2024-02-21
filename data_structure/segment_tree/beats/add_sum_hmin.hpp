#ifndef _BEAATS_ADD_SUM_HMIN_H_
#define _BEAATS_ADD_SUM_HMIN_H_
#include "chmin_add.hpp"
#include "../binary_indexed_tree_range_add.hpp"

// 区間加算-区間historic_minの和　
// ならしO(log^2N)
template<typename Val, typename ValSum>
struct beats_add_sum_hmin{
  beats_chmin_add<Val, ValSum> seg;
  binary_indexed_tree_range_add<ValSum> bit;
  beats_add_sum_hmin(const std::vector<ValSum> &v): seg(v.size(), 0), bit(v){}
  void update_add(int l, int r, Val x){
    seg.update_add(l, r, -x);
    bit.update(l, r, x);
    if(x < 0) seg.update_chmin(l, r, 0);
  }
  Val get(int k){
    return bit.query(k, k + 1);
  }
  ValSum query_sum(int l, int r){
    return bit.query(l, r);
  }
  Val get_hmin(int k){
    return seg.get(k) + bit.query(k, k + 1);
  }
  ValSum query_sum_hmin(int l, int r){
    return seg.query_sum(l, r) + bit.query(l, r);
  }
};

#endif