#ifndef _BEAATS_ADD_SUM_HSUM_H_
#define _BEAATS_ADD_SUM_HSUM_H_
#include "../binary_indexed_tree_range_add.hpp"

// 区間加算-区間historic_sumの和　
// O(logN)
template<typename Val, typename ValSum>
struct beats_add_sum_hsum{
  int qid;
  binary_indexed_tree_range_add<ValSum> hsum, bit;
  beats_add_sum_hsum(const std::vector<ValSum> &v): qid(1), hsum(v.size()), bit(v){}
  void update_add(int l, int r, Val x){
    hsum.update(l, r, -(ValSum)qid * x);
    bit.update(l, r, x);
    qid++;
  }
  Val get(int k){
    return bit.query(k, k + 1);
  }
  ValSum query_sum(int l, int r){
    return bit.query(l, r);
  }
  Val get_hmin(int k){
    return hsum.get(k) + (ValSum)qid * bit.query(k, k + 1);
  }
  ValSum query_sum_hmin(int l, int r){
    return hsum.query_sum(l, r) + (ValSum)qid * bit.query(l, r);
  }
};

#endif