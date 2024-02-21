#ifndef _BEAATS_ADD_MAX_HMIN_H_
#define _BEAATS_ADD_MAX_HMIN_H_

#include "chmin_add_maxsum.hpp"
// 区間加算-区間historic_minのmax
// ならしO(log^2N)
template<typename Val>
struct beats_add_max_hmin{
  beats_chmin_add_maxsum<Val> seg;
  beats_add_max_hmin(const std::vector<Val> &v){
    int n = v.size();
    std::vector<std::pair<Val, Val>> tmp(n);
    for(int i = 0; i < n; i++) tmp[i] = {0, v[i]};
    seg = beats_chmin_add_maxsum<Val>(tmp);
  }
  void update_add(int l, int r, Val x){
    seg.update_add(l, r, -x, x);
    if(x < 0) seg.update_chmin_x(l, r, 0);
  }
  Val get(int k){
    return seg.get(k).second;
  }
  Val get_hmin(int k){
    auto tmp = seg.get(k);
    return tmp.first + tmp.second;
  }
  Val query_max_hmin(int l, int r){
    return seg.query_maxsum(l, r);
  }
};

#endif