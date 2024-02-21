#ifndef _COMPRESS_SUBSET_SUM_H_
#define _COMPRESS_SUBSET_SUM_H_
#include <vector>
#include "../data_structure/basic/binary_heap.hpp"

// Sを作れるか判定する場合, Sより大きいものは捨てられる
// 同じものが3つ以上ある場合, 2つまとめていい
template<typename T>
std::vector<T> compress_subset_sum(const std::vector<T> &a){
  binary_heap<T> h(a);
  std::vector<T> ret;
  while(!h.empty()){
    T x = h.pop_min();
    if(h.empty() || h.min() > x){
      ret.push_back(x);
    }else if(h.min() == x){
      h.pop_min();
      if(h.empty() || h.min() != x){
        h.push(x);
        ret.push_back(x);
      }else{
        h.push(2 * x);
      }
    }
  }
  return ret;
}
#endif
