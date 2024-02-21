#ifndef _SWAG_H_
#define _SWAG_H_
#include <vector>
#include <deque>
#include "../../algebraic_structure/monoid.hpp"

template<typename monoid>
struct swag{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  std::deque<Val> q;
  std::vector<Val> sum_left, sum_right;
  swag(){
    sum_left.push_back(id());
    sum_right.push_back(id());
  }
  int size(){
    return q.size();
  }
  bool empty(){
    return q.size() == 0;
  }
  // 空の場合は単位元
  Val query_all(){
    return merge(sum_left.back(), sum_right.back());
  }
  void push_back(Val x){
    q.push_back(x);
    sum_right.push_back(merge(sum_right.back(), x));
  }
  void pop_back(){
    assert(size());
    if(sum_right.size() > 1){
      sum_right.pop_back();
      q.pop_back();
      return;
    }
    int new_r_size = (q.size() + 1) / 2;
    int new_l_size = q.size() - new_r_size;
    q.pop_back();
    for(int i = new_l_size; i < q.size(); i++) sum_right.push_back(merge(sum_right.back(), q[i]));
    sum_left.resize(new_l_size + 1);
    for(int i = 0; i < new_l_size; i++) sum_left[i + 1] = merge(q[new_l_size - 1 - i], sum_left[i]);
  }
  void push_front(Val x){
    q.push_front(x);
    sum_left.push_back(merge(x, sum_left.back()));
  }
  void pop_front(){
    assert(size());
    if(sum_left.size() > 1){
      sum_left.pop_back();
      q.pop_front();
      return;
    }
    int new_l_size = (q.size() + 1) / 2;
    int new_r_size = q.size() - new_l_size;
    for(int i = 0; i < new_l_size - 1; i++) sum_left.push_back(merge(q[new_l_size - 1 - i], sum_left[i]));
    sum_right.resize(new_r_size + 1);
    for(int i = 0; i < new_r_size; i++) sum_right[i + 1] = merge(sum_right[i], q[new_l_size + i]);
    q.pop_front();
  }
  Val back(){
    return q.back();
  }
  Val front(){
    return q.front();
  }
};
#endif