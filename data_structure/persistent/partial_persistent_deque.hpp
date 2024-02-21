#ifndef _PARTIAL_PERSISTENT_DEQUE_H_
#define _PARTIAL_PERSISTENT_DEQUE_H_
#include "partial_persistent_variable.hpp"
#include "partial_persistent_array.hpp"
#include <vector>
#include <cassert>
#include <iostream>

template<typename Val, Val error_val = 0>
struct partial_persistent_deque{
private:
  partial_persistent_variable<int> l, r;
  int buffer;
  partial_persistent_array<Val> v;
public:
  partial_persistent_deque(){}
  partial_persistent_deque(int _buffer): l(0), r(0), buffer(_buffer + 1), v(_buffer + 1, error_val){}

  int size(int t){
    int r_t = r.get(t), l_t = l.get(t);
    return r_t < l_t ? r_t - l_t + buffer : r_t - l_t;
  }
  int size_new(){
    int r_t = r.get_new(), l_t = l.get_new();
    return r_t < l_t ? r_t - l_t + buffer : r_t - l_t;
  }
  int time(){
    return v.time();
  }
  void push_back(Val x){
    int r_new = r.get_new();
    assert(size_new() < buffer - 1);
    v.next();
    v.set(r_new, x);
    r.set(v.time(), (r_new + 1 == buffer ? 0 : r_new + 1));
  }
  void pop_back(){
    int r_new = r.get_new();
    assert(size_new());
    v.next();
    r.set(v.time(), (r_new == 0 ? buffer - 1 : r_new - 1));
  }
  Val back(int t){
    assert(size(t));
    int r_t = r.get(t);
    r_t = (r_t == 0 ? buffer - 1 : r_t - 1);
    return v.get(t, r_t);
  }
  Val back_new(){
    assert(size_new());
    int r_t = r.get_new();
    r_t = (r_t == 0 ? buffer - 1 : r_t - 1);
    return v.get_new(r_t);
  }
  void push_front(Val x){
    int l_new = l.get_new();
    assert(size_new() < buffer - 1);
    v.next();
    l_new = (l_new == 0 ? buffer - 1 : l_new - 1);
    v.set(l_new, x);
    l.set(v.time(), l_new);
  }
  void pop_front(){
    int l_new = l.get_new();
    assert(size_new());
    v.next();
    l.set(v.time(), (l_new + 1 == buffer ? 0 : l_new + 1));
  }
  Val front(int t){
    assert(size(t));
    return v.get(t, l.get(t));
  }
  Val front_new(){
    assert(size_new());
    return v.get_new(l.get_new());
  }
  Val get(int t, int k){
    int l_t = l.get(t);
    assert(size(t) > k);
    return v.get(t, (l_t + k >= buffer ? l_t + k - buffer : l_t + k));
  }
  Val get_new(int k){
    int l_t = l.get_new();
    assert(size_new() > k);
    return v.get_new((l_t + k >= buffer ? l_t + k - buffer : l_t + k));
  }
  void set_new(int k, Val x){
    int l_t = l.get_new();
    assert(size_new() > k);
    v.next();
    v.set((l_t + k >= buffer ? l_t + k - buffer : l_t + k), x);
  }
  void rollback(int t){
    l.rollback(t);
    r.rollback(t);
    v.rollback(t);
  }
  void print(int t){
    int s_t = size(t);
    std::cout << "time: " << t << ", size: " << s_t << '\n';
    for(int i = 0; i < s_t; i++){
      if(i == s_t - 1) std::cout << get(t, i) << '\n';
      else std::cout << get(t, i) << ' ';
    }
  }
};
#endif