#ifndef _ROLLBACK_DEQUE_US_H_
#define _ROLLBACK_DEQUE_US_H_
#include <vector>

template<typename T, int qmax>
struct stack_unsafe{
  int s = 0;
  T v[qmax];
  inline int size(){
    return s;
  }
  inline void push(T x){
    v[s++] = x;
  }
  inline void pop(){
    s--;
  }
  inline T top(){
    return v[s - 1];
  }
  void reset(){
    s = 0;
  }
};
// qmax := 行われる操作の上限
template<typename T, int qmax>
struct rollback_deque_unsafe{
  static constexpr T z = 0;// Tの零元
  stack_unsafe<std::pair<int, T>, qmax> h;
  T v[qmax * 2 + 1];
  int l = qmax, r = qmax;

  inline int size(){
    return r - l;
  }
  inline T front(){
    return v[l];
  }
  inline T back(){
    return v[r - 1];
  }
  inline void set(int k, T x){
    h.push({k, v[l + k]});
    v[l + k] = x;
  }
  inline T get(int k){
    return v[l + k];
  }
  inline void bookmark(){
    h.push({-1, 0});
  }
  inline void pop_front(){
    h.push({-2, v[l++]});
  }
  inline void push_front(T x){
    h.push({-3, z});
    v[--l] = x;
  }
  inline void pop_back(){
    h.push({-4, v[--r]});
  }
  inline void push_back(T x){
    h.push({-5, z});
    v[r++] = x;
  }
  inline void process(int k, T x){
    switch (k){
      case -1: return;
      case -2: v[k + (--l)] = x;return;
      case -3: l++;return;
      case -4: v[k + (r++)] = x;return;
      case -5: r--;return;
      default: v[k + l] = x;return;
    }
  }
  // 1つ前のbookmarkまで戻して消す
  inline void rollback_bookmark(){
    while(h.size()){
      auto [k, x] = h.top();
      h.pop();
      process(k, x);
      if(k == -1) break;
    }
  }
  // bookmarkを無視して1つ戻す
  inline void rollback(){
    while(h.size()){
      auto [k, x] = h.top();
      h.pop();
      process(k, x);
      if(k != -1) break;
    }
  }
  void reset(){
    l = r = qmax;
    h.reset();
  }
};
#endif