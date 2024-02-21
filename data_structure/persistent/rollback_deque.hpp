#ifndef _ROLLBACK_DEQUE_H_
#define _ROLLBACK_DEQUE_H_
#include <vector>
#include <deque>
#include <stack>

template<typename T>
struct rollback_deque{
private:
  static constexpr T z = 0;// Tの零元
  std::deque<T> v;
  // {クエリタイプ, 以前の値}
  // 0以上: 書き込み
  // -1: bookmark
  // -2: pop_front
  // -3: push_front
  // -4: pop_back
  // -5: push_back
  std::stack<std::pair<int, T>> h;
public:
  rollback_deque(){}
  rollback_deque(int n, T val): v(n, val){}
  rollback_deque(const std::vector<T> &_v) : v(_v){}
  rollback_deque(const std::deque<T> &_v) : v(_v){}

  int size(){
    return v.size();
  }
  T front(){
    return v.front();
  }
  T back(){
    return v.back();
  }
  void set(int k, T x){
    h.push({k, v[k]});
    v[k] = x;
  }
  T get(int k){
    return v[k];
  }
  void bookmark(){
    h.push({-1, 0});
  }
  void pop_front(){
    h.push({-2, v.front()});
    v.pop_front();
  }
  void push_front(T x){
    h.push({-3, z});
    v.push_front(x);
  }
  void pop_back(){
    h.push({-4, v.back()});
    v.pop_back();
  }
  void push_back(T x){
    h.push({-5, z});
    v.push_back(x);
  }
  void process(int k, T x){
    switch (k){
      case -1: return;
      case -2: v.push_front(x);return;
      case -3: v.pop_front();return;
      case -4: v.push_back(x);return;
      case -5: v.pop_back();return;
      default: v[k] = x;return;
    }
  }
  // 1つ前のbookmarkまで戻して消す
  void rollback_bookmark(){
    while(!h.empty()){
      auto [k, x] = h.top();
      h.pop();
      process(k, x);
      if(k == -1) break;
    }
  }
  // bookmarkを無視して1つ戻す
  void rollback(){
    while(!h.empty()){
      auto [k, x] = h.top();
      h.pop();
      process(k, x);
      if(k != -1) break;
    }
  }
};
#endif