#ifndef _ROLLBACK_ARRAY_H_
#define _ROLLBACK_ARRAY_H_
#include <vector>
#include <stack>

template<typename T>
struct rollback_array{
private:
  static constexpr T z = 0;// Tの零元
  std::vector<T> v;
  std::stack<std::pair<int, T>> h;// {書き換えられたインデックス, 以前の値}
public:
  rollback_array(){}
  rollback_array(int n, T val): v(n, val){}
  rollback_array(const std::vector<T> &_v) : v(_v){}

  int size(){
    return v.size();
  }
  void set(int k, T x){
    h.push({k, v[k]});
    v[k] = x;
  }
  T get(int k){
    return v[k];
  }
  void bookmark(){
    h.push({-1, z});
  }
  // 1つ前のbookmarkまで戻して消す
  void rollback_bookmark(){
    while(!h.empty()){
      auto [k, x] = h.top();
      h.pop();
      if(k == -1) break;
      v[k] = x;
    }
  }
  // bookmarkを無視して1つ戻す
  void rollback(){
    while(!h.empty()){
      auto [k, x] = h.top();
      h.pop();
      if(k != -1){
        v[k] = x;
        break;
      }
    }
  }
};
#endif