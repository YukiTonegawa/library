#ifndef _SIMPLE_DATA_STRUCTURE_H_
#define _SIMPLE_DATA_STRUCTURE_H_
#include <algorithm>
#include <array>

template<typename Val>
struct top2{
  Val a = -1, b = -1;
  inline void push(Val x){
    if(a < x) std::swap(a, x);
    b = max(b, x);
  }
  inline Val get(){
    return a;
  }
};
template<typename Val, typename T>
struct top2p{
  Val a = -1, b = -1;
  T ai = -1, bi = -1;
  inline void push(Val x, T xi){
    if(a < x){
      std::swap(a, x);
      std::swap(ai, xi);
    }
    if(b < x){
      std::swap(b, x);
      std::swap(bi, xi);
    }
  }
  inline std::pair<Val, T> get(){
    return {a, ai};
  }
};
template<typename Val, int k>
struct topk{
  std::array<Val, k> v;
  topk(){v.fill(-1);}
  inline void push(Val x){
    for(int i = 0; i < k; i++){
      if(v[i] < x) std::swap(v[i], x);
    }
  }
  inline Val get(int i){
    return v[i];
  }
};
template<typename Val, typename T, int k>
struct topkp{
  std::array<Val, k> v;
  std::array<T, k> t;
  topkp(){v.fill(-1);t.fill(-1);}
  inline void push(Val x, T y){
    for(int i = 0; i < k; i++){
      if(v[i] < x){
        std::swap(v[i], x);
        std::swap(t[i], y);
      }
    }
  }
  inline std::pair<Val, T> get(int i){
    return {v[i], t[i]};
  }
};
#endif