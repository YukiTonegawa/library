#ifndef _PERSISTENT_STACK_H_
#define _PERSISTENT_STACK_H_
#include <vector>
#include <algorithm>
#include <cassert>
// 空のqueueをバージョン0とする
template<typename T>
struct persistent_stack{
  std::vector<T> val;
  std::vector<int> id{0}, len{0}, par{0};
  persistent_stack(){}
  void push(int prev, T x){
    val.push_back(x);
    id.push_back(par.size());
    par.push_back(id[prev]);
    len.push_back(len[prev] + 1);
  }
  T pop(int prev){
    assert(len[prev] > 0);
    id.push_back(par[id[prev]]);
    len.push_back(len[prev] - 1);
    return val[id[prev] - 1];
  }
  T back(int v){
    assert(len[v] > 0);
    return val[id[prev] - 1];
  }
};

#endif