#ifndef _PERSISTENT_QUEUE_H_
#define _PERSISTENT_QUEUE_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include "../../graph/dynamic_graph/incremental_tree.hpp"
// 空のqueueをバージョン0とする
template<typename T>
struct persistent_queue{
  incremental_tree t;
  std::vector<T> val;
  std::vector<int> id{0}, len{0};
  persistent_queue(){}
  void push(int prev, T x){
    val.push_back(x);
    id.push_back(t.make_node(id[prev]));
    len.push_back(len[prev] + 1);
  }
  T pop(int prev){
    assert(len[prev] > 0);
    id.push_back(id[prev]);
    len.push_back(len[prev] - 1);
    int v = t.la(id[prev], len[prev] - 1);
    assert(v);
    return val[v - 1];
  }
  T front(int v){
    assert(len[v] > 0);
    v = t.la(id[v], len[v] - 1);
    assert(v);
    return val[v - 1];
  }
  T access(int v, int k){
    assert(len[v] > k);
    v = t.la(id[v], len[v] - k - 1);
    assert(v);
    return val[v - 1];
  }
};
#endif