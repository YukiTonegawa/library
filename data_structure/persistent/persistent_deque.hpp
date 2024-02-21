#ifndef _PERSISTENT_QUEUE_H_
#define _PERSISTENT_QUEUE_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include "../../graph/dynamic_graph/incremental_tree.hpp"

template<typename T>
struct persistent_deque{
  incremental_tree t;
  std::vector<T> val;
  std::vector<int> first{0}, last{0}, len{0};
  persistent_deque(){}
  void push_back(int prev, T x){
    val.push_back(x);
    first.push_back(first[prev]);
    last.push_back(t.make_node(last[prev]));
    len.push_back(len[prev] + 1);
  }
  void push_front(int prev, T x){
    val.push_back(x);
    first.push_back(t.make_node(first[prev]));
    last.push_back(last[prev]);
    len.push_back(len[prev] + 1);
  }
  T pop_front(int prev){
    assert(len[prev]);
    if(len[prev] == 1){      
      first.push_back(0);
      last.push_back(0);
      len.push_back(0);
      return first[prev] ? val[first[prev] - 1] : val[last[prev] - 1];
    }
    int l = t.lca(first[prev], last[prev]);
    last.push_back(last[prev]);
    len.push_back(len[prev] - 1);
    if(first[prev] != l){
      // first != lかつfirstが根でない
      first.push_back(t.parent(prev));
      return val[first[prev] - 1];
    }else{
      first.push_back(t.la(last[prev], len[prev] - 2));
      return val[t.parent(first.back()) - 1];
    }
  }
  T pop_back(int prev){
    assert(len[prev]);
    if(len[prev] == 1){      
      first.push_back(0);
      last.push_back(0);
      len.push_back(0);
      return first[prev] ? val[first[prev] - 1] : val[last[prev] - 1];
    }
    int l = t.lca(first[prev], last[prev]);
    first.push_back(first[prev]);
    len.push_back(len[prev] - 1);
    if(last[prev] != l){
      // last != lかつlastが根でない
      last.push_back(t.parent(prev));
      return val[last[prev] - 1];
    }else{
      last.push_back(t.la(first[prev], len[prev] - 2));
      return val[t.parent(last.back()) - 1];
    }
  }
  T back(int v){
    assert(len[v]);
    if(last[v]) return val[last[v]];
    else return val[t.la(first[v], len[v] - 1) - 1];
  }
  T front(int v){
    assert(len[v]);
    if(first[v]) return val[first[v]];
    else return val[t.la(last[v], len[v] - 1) - 1];
  }
  T access(int v, int k){
    assert(0 <= k && k < len[v]);
    int l = t.lca(first[v], last[v]);
    if(k < t.depth(first[v]) - t.depth(l)) return val[t.la(first[v], k) - 1];
    k = len[v] - 1 - k;
    return val[t.la(last[v], k) - 1];
  }
};
#endif