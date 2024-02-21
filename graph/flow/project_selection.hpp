#ifndef _PROJECT_SELECTION_H_
#define _PROJECT_SELECTION_H_

#include "max_flow.hpp"

// https://ei1333.github.io/luzhiled/snippets/memo/project-selection.html
template<typename Cost>
struct project_selection{
  int n, s, t;
  mf_graph<Cost> mf;
  Cost offset;
  project_selection(int _n): n(_n), s(_n), t(_n + 1), mf(_n + 2), offset(0){}
  // iがfのとき+x
  // x < 0でもok
  void add_if(int i, bool f, Cost x){
    assert(0 <= i && i < n);
    if(x < 0){
      if(f) mf.add_edge(s, i, -x);
      else mf.add_edge(i, t, -x);
    }
    if(x > 0){
      offset += x;
      if(f) mf.add_edge(i, t, x);
      else mf.add_edge(s, i, x);
    }
  }
  // i, jが1のとき+x
  // x >= 0
  void add_and(int i, int j, Cost x){
    assert(x >= 0);
    assert(0 <= i && i < n);
    assert(0 <= j && j < n);
    if(x == 0) return;
    offset += x;
    int k = mf.add_vertex();
    mf.add_edge(k, t, x);
    mf.add_edge(i, k, x);
    mf.add_edge(j, k, x);
  }
  // i, jが0のとき+x
  // x >= 0
  void add_nor(int i, int j, Cost x){
    assert(x >= 0);
    assert(0 <= i && i < n);
    assert(0 <= j && j < n);
    if(x == 0) return;
    offset += x;
    int k = mf.add_vertex();
    mf.add_edge(s, k, x);
    mf.add_edge(k, i, x);
    mf.add_edge(k, j, x);
  }
  Cost max_score(){
    return offset - mf.flow(s, t);
  }
};
#endif
