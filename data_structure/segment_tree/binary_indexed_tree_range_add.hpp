#ifndef _BINARY_INDEXED_TREE_RANGE_ADD_H_
#define _BINARY_INDEXED_TREE_RANGE_ADD_H_
#include <vector>
#include <iostream>

template<typename Val>
struct binary_indexed_tree_range_add{
  int M;
  std::vector<std::pair<Val, Val>> sum;
  binary_indexed_tree_range_add(){}
  binary_indexed_tree_range_add(int N): M(N), sum(M + 1, {0, 0}){}
  binary_indexed_tree_range_add(const std::vector<Val> &v): M(v.size()), sum(M + 1, {0, 0}){
    for(int i = 0; i < M; i++) sum[i + 1].first = v[i];
    for(int i = 1; i <= M; i++){
      int nxt = i + (i & (-i));
      if(nxt <= M) sum[nxt].first += sum[i].first;
    }
  }
  void update(int k, Val x){
    for(int i = k + 1; i <= M; i += (i & (-i))) sum[i].first += x;
  }
  void update(int l, int r, Val x){
    Val a = x * l, b = x * r;
    for(int i = l + 1; i <= M;i += (i & (-i))){
      sum[i].first -= a;
      sum[i].second += x;
    }
    for(int i = r + 1; i <= M; i += (i & (-i))){
      sum[i].first += b;
      sum[i].second -= x;
    }
  }
  Val query(int r){
    Val a = 0, b = 0;
    for(int i = r; i > 0; i -= (i & (-i))){
      a += sum[i].first;
      b += sum[i].second;
    }
    return a + (b * r);
  }
  Val query(int l, int r){
    if(l >= r) return 0;
    return query(r) - query(l);
  }
};
#endif