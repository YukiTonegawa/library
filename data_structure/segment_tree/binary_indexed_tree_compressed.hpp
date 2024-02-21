#ifndef _BINARY_INDEXED_TREE_COMPRESSED_H_
#define _BINARY_INDEXED_TREE_COMPRESSED_H_
#include <vector>
#include <iostream>
template<typename Idx, typename Val>
struct binary_indexed_tree_compressed{
  int M;
  std::vector<std::pair<Val, Val>> sum;
  binary_indexed_tree_compressed(){}
  binary_indexed_tree_compressed(int N): M(N), sum(N + 1, {0, 0}){}
  binary_indexed_tree_compressed(const std::vector<Val> &v): M(v.size()), sum(M + 1, {0, 0}){
    for(int i = 0; i < M; i++) sum[i + 1].first = v[i];
    for(int i = 1; i <= M; i++){
      int nxt = i + (i & (-i));
      if(nxt <= M) sum[nxt].first += sum[i].first;
    }
  }
  void update(int kc, Val x){
    for(int i = kc + 1; i <= M; i += (i & (-i))) sum[i].first += x;
  }
  // [lc, rc) (元の数直線上では[l, r))にzを加算
  void update(Idx l, Idx r, int lc, int rc, Val z){
    Val a = l * z, b = r * z;
    for(int i = lc + 1; i <= M; i += (i & (-i))){
      sum[i].first -= a;
      sum[i].second += z;
    }
    for(int i = rc + 1; i <= M; i += (i & (-i))){
      sum[i].first += b;
      sum[i].second -= z;
    }
  }
  Val query(Idx r, int rc){
    Val a = 0, b = 0;
    for(int i = rc; i > 0; i -= (i & (-i))){
      a += sum[i].first;
      b += sum[i].second;
    }
    return a + (b * r);
  }
  // [lc, rc) (元の数直線上では[l, r))の和
  Val query(Idx l, Idx r, int lc, int rc){
    if(lc >= rc) return 0;
    return query(r, rc) - query(l, lc);
  }
};
#endif