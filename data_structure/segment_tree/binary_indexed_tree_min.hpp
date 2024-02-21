#ifndef _BINARY_INDEXED_TREE_MIN_H_
#define _BINARY_INDEXED_TREE_MIN_H_
#include <vector>
#include <algorithm>

template<typename Val, Val (*id)()>
struct binary_indexed_tree_min{
  int M;
  std::vector<Val> sum;
  binary_indexed_tree_min(int N): M(N), sum(M + 1, id()){assert(N > 0);}
  void reset(){
    std::fill(sum.begin(), sum.end(), id());
  }
  void update(int k, Val x){
    for(int i = k + 1; i <= M; i += (i & (-i))){
      if(sum[i] > x) sum[i] = x;
      else return;
    }
  }
  // min([0, r))
  Val query(int r){
    Val res = id();
    for(int k = r; k > 0; k -= (k & (-k))) res = std::min(res, sum[k]);
    return res;
  }
};
#endif