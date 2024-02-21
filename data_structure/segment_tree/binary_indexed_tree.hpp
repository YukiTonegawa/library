#ifndef _BINARY_INDEXED_TREE_H_
#define _BINARY_INDEXED_TREE_H_
#include <vector>
#include <cassert>
template<typename Val>
struct binary_indexed_tree{
  int M, H;
  std::vector<Val> sum;
  binary_indexed_tree(){}
  binary_indexed_tree(int N): M(N), H(31 - __builtin_clz(M)), sum(M + 1 , 0){}
  binary_indexed_tree(const std::vector<Val> &v): M(v.size()), H(31 - __builtin_clz(M)), sum(1){
    sum.insert(sum.begin() + 1, v.begin(), v.end());
    for(int i = 1; i <= M; i++){
      int nxt = i + (i & (-i));
      if(nxt <= M) sum[nxt] += sum[i];
    }
  }
  void update(int k, Val x){
    for(int i = k + 1; i <= M; i += (i & (-i))) sum[i] += x;
  }
  Val query(int r){
    Val ret = 0;
    for(int k = r; k > 0; k -= (k & (-k))) ret += sum[k];
    return ret;
  }
  Val query(int l, int r){
    return query(r) - query(l);
  }
  // sum[0, k]がx以上になるような最小のkとsum[0, k], 無い場合は{M, 総和}
  // sumが単調非減少であることが必要
  using p = std::pair<int, Val>;
  p lower_bound(Val x){
    int v = 1 << H, h = H;
    Val s = 0, t = 0;
    while(h--){
      if(M < v) v -= 1 << h;
      else if(x <= s + sum[v]) t = s + sum[v], v -= 1 << h;
      else s += sum[v], v += 1 << h;
    }
    if(v == M + 1) return {M, s};
    return (x <= s + sum[v] ? p{v - 1, s + sum[v]} : p{v, t});
  }
};
#endif