#ifndef _PARALLEL_UNION_FIND_H_
#define _PARALLEL_UNION_FIND_H_
#include <vector>
#include <numeric>
#include <cassert>

struct parallel_union_find{
private:
  int N, H;
  std::vector<int> sz, par;
  void push(int h, int l1, int l2){
    int offset = h * N;
    bool f = unite(offset + l1, offset + l2);
    if(!f || !h) return;
    int shift = 1 << (h - 1);
    push(h - 1, l1, l2);
    push(h - 1, l1 + shift, l2 + shift);
  }
public:
  parallel_union_find(int N): N(N), H(0){
    int M = 1;
    while(M <= N) M <<= 1, H++;
    sz.resize(N * H, 1);
    par.resize(N * H);
    std::iota(par.begin(), par.end(), 0);
  }
  // unite(l1 + i, l2 + i), 0 <= i < k
  void parallel_unite(int l1, int l2, int k){
    assert(0 <= k && l1 + k <= N && l2 + k <= N);
    for(int h = H - 1, offset = 0; h >= 0; h--){
      int w = 1 << h;
      if(k >= w){
        push(h, l1 + offset, l2 + offset);
        k -= w;
        offset += w;
      }
    }
  }
  int find(int a){
    if(par[a] == a) return a;
    return par[a] = find(par[a]);
  }
  bool unite(int a, int b){
    a = find(a);
    b = find(b);
    if(a == b) return false;
    if(sz[a] > sz[b]) std::swap(a, b);
    par[a] = b;
    sz[b] += sz[a];
    return true;
  }
  int size(int a){
    assert(a < N);
    return sz[find(a)];
  }
  bool same(int a, int b){
    assert(a < N && b < N);
    return find(a) == find(b);
  }
};
#endif