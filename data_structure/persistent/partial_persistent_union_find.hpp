#ifndef _PARTIAL_PERSISTENT_UNION_FIND_H_
#define _PARTIAL_PERSISTENT_UNION_FIND_H_
#include <vector>
#include <tuple>
struct partial_persistent_union_find{
  int N, T;
  std::vector<std::tuple<int, int, int>> data; // {親, 辺が張られた時刻, その時点での親側のサイズ}
  static constexpr int inf_time = 1 << 30;
  partial_persistent_union_find(int n): N(n), T(0), data(n + 1, {n, inf_time, 0}){}
  // 時刻tでの根, サイズ
  std::pair<int, int> find(int v, int t){
    int s = 1;
    while(std::get<1>(data[v]) <= t){
      auto [x, y, z] = data[v];
      v = x, s = z;
    }
    return {v, s};
  }
  void unite(int a, int b){
    T++;
    auto [A, Asz] = find(a, T);
    auto [B, Bsz] = find(b, T);
    if(A == B) return;
    if(Asz > Bsz) std::swap(A, B), std::swap(Asz, Bsz);
    data[A] = {B, T, Asz + Bsz};
  }
  int size(int a, int t){
    return find(a, t).second;
  }
  bool same(int a, int b, int t){
    return find(a, t).first == find(b, t).first;
  }
  // a, bが連結になった時刻(a = bの場合0, 非連結な場合-1)
  int connect_time(int a, int b){
    int t = 0;
    while(a != b){
      auto [A, At, Asz] = data[a];
      auto [B, Bt, Bsz] = data[b];
      if(a != N && At <= Bt) a = A, t = At;
      else b = B, t = Bt;
    }
    return t == inf_time ? -1 : t;
  }
};
#endif