#ifndef _DUAL_SEGMENT_TREE_H_
#define _DUAL_SEGMENT_TREE_H_
#include <cassert>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#include "../../algebraic_structure/monoid.hpp"
#include "../../algebraic_structure/abelian_group.hpp"

// 可換群
template<typename abelian_group>
struct dual_segment_tree_abelian_group{
  using Val = typename abelian_group::Val;
  static constexpr auto id = abelian_group::id;
  static constexpr auto inv = abelian_group::inv;
  static constexpr auto merge = abelian_group::merge;
  int M = 1;
  std::vector<Val> val;
  std::vector<Val> lazy;
  dual_segment_tree_abelian_group(){}
  dual_segment_tree_abelian_group(int N): M(N), val(M, id()), lazy(M + 1 , id()){assert(N > 0);}
  dual_segment_tree_abelian_group(const std::vector<Val> &v): M(v.size()), val(v), lazy(M + 1 , id()){}
  void set(int k, Val x){
    Val y = val[k];
    Val z = get(k);
    val[k] = merge(x, inv(merge(inv(y), z)));
  }
  Val get(int k){
    Val res = val[k];
    for(int i = k + 1; i <= M; i += (i & (-i))) res = merge(res, lazy[k]);
    return res;
  }
  void update(int r, Val x){
    Val ret = 0;
    for(int k = r; k > 0; k -= (k & (-k))) lazy[k] = merge(lazy[k], x);
  }
  void update(int l, int r, Val x){
    update(r, x);
    update(l, inv(x));
  }
};
template<typename monoid>
struct dual_segment_tree_monoid{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  int N, M, H;
  std::vector<Val> val;
  std::vector<Val> lazy;
private:
  int ceil_pow2(int y){
    int x = 0;
    while ((1U << x) < (unsigned int)(y)) x++;
    return x;
  };
  void push_down(int k){
    if(k >= M - 1 || lazy[k] == id()) return;
    lazy[k * 2 + 1] = merge(lazy[k * 2 + 1], lazy[k]);
    lazy[k * 2 + 2] = merge(lazy[k * 2 + 2], lazy[k]);
    lazy[k] = id();
  }
  void set(int a, Val x, int k){
    for(int i = H - 1; i >= 0; i--){
      push_down(k);
      k = k * 2 + 1 + ((a >> i) & 1);
    }
    lazy[k] = id();
    val[a] = x;
  }
  Val get(int a, int k){
    Val res = id();
    for(int i = H - 1; i >= 0; i--){
      res = merge(lazy[k], res);
      k = k * 2 + 1 + ((a >> i) & 1);
    }
    return merge(val[a], merge(lazy[k], res));
  }
  void update(int a, int b, Val x, int k, int l, int r){
    if(r <= a || b <= l) return;
    if(a <= l && r <= b){
      lazy[k] = merge(lazy[k], x);
      return;
    }
    push_down(k);
    update(a, b, x, k * 2 + 1, l, (l + r) / 2);
    update(a, b, x, k * 2 + 2, (l + r) / 2, r);
  }
public:
  dual_segment_tree_monoid(){}
  dual_segment_tree_monoid(int n, Val v = id()): N(n), M(1 << ceil_pow2(N)), H(31 -__builtin_clz(M)), val(n, v), lazy(2 * M - 1, id()){}
  dual_segment_tree_monoid(const std::vector<Val> &v): N(v.size()), M(1 << ceil_pow2(N)), H(31 -__builtin_clz(M)), val(v), lazy(2 * M - 1, id()){}
  void set(int k, Val x){set(k, x, 0);}
  Val get(int k){return get(k, 0);}
  void update(int a, int b, Val x){update(a, b, x, 0, 0, M);}
  void update_all(Val x){update(0, M, x, 0, 0, M);}
  std::vector<Val> to_list(){
    for(int i = 0; i < M - 1; i++) push_down(i);
    std::vector<Val> res = val;
    for(int i = 0; i < N; i++) res[i] = merge(res[i], lazy[M - 1 + i], i, i + 1);
    return res;
  }
};
#endif