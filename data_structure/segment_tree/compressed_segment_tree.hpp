#ifndef _COMPRESSED_SEGMENT_TREE_H_
#define _COMPRESSED_SEGMENT_TREE_H_
#include <vector>
#include <cassert>
#include <algorithm>
#include "segment_tree.hpp"

template<typename monoid>
struct compressed_segment_tree{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  using I = int;
  std::vector<I> Idx;
  int lb(I i){return std::lower_bound(Idx.begin(), Idx.end(), i) - Idx.begin();}
  segment_tree<monoid> seg;
  compressed_segment_tree(){}
  //　{点, 値}を与える
  compressed_segment_tree(const std::vector<std::pair<I, Val>> &v){
    int n = v.size();
    for(int i = 0; i < n; i++) Idx.push_back(v[i].first);
    std::sort(Idx.begin(), Idx.end());
    Idx.erase(std::unique(Idx.begin(), Idx.end()), Idx.end());
    int m = Idx.size();
    std::vector<Val> tmp(m, id());
    for(int i = 0; i < n; i++){
      int j = lb(v[i].first);
      tmp[j] = merge(tmp[j], v[i].second);
    }
    seg = segment_tree<monoid>(tmp);
  }
  void set(I k, Val x){
    int i = lb(k);
    assert(i < Idx.size() && Idx[i] == k);
    seg.set(i, x);
  }
  Val get(I k){
    int i = lb(k);
    assert(i < Idx.size() && Idx[i] == k);
    return seg.get(i);
  }
  Val query(I l, I r){
    return seg.query(lb(l), lb(r));
  }
  Val query_all(){
    return seg.query_all();
  }
  // f(sum[l, r])がtrueになる最左のr. ない場合は-1
  template<typename F>
  I bisect_from_left(I l, const F &f){
    int ret = seg.bisect_from_left(lb(l), f);
    if(ret == -1) return -1;
    return Idx[ret];
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合は-1
  template<typename F>
  I bisect_from_right(I r, const F &f){
    int ret = seg.bisect_from_right(lb(r), f);
    if(ret == -1) return -1;
    return Idx[ret];
  }
};

template<typename monoid>
struct compressed_lazy_segment_tree{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto propagate_lazy = monoid::propagate;
  static constexpr auto apply = monoid::apply;
  using I = int;
private:
  std::vector<I> Idx;
  int lb(I i){return std::lower_bound(Idx.begin(), Idx.end(), i) - Idx.begin();}
  int N, M;
  struct node{
    Val sum;
    Lazy lazy;
    I len;
    node(Val val = id()): sum(val), lazy(id_lazy()){}
  };
  static constexpr int ceil_pow2(int n){
    int m = 1;
    while(m < n) m <<= 1;
    return m;
  }
  std::vector<node> V;
  inline void push_down(int k, int l, int r){
    if(V[k].lazy == id_lazy()) return;
    int mid = (l + r) >> 1;
    propagate(k * 2 + 1, l, mid, V[k].lazy);
    propagate(k * 2 + 2, mid, r, V[k].lazy);
    V[k].lazy = id_lazy();
  }
  inline void propagate(int k, int l, int r, Lazy x){
    V[k].sum = apply(V[k].sum, x, Idx[l], Idx[l] + V[k].len);
    V[k].lazy = propagate_lazy(V[k].lazy, x);
  }
  Val set(int a, Val x, int k, int l, int r){
    if(r - l == 1) return V[k].sum = x;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid){
      return V[k].sum = merge(set(a, x, k * 2 + 1, l, mid), V[k * 2 + 2].sum);
    }else{
      return V[k].sum = merge(V[k * 2 + 1].sum, set(a, x, k * 2 + 2, mid, r));
    }
  }
  Val get(int a, int k, int l, int r){
    if(r - l == 1) return V[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid) return get(a, k * 2 + 1, l, mid);
    else return get(a, k * 2 + 2, mid, r);
  }
  Val update(int a, int b, Lazy x, int k, int l, int r){
    if(r <= a || b <= l) return V[k].sum;
    if(a <= l && r <= b){
      propagate(k, l, r, x);
      return V[k].sum;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return V[k].sum = merge(update(a, b, x, k * 2 + 1, l, mid), update(a, b, x, k * 2 + 2, mid, r));
  }
  Val query(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return id();
    if(a <= l && r <= b) return V[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return merge(query(a, b, k * 2 + 1, l, mid), query(a, b, k * 2 + 2, mid, r));
  }
public:
  compressed_lazy_segment_tree(): N(0), M(0){}

  // 更新, 取得される可能性のある区間[l, r)を与える
  // 点更新, 点取得は[l, l + 1)の区間更新と見なす
  compressed_lazy_segment_tree(const std::vector<std::pair<I, I>> &v){
    int n = v.size();
    for(int i = 0; i < n; i++){
      if(v[i].first >= v[i].second) continue;
      Idx.push_back(v[i].first);
      Idx.push_back(v[i].second);
    }
    std::sort(Idx.begin(), Idx.end());
    Idx.erase(std::unique(Idx.begin(), Idx.end()), Idx.end());
    N = Idx.size();
    M = ceil_pow2(N);
    V.resize(2 * M - 1, node());
    for(int i = 0; i < M; i++) V[M - 1 + i].len = (i + 1 >= N ? 0 : Idx[i + 1] - Idx[i]);
    for(int i = M - 2; i >= 0; i--) V[i].len = V[i * 2 + 1].len + V[i * 2 + 2].len;
  }
  // val[k] <- x
  Val set(I k, Val x){
    int i = lb(k);
    assert(i < Idx.size() && Idx[i] == k);
    assert(i + 1 < Idx.size() && Idx[i + 1] == k + 1);
    return set(i, x, 0, 0, M);
  }
  // val[k]
  Val get(I k){
    int i = lb(k);
    assert(i < Idx.size() && Idx[i] == k);
    assert(i + 1 < Idx.size() && Idx[i + 1] == k + 1);
    return get(i, 0, 0, M);
  }
  // sum[l, r)
  Val query(I l, I r){
    if(l >= r) return id();
    int l2 = lb(l), r2 = lb(r);
    assert(l2 < Idx.size() && Idx[l2] == l);
    assert(r2 < Idx.size() && Idx[r2] == r);
    return query(l2, r2, 0, 0, M);
  }
  Val query_all(){
    return V[0].sum;
  }
  // apply([l, r), x)
  void update(I l, I r, Lazy x){
    if(l >= r) return;
    int l2 = lb(l), r2 = lb(r);
    assert(l2 < Idx.size() && Idx[l2] == l);
    assert(r2 < Idx.size() && Idx[r2] == r);
    update(l2, r2, x, 0, 0, M);
  }
};
#endif