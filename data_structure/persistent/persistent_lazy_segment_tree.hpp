#ifndef _PERSISTENT_LAZY_SEGMENT_TREE_H_
#define _PERSISTENT_LAZY_SEGMENT_TREE_H_

#include <vector>
#include <cassert>
#include <iostream>
#include "../../algebraic_structure/monoid.hpp"

template<typename monoid>
struct persistent_lazy_segment_tree_iter;

template<typename monoid>
struct persistent_lazy_segment_tree{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
private:
  using node = persistent_lazy_segment_tree<monoid>;
  node *l, *r;
  Val sum;
  Lazy lazy;
  persistent_lazy_segment_tree(){}
  persistent_lazy_segment_tree(Val val): l(nullptr), r(nullptr), sum(val), lazy(id_lazy()){}
  static node *make_node(Val x = id()){
    return new node(x);
  }
  static node *copy_node(node *v){
    if(!v) return new node(id());
    return new node(*v);
  }
  static void propagate(node *v, Lazy x, int l, int r){
    v->sum = apply(v->sum, x, l, r);
    v->lazy = propagate_lazy(v->lazy, x);
  }
  static void push_down(node *v, int l, int r){
    if(!v || v->lazy == id_lazy()) return;
    int mid = (l + r) >> 1;
    v->l = copy_node(v->l);
    v->r = copy_node(v->r);
    propagate(v->l, v->lazy, l, mid);
    propagate(v->r, v->lazy, mid, r);
    v->lazy = id_lazy();
  }
  static node *eval(node *v){
    v->sum = merge(v->l ? v->l->sum : id(), v->r ? v->r->sum : id());
    return v;
  }
  static node *set(node *v, int k, Val x, int l, int r){
    if(r - l == 1){
      v = copy_node(v);
      v->sum = x;
      return v;
    }
    push_down(v, l, r);
    v = copy_node(v);
    int mid = ((long long)l + r) >> 1;
    if(mid <= k) v->r = set(v->r, k, x, mid, r);
    else v->l = set(v->l, k, x, l, mid);
    return eval(v);
  }
  static node *update(node *v, int a, int b, Lazy x, int l, int r){
    if(r <= a || b <= l) return v;
    if(a <= l && r <= b){
      v = copy_node(v);
      propagate(v, x, l, r);
      return v;
    }
    push_down(v, l, r);
    v = copy_node(v);
    int mid = ((long long)l + r) >> 1;
    v->r = update(v->r, a, b, x, mid, r);
    v->l = update(v->l, a, b, x, l, mid);
    return eval(v);
  }
  static Val get(node *v, int k, int l, int r){
    if(!v) return id();
    if(r - l == 1) return v->sum;
    push_down(v, l, r);
    int mid = ((long long)l + r) >> 1;
    if(mid <= k) return get(v->r, k, mid, r);
    else return get(v->l, k, l, mid);
  }
  static Val query(node *v, int a, int b, int l, int r){
    if(!v || b <= l || r <= a) return id();
    if(a <= l && r <= b) return v->sum;
    push_down(v, l, r);
    int mid = ((long long)l + r) >> 1;
    return merge(query(v->l, a, b, l, mid), query(v->r, a, b, mid, r));
  }
  template<typename F>
  static std::pair<int, Val> bisect_from_left(node *v, const int l, int a, int b, const F &f, Val ok){
    if(b <= l) return {-1, ok};
    push_down(v, a, b);
    if(l <= a){
      Val m = merge(ok, v->sum);
      if(!f(m)) return {-1, m};
      if(b - a == 1) return {a, m};
    }
    std::pair<int, Val> x{-1, id()};
    int mid = (a + b) >> 1;
    if(v->l) x = bisect_from_left(v->l, l, a, mid, f, ok);
    if(x.first != -1) return x;
    if(v->r) x = bisect_from_left(v->r, l, mid, b, f, ok);
    return x;
  }
  template<typename F>
  static std::pair<int, Val> bisect_from_right(node *v, const int r, int a, int b, const F &f, Val ok){
    if(r < a) return {-1, ok};
    push_down(v, a, b);
    if(b <= r + 1){
      Val m = merge(v->sum, ok);
      if(!f(m)) return {-1, m};
      if(b - a == 1) return {a, m};
    }
    std::pair<int, Val> x{-1, id()};
    int mid = (a + b) >> 1;
    if(v->r) x = bisect_from_right(v->r, r, mid, b, f, ok);
    if(x.first != -1) return x;
    if(v->l) x = bisect_from_right(v->l, r, a, mid, f, ok);
    return x;
  }

  friend persistent_lazy_segment_tree_iter<monoid>;
};

template<typename monoid>
struct persistent_lazy_segment_tree_iter{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
private:
  using node = persistent_lazy_segment_tree<monoid>;
  using iter = persistent_lazy_segment_tree_iter<monoid>;
  int lx, rx;
  node *root;
  persistent_lazy_segment_tree_iter(int minx, int maxx, node *v): lx(minx), rx(maxx), root(v){}
public:
  persistent_lazy_segment_tree_iter(int minx, int maxx): lx(minx), rx(maxx){
    assert(lx < rx);
    root = node::make_node();
  }
  iter set(int k, Val x){
    assert(lx <= k && k < rx);
    return iter(lx, rx, node::set(root, k, x, lx, rx));
  }
  iter update(int l, int r, Lazy x){
    assert(lx <= l && r <= rx);
    return iter(lx, rx, node::update(root, l, r, x, lx, rx));
  }
  Val get(int k){
    assert(lx <= k && k < rx);
    return node::get(root, k, lx, rx);
  }
  Val query(int l, int r){
    assert(lx <= l && r <= rx);
    return node::query(root, l, r, lx, rx);
  }
  // f(sum[l, r])がtrueになる最左のr. ない場合は-1
  template<typename F>
  int bisect_from_left(int l, const F &f){
    assert(root && lx <= l && l < rx);
    assert(!f(id()));
    return node::bisect_from_left(root, l, lx, rx, f, id()).first;
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合は-1
  template<typename F>
  int bisect_from_right(int r, const F &f){
    assert(root && lx <= r && r < rx);
    assert(!f(id()));
    return node::bisect_from_right(root, r, lx, rx, f, id()).first;
  }
};
#endif