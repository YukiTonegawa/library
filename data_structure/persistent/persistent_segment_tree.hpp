#ifndef _PERSISTENT_SEGMENT_TREE_H_
#define _PERSISTENT_SEGMENT_TREE_H_
#include <vector>
#include <numeric>
#include <limits>
#include <cassert>
#include "../../algebraic_structure/monoid.hpp"

template<typename monoid>
struct persistent_segment_tree_iter;

template<typename monoid>
struct persistent_segment_tree{
private:
  using node = persistent_segment_tree<monoid>;
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;

  node *l, *r;
  Val sum;
  persistent_segment_tree(){}
  persistent_segment_tree(Val val): l(nullptr), r(nullptr), sum(val){}
  static node *make_node(Val x = id()){
    return new node(x);
  }
  static node *copy_node(node *v){
    if(!v) return new node(id());
    return new node(*v);
  }
  static node *eval(node *v){
    v->sum = merge(v->l ? v->l->sum : id(), v->r ? v->r->sum : id());
    return v;
  }
  static node *set(node *v, int k, Val x, int l, int r){
    v = copy_node(v);
    if(r - l == 1){
      v->sum = x;
      return v;
    }
    int mid = ((long long)l + r) >> 1;
    if(mid <= k) v->r = set(v->r, k, x, mid, r);
    else v->l = set(v->l, k, x, l, mid);
    return eval(v);
  }
  static Val get(node *v, int k, int l, int r){
    if(!v) return id();
    if(r - l == 1) return v->sum;
    int mid = ((long long)l + r) >> 1;
    if(mid <= k) return get(v->r, k, mid, r);
    else return get(v->l, k, l, mid);
  }
  static Val query(node *v, int a, int b, int l, int r){
    if(!v || b <= l || r <= a) return id();
    if(a <= l && r <= b) return v->sum;
    int mid = ((long long)l + r) >> 1;
    return merge(query(v->l, a, b, l, mid), query(v->r, a, b, mid, r));
  }
  template<typename F>
  static std::pair<int, Val> bisect_from_left(node *v, const int l, int a, int b, const F &f, Val ok){
    if(b <= l) return {-1, ok};
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
  friend persistent_segment_tree_iter<monoid>;
};

template<typename monoid>
struct persistent_segment_tree_iter{
private:
  using node = persistent_segment_tree<monoid>;
  using iter = persistent_segment_tree_iter<monoid>;
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;

  int lx, rx;
  node *root;
  persistent_segment_tree_iter(int minx, int maxx, node *v): lx(minx), rx(maxx), root(v){}
public:
  persistent_segment_tree_iter(): root(nullptr){}
  persistent_segment_tree_iter(int minx, int maxx): lx(minx), rx(maxx){
    assert(lx < rx);
    root = node::make_node();
  }
  iter set(int k, Val x){
    assert(lx <= k && k < rx);
    return iter(lx, rx, node::set(root, k, x, lx, rx));
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