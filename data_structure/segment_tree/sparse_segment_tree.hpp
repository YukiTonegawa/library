#ifndef _SPARSE_SEGMENT_TREE_H_
#define _SPARSE_SEGMENT_TREE_H_
#include <vector>
#include <cassert>
#include <limits>
#include "../../algebraic_structure/monoid.hpp"

template<typename monoid>
struct sparse_segment_tree{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
private:
  using I = int;
  struct node{
    I idx;
    node *l, *r;
    Val val, sum;
    node(I idx, Val val): idx(idx), l(nullptr), r(nullptr), val(val), sum(val){}
  };
  node *make_node(I idx,  Val val = id()){
    return new node(idx, val);
  }
  void eval(node *v){
    v->sum = v->val;
    if(v->l) v->sum = merge(v->l->sum, v->sum);
    if(v->r) v->sum = merge(v->sum, v->r->sum);
  }
  void set(node* &v, I k, Val x, I l, I r){
    if(!v){
      v = make_node(k, x);
      return;
    }
    if(v->idx == k){
      v->val = x;
      eval(v);
      return;
    }
    I mid = ((long long)l + r) >> 1;
    if(k < mid){
      if(v->idx < k) std::swap(v->idx, k), std::swap(v->val, x);
      set(v->l, k, x, l, mid);
    }else{
      if(v->idx > k) std::swap(v->idx, k), std::swap(v->val, x);
      set(v->r, k, x, mid, r);
    }
    eval(v);
  }
  Val get(node *v, I k, I l, I r){
    if(!v) return id();
    if(v->idx == k) return v->val;
    I mid = ((long long)l + r) >> 1;
    if(k < mid) return get(v->l, k, l, mid);
    else return get(v->r, k, mid, r);
  }
  void update(node* &v, I k, Val x, I l, I r){
    if(!v){
      v = make_node(k, x);
      return;
    }
    if(v->idx == k){
      v->val = merge(v->val, x);
      eval(v);
      return;
    }
    I mid = ((long long)l + r) >> 1;
    if(k < mid){
      if(v->idx < k) std::swap(v->idx, k), std::swap(v->val, x);
      update(v->l, k, x, l, mid);
    }else{
      if(v->idx > k) std::swap(v->idx, k), std::swap(v->val, x);
      update(v->r, k, x, mid, r);
    }
    eval(v);
  }
  Val query(node *v, I a, I b, I l, I r){
    if(!v || b <= l || r <= a) return id();
    if(a <= l && r <= b) return v->sum;
    I mid = ((long long)l + r) >> 1;
    Val ret = query(v->l, a, b, l, mid);
    if(a <= v->idx && v->idx < b) ret = merge(ret, v->val);
    return merge(ret, query(v->r, a, b, mid, r));
  }
  template<typename F>
  I bisect_from_left(node *v, const I k, F &f, I l, I r, Val &ok){
    if(!v || r <= k) return maxx;
    Val m = merge(ok, v->sum);
    if(k <= l && !f(m)){
      ok = m;
      return maxx;
    }
    I mid = ((long long)l + r) >> 1;
    I x = bisect_from_left(v->l, k, f, l, mid, ok);
    if(x != maxx) return x;
    if(k <= v->idx){
      ok = merge(ok, v->val);
      if(f(ok)) return v->idx;
    }
    return bisect_from_left(v->r, k, f, mid, r, ok);
  }
  template<typename F>
  I bisect_from_right(node *v, const I k, F &f, I l, I r, Val ok){
    if(!v || k < l) return maxx;
    Val m = merge(v->sum, ok);
    if(r <= k + 1 && !f(m)){
      ok = m;
      return maxx;
    }
    I mid = ((long long)l + r) >> 1;
    I x = bisect_from_right(v->r, k, f, mid, r, ok);
    if(x != maxx) return x;
    if(v->idx <= k){
      ok = merge(v->val, ok);
      if(f(ok)) return v->idx;
    }
    return bisect_from_right(v->l, k, f, l, mid, ok);
  }
  node *root;
  I minx, maxx;
public:
  sparse_segment_tree(): root(nullptr), minx(0), maxx(0){}
  sparse_segment_tree(I minx, I maxx): root(nullptr), minx(minx), maxx(maxx){}

  void set(I k, Val x){
    assert(minx <= k && k < maxx);
    set(root, k, x, minx, maxx);
  }
  Val get(I k){
    assert(minx <= k && k < maxx);
    return get(root, k, minx, maxx);
  }
  void update(I k, Val x){
    assert(minx <= k && k < maxx);
    update(root, k, x, minx, maxx);
  }
  Val query(I l, I r){
    assert(minx <= l && r <= maxx);
    return query(root, l, r, minx, maxx);
  }
  Val query_all(){
    return root ? root->sum : id();
  }
  // f(sum[l, r])がtrueになる最左のr. ない場合はmaxx
  template<typename F>
  I bisect_from_left(I l, const F &f){
    Val x = id();
    return bisect_from_left(root, l, f, minx, maxx, x);
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合はmaxx
  template<typename F>
  I bisect_from_right(I r, const F &f){
    Val x = id();
    return bisect_from_right(root, r, f, minx, maxx, x);
  }
};
#endif
