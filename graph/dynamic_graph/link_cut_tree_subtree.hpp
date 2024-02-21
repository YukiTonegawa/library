#ifndef _LINK_CUT_TREE_SUBTREE_H_
#define _LINK_CUT_TREE_SUBTREE_H_

#include <vector>
#include <algorithm>
#include <cassert>
#include "../../algebraic_structure/abelian_group.hpp"

template<typename abelian_group>
struct link_cut_tree_subtree{
  using Val = typename abelian_group::Val;
  using Lazy = typename abelian_group::Lazy;
  static constexpr auto id = abelian_group::id;
  static constexpr auto id_lazy = abelian_group::id_lazy;
  static constexpr auto inv = abelian_group::inv;
  static constexpr auto inv_lazy = abelian_group::inv_lazy;
  static constexpr auto merge = abelian_group::merge;
  static constexpr auto propagate_lazy = abelian_group::propagate;
  static constexpr auto apply = abelian_group::apply;
  struct node{
    node *l, *r, *p;
    int sz, szsub, label;
    bool flip;
    Val val, valsub, sum;
    Lazy lazysum, lazyhist;
    node(Val _val = id(), int label = -1):
    l(nullptr), r(nullptr), p(nullptr), sz(1), szsub(0), label(label), flip(false),
    val(_val), valsub(id()), sum(_val), lazysum(id_lazy()), lazyhist(lazysum){}
    bool is_root(){return !p || (p->l != this && p->r != this);}
  };
  static node *make_node(Val val, int label = -1){return new node(val, label);}
  link_cut_tree_subtree(){}
private:
  static void update(node *v){
    v->sz = 1 + v->szsub;
    v->sum = v->val + v->valsub;
    if(v->l){
      v->sz += v->l->sz;
      v->sum = merge(v->l->sum, v->sum);
    }
    if(v->r){
      v->sz += v->r->sz;
      v->sum = merge(v->sum, v->r->sum);
    }
  }
  static void propagate(node *v, Lazy x){
    v->lazysum = propagate_lazy(v->lazysum, x);
    v->val = apply(v->val, x, 0, 1);
    v->valsub = apply(v->valsub, x, 0, v->szsub);
    v->sum = apply(v->sum, x, 0, v->sz);
  }
  static void flip(node *v){
    std::swap(v->l, v->r);
    v->flip ^= 1;
  }
  static void push_down(node *v){
    if(v->flip){
      if(v->l) flip(v->l);
      if(v->r) flip(v->r);
      v->flip = false;
    }
    if(v->l) fetch(v->l);
    if(v->r) fetch(v->r);
  }
  static void modify(node *v){
    if(!v || !v->p) return;
    v->lazyhist = v->p->lazysum;
  }
  static void fetch(node *v){
    node *p = v->p;
    if(!p || p->lazysum == v->lazyhist) return;
    Lazy diff = propagate_lazy(p->lazysum, inv_lazy(v->lazyhist));
    propagate(v, diff);
    v->lazyhist = p->lazysum;
  }
  static void rotate_right(node *v){
    node *p = v->p, *pp = p->p;
    if((p->l = v->r)) v->r->p = p, modify(v->r);
    v->r = p, p->p = v;
    modify(p);
    update(p), update(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      modify(v);
      update(pp);
    }
  }
  static void rotate_left(node *v){
    node *p = v->p, *pp = p->p;
    if((p->r = v->l)) v->l->p = p, modify(v->l);
    v->l = p, p->p = v;
    modify(p);
    update(p), update(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      modify(v);
      update(pp);
    }
  }
  static void splay(node *v){
    push_down(v);
    while(!v->is_root()){
      node *p = v->p;
      if(p->is_root()){
        fetch(p);
        push_down(p), push_down(v);
        if(p->l == v) rotate_right(v);
        else rotate_left(v);
      }else{
        node *pp = p->p;
        fetch(pp);
        push_down(pp), push_down(p), push_down(v);
        if(pp->l == p){
          if(p->l == v) rotate_right(p);
          else rotate_left(v);
          rotate_right(v);
        }else{
          if(p->r == v) rotate_left(p);
          else rotate_right(v);
          rotate_left(v);
        }
      }
    }
  }
public:
  static node *get_root(node *v){
    expose(v);
    while(v->l){
      push_down(v);
      v = v->l;
    }
    splay(v);
    return v;
  }
  static bool is_same(node *u, node *v){
    if(!u || !v) return false;
    return get_root(u) == get_root(v);
  }
  // 非連結であることが保証されている場合
  static void _link(node *p, node *c){
    evert(c);
    expose(p);
    c->p = p;
    p->r = c;
    modify(c);
    update(p);
  }
  // cの親との辺を切る, false: 辺を切れなかった
  static void cut_from_parent(node *c){
    expose(c);
    node *p = c->l;
    if(p == nullptr) return;
    c->l = p->p = nullptr;
    update(c);
  }
  // 辺(a, b)があることが保証されている場合
  static void _cut(node *u, node *v){
    evert(u);
    cut_from_parent(v);
  }
  // vを根にする
  static node *evert(node *v){
    expose(v);
    flip(v);
    push_down(v);
    return v;
  }
  static node *expose(node *v){
    node *c = nullptr;
    for(node *u = v; u; u = u->p){
      splay(u);
      if(u->r){
        fetch(u->r);
        u->szsub += u->r->sz;
        u->valsub = merge(u->valsub, u->r->sum);
      }
      if(c){
        fetch(c);
        u->szsub -= c->sz;
        u->valsub = merge(u->valsub, inv(c->sum));
      }
      u->r = c;
      update(u);
      c = u;
    }
    splay(v);
    return c;
  }
  // val[v] = x
  static void set(node *v, Val x){
    expose(v);
    v->val = x;
    update(v);
  }
  // val[v]
  static Val get(node *v){
    expose(v);
    return v->val;
  }
  // 呼ぶ前に根をevertしておく
  static void update_subtree(node *v, Lazy x){
    node *p = get_parent(v);
    if(p){
      evert(p);
      cut_from_parent(v);
      propagate(v, x);
      _link(p, v);
    }else{
      propagate(v, x);
    }
  }
  // 呼ぶ前に根をevertしておく
  static Val query_subtree(node *v){
    expose(v);
    return merge(v->val, v->valsub);
  }
  // 部分木のサイズ
  static int size_subtree(node *v){
    expose(v);
    return 1 + v->szsub;
  }
  static int depth(node *v){
    expose(v);
    return v->sz - 1 - v->szsub;
  }
  // 連結なことが保証されている場合
  static node *_lca(node *u, node *v){
    expose(u);
    return expose(v);
  }
  // 連結であることが保証されている場合
  static int _dist(node *u, node *v){
    return depth(u) + depth(v) - 2 * depth(_lca(u, v));
  }

  static node *lca(node *u, node *v){
    if(!is_same(u, v)) return nullptr;
    return _lca(u, v);
  }
  // depth(v) < kの場合根を返す
  static node *la(node *v, int k){
    expose(v);
    k = std::max(0, v->sz - 1 - k);
    while(true){
      push_down(v);
      int lsz = v->l ? v->l->sz : 0;
      if(lsz == k) break;
      if(lsz > k) v = v->l;
      else v = v->r, k -= lsz + 1;
    }
    splay(v);
    return v;
  }
  // 非連結な場合 -1
  static int dist(node *u, node *v){
    if(!is_same(u, v)) return -1;
    return _dist(u, v);
  }
  static bool is_there_edge(node *u, node *v){
    evert(u);
    return u == get_parent(v);
  }
  static node *get_parent(node *v){
    expose(v);
    if(!v->l) return nullptr;// aが根
    push_down(v);
    v = v->l;
    while(v->r){
      push_down(v);
      v = v->r;
    }
    splay(v);
    return v;
  }
};
#endif