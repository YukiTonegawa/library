#ifndef _LINK_CUT_TREE_PLAIN_H_
#define _LINK_CUT_TREE_PLAIN_H_

#include <vector>
#include <algorithm>
#include <cassert>

template<typename Label = int>
struct link_cut_tree_plain{
  struct node{
    node *l, *r, *p;
    int sz;
    Label label;
    bool flip;
    node(Label label = -1):
    l(nullptr), r(nullptr), p(nullptr), sz(1), label(label), flip(false){}
    bool is_root(){
      return !p || (p->l != this && p->r != this);
    }
  };
  static node *make_node(Label label = -1){return new node(label);}
  link_cut_tree_plain(){}
private:
  static void update(node *v){
    v->sz = 1 + (v->l ? v->l->sz : 0) + (v->r ? v->r->sz : 0);
  }
  static void push_down(node *v){
    if(v->flip){
      if(v->l) flip(v->l);
      if(v->r) flip(v->r);
      v->flip = false;
    }
  }
  static void flip(node *v){
    std::swap(v->l, v->r);
    v->flip ^= 1;
  }
  static void rotate_right(node *v){
    node *p = v->p, *pp = p->p;
    if((p->l = v->r)) v->r->p = p;
    v->r = p, p->p = v;
    update(p), update(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      update(pp);
    }
  }
  static void rotate_left(node *v){
    node *p = v->p, *pp = p->p;
    if((p->r = v->l)) v->l->p = p;
    v->l = p, p->p = v;
    update(p), update(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      update(pp);
    }
  }
  static void splay(node *v){
    push_down(v);
    while(!v->is_root()){
      node *p = v->p;
      if(p->is_root()){
        push_down(p), push_down(v);
        if(p->l == v) rotate_right(v);
        else rotate_left(v);
      }else{
        node *pp = p->p;
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
  // 連結なことが保証されている場合
  static node *_lca(node *u, node *v){
    expose(u);
    return expose(v);
  }
  // 連結であることが保証されている場合
  static int _dist(node *u, node *v){
    return depth(u) + depth(v) - 2 * depth(_lca(u, v));
  }
  // 非連結であることが保証されている場合
  static void _link(node *p, node *c){
    evert(c);
    expose(p);
    c->p = p;
    p->r = c;
    update(p);
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
      u->r = c;
      update(u);
      c = u;
    }
    splay(v);
    return c;
  }
  static node *get_root(node *v){
    expose(v);
    while(v->l){
      push_down(v);
      v = v->l;
    }
    splay(v);
    return v;
  }
  // 同じ連結成分か
  static bool is_same(node *u, node *v){
    if(!u || !v) return false;
    return get_root(u) == get_root(v);
  }
  // 0-indexed
  static int depth(node *v){
    expose(v);
    return v->sz - 1;
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
  // s-tパスに含まれるdist(s, t) + 1個のノードのうちk番目(0-indexed)
  // 非連結な場合, kが条件を満たさない場合 nullptr
  static node *kth_node_of_path(node *s, node *t, int k){
    if(!is_same(s, t) || k < 0) return nullptr;
    if(k == 0) return s;
    node *l = _lca(s, t);
    int ds = depth(s), dt = depth(t), dl = depth(l), dist = ds + dt - 2 * dl;
    if(k <= ds - dl) return la(s, k);
    else if(k <= dist) return la(t, dist - k);
    else return nullptr;
  }
  // 非連結な場合 -1
  static int dist(node *u, node *v){
    if(!is_same(u, v)) return -1;
    return _dist(u, v);
  }
  // false: 辺を繋げなかった
  static bool link(node *p, node *c){
    if(is_same(p, c)) return false;
    _link(p, c);
    return true;
  }
  // cの親との辺を切る
  static void cut_from_parent(node *c){
    expose(c);
    node *p = c->l;
    if(p == nullptr) return;
    c->l = p->p = nullptr;
    update(c);
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
  static bool is_there_edge(node *u, node *v){
    evert(u);
    return u == get_parent(v);
  }
};
#endif
