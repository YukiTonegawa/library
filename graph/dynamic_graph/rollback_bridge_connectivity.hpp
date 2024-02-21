#ifndef _ROLLBACK_BRIDGE_CONNECTIVITY_H_
#define _ROLLBACK_BRIDGE_CONNECTIVITY_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <tuple>

struct rollback_bridge_connectivity{
  static constexpr int node_threshold = 1 << 30; // これ以上ならノード　
  struct node{
    node *l, *r, *p;
    int label;
    int dup, mindup, mindupcnt, lazydup;// 重複度, 最小値, 最小値の数, 遅延
    long long dupsum, dupsumsub;
    int sz, szall, szsub;
    bool flip;
    node(int d = node_threshold, int id = -1): l(nullptr), r(nullptr), p(nullptr), label(id), 
    dup(d), mindup(d), mindupcnt(1), lazydup(0), dupsum(0), dupsumsub(0), sz(!bool(d)), szall(sz), szsub(0), flip(false){}
    bool is_root(){
      return !p || (p->l != this && p->r != this);
    }
  };
  node *make_node(int label = -1){return new node(node_threshold, label);}
  rollback_bridge_connectivity(){}
private:
  node *make_edge(){return new node(0);}
  inline int chmin(int &x, int y){
    if(x == y) return 1;
    else if(x < y) return 2;
    x = y;
    return 0;
  }
  void update(node *v){
    v->sz = (v->dup >= node_threshold ? 0 : 1);
    v->szall = 1 + v->szsub;
    v->dupsum = (v->dup >= node_threshold ? 0 : v->dup) + v->dupsumsub;
    v->mindup = v->dup;
    v->mindupcnt = 1;
    if(v->l){
      v->sz += v->l->sz;
      v->szall += v->l->szall;
      v->dupsum += v->l->dupsum;
      int t = chmin(v->mindup, v->l->mindup);
      if(t == 1) v->mindupcnt += v->l->mindupcnt;
      else if(t == 0) v->mindupcnt = v->l->mindupcnt;
    }
    if(v->r){
      v->sz += v->r->sz;
      v->szall += v->r->szall;
      v->dupsum += v->r->dupsum;
      int t = chmin(v->mindup, v->r->mindup);
      if(t == 1) v->mindupcnt += v->r->mindupcnt;
      else if(t == 0) v->mindupcnt = v->r->mindupcnt;
    }
  }
  // 遅延の伝播
  void push_down(node *v){
    if(v->lazydup != 0){
      if(v->l) propagate_dup(v->l, v->lazydup);
      if(v->r) propagate_dup(v->r, v->lazydup);
      v->lazydup = 0;
    }
    if(v->flip){
      if(v->l) flip(v->l);
      if(v->r) flip(v->r);
      v->flip = false;
    }
  }
  inline void propagate_dup(node *v, int x){
    v->lazydup += x;
    v->dup += x;
    v->mindup += x;
    v->dupsum += (long long)x * v->sz;
  }
  inline void flip(node *v){
    std::swap(v->l, v->r);
    v->flip ^= 1;
  }
  void rotate_right(node *v){
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
  void rotate_left(node *v){
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
  void splay(node *v){
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
  node *expose(node *v){
    node *c = nullptr;
    for(node *u = v; u; u = u->p){
      splay(u);
      if(c) heavy_link(u, c);
      if(u->r) heavy_cut(u, u->r);
      u->r = c;
      update(u);
      c = u;
    }
    splay(v);
    return c;
  }
  void heavy_link(node *p, node *c){
    p->szsub -= c->szall;
    p->dupsumsub -= c->dupsum;
  }
  void heavy_cut(node *p, node *c){
    p->szsub += c->szall;
    p->dupsumsub += c->dupsum;
  }
  node *evert(node *v){
    expose(v);
    flip(v);
    push_down(v);
    return v;
  }
  void _link(node *p, node *c){
    evert(c);
    expose(p);
    c->p = p;
    p->r = c;
    update(p);
  }
  bool _cut(node *c){
    expose(c);
    node *p = c->l;
    if(p == nullptr) return false;
    c->l = p->p = nullptr;
    update(c);
    return true;
  }
  // 辺(a, b)があることが保証されている場合
  bool _cut(node *u, node *v){
    evert(u);
    return _cut(v);
  }
  node *get_root(node *v){
    expose(v);
    while(v->l){
      push_down(v);
      v = v->l;
    }
    splay(v);
    return v;
  }
public:
  // 2頂点の関係を以下のように定義
  // 0 : 非連結, 1 : 連結かつ非2辺連結 2 : 2辺連結
  using info = std::tuple<node*, node*, node*, char, int>;
  std::vector<info> history;

  // 頂点u, vが同じ連結成分にあるか
  bool is_same(node *u, node *v){
    assert(u && v);
    return get_root(u) == get_root(v);
  }
  // 0: 非連結, 1: 連結, 2: 2辺連結
  int relation(node *u, node *v){
    assert(u && v);
    evert(u);
    if(get_root(u) != get_root(v)) return 0;
    expose(v);
    return v->mindup ? 2 : 1;
  }
  // uを含む連結成分が木か
  bool is_tree(node *u){
    assert(u);
    expose(u);
    return u->dupsum == 0;
  }
  // uを含む連結成分のサイズ
  int components_size(node *u){
    assert(u);
    expose(u);
    return u->szall + 1;
  }
  // uが2辺連結成分に含まれるか
  bool is_two_edge_connected(node *u){
    return u->dup > node_threshold ? true : false;
  }
  // link直前の関係, 橋の増減
  std::pair<int, int> link(node *u, node *v){
    assert(u && v);
    evert(u);
    if(get_root(u) != get_root(v)){
      node *e = make_edge();
      expose(v);
      // v -> e -> u
      v->r = e;
      e->p = v;
      e->r = u;
      u->p = e;
      update(e);
      update(v);
      history.push_back({u, v, e, 0, 1});
      return {0, +1};
    }else{
      expose(v);
      int m = v->mindup, b = v->mindupcnt;
      propagate_dup(v, 1);
      if(m >= 1){
        history.push_back({u, v, nullptr, 2, 0});
        return {2, 0};
      }else{
        history.push_back({u, v, nullptr, 1, -b});
        return {1, -b};
      }
    }
  }
  // 橋の増減
  int rollback(){
    assert(!history.empty());
    auto [u, v, e, t, b] = history.back();
    history.pop_back();
    if(t == 0){
      evert(u);
      expose(v);
      v->l = v->l->p = nullptr;
      update(v);
      splay(e);
      e->l = e->l->p = nullptr;
      return -b;
    }else{
      evert(u);
      expose(v);
      propagate_dup(v, -1);
      return -b;
    }
  }
};
#endif