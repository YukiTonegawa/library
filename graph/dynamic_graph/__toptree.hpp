#ifndef _TOPTREE_H_
#define _TOPTREE_H_
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>

// パス, 部分木に対して更新, 求値ができる
template<typename top_tree_structure, typename Val, typename Lazy>
struct toptree{
  struct node{
    node *l, *r, *p, *light_top, *l2, *r2, *p2;
    int subtree_size, heavy_size, light_size;
    bool flip;
    Val val, subtree_sum, heavy_sum, light_sum;
    Lazy lazy_subtree, lazy_heavy;
    node(Val _val = top_tree_structure::template id<Val>()):
    l(nullptr), r(nullptr), p(nullptr), light_top(nullptr), l2(nullptr), r2(nullptr), p2(nullptr),
    subtree_size(1), heavy_size(1), light_size(1), flip(false),
    val(_val), subtree_sum(_val), heavy_sum(_val), light_sum(_val),
    lazy_subtree(top_tree_structure::template id2<Lazy>()), lazy_heavy(lazy_subtree){}
    bool is_root_heavy(){return !p || (p->l != this && p->r != this);}
    bool is_root_light(){return !p2 || (p2->l2 != this && p2->r2 != this);}
  };
  static node *make_node(Val val){return new node(val);}
  toptree(){}
private:
  static void update_heavy(node *v){
    v->subtree_size = v->heavy_size = 1;
    v->subtree_sum = v->heavy_sum = v->val;
    Val left_subtree = top_tree_structure::template id1<Val>();
    Val right_subtree = left_subtree;
    Val light = left_subtree;
    Val left_heavy = left_subtree;
    Val right_heavy = left_subtree;
    if(v->l){
      v->subtree_size += v->l->subtree_size;
      v->heavy_size += v->l->heavy_size;
      left_subtree = v->l->subtree_sum;
    }
    if(v->r){
      v->subtree_size += v->r->subtree_size;
      v->heavy_size += v->r->heavy_size;
      right_subtree = v->r->subtree_sum;
    }
    if(v->light_top){
      v->subtree_size += v->light_top->light_size;
      light = v->light_top->light_sum;
    }
    top_tree_structure::template merge_subtree<Val>(v->subtree_sum, left_subtree, right_subtree, light);
    top_tree_structure::template merge_subtree<Val>(v->subtree_sum, left_heavy, right_heavy, top_tree_structure::template id1<Val>());
  }
  static void update_light(node *v){
    v->light_sum = v->subtree_sum;
    v->light_size = v->subtree_size;
    if(v->l2){
      v->light_size += v->l2->light_size;
      top_tree_structure::template merge_light<Val>(v->light_sum, v->l2->light_sum);
    }
    if(v->r2){
      v->light_size += v->r2->light_size;
      top_tree_structure::template merge_light<Val>(v->light_sum, v->r2->light_sum);
    }
  }
  static void push_down(node *v){
    if(v->lazy_subtree != top_tree_structure::template id2<Lazy>()){
      if(v->l) propagate_subtree(v->l, v->lazy_subtree);
      if(v->r) propagate_subtree(v->r, v->lazy_subtree);
      if(v->light_top) propagate_subtree(v->light_top, v->lazy_subtree);
      if(v->l2) propagate_subtree(v->l2, v->lazy_subtree);
      if(v->r2) propagate_subtree(v->r2, v->lazy_subtree);
      v->lazy_subtree = top_tree_structure::template id2<Lazy>();
    }
    if(v->lazy_heavy != top_tree_structure::template id2<Lazy>()){
      if(v->l) propagate_heavy(v->l, v->lazy_heavy);
      if(v->r) propagate_heavy(v->r, v->lazy_heavy);
      v->lazy_heavy = top_tree_structure::template id2<Lazy>();
    }
    if(v->flip){
      if(v->l) flip(v->l);
      if(v->r) flip(v->r);
      v->flip = false;
    }
  }
  static void propagate_subtree(node *v, Lazy x){
    top_tree_structure::template propagate<Lazy>(v->lazy_subtree, x);
    top_tree_structure::template apply<Val, Lazy>(v->val, x, 1);
    top_tree_structure::template apply<Val, Lazy>(v->subtree_sum, x, v->subtree_size);
    top_tree_structure::template apply<Val, Lazy>(v->heavy_sum, x, v->heavy_size);
    top_tree_structure::template apply<Val, Lazy>(v->light_sum, x, v->light_size);
  }
  static void propagate_heavy(node *v, Lazy x){
    top_tree_structure::template propagate<Lazy>(v->lazy_heavy, x);
    top_tree_structure::template apply<Val, Lazy>(v->val, x, 1);
    top_tree_structure::template apply<Val, Lazy>(v->subtree_sum, x, v->heavy_size);
    top_tree_structure::template apply<Val, Lazy>(v->heavy_sum, x, v->heavy_size);
    top_tree_structure::template apply<Val, Lazy>(v->light_sum, x, v->heavy_size);
  }
  static void flip(node *v){
    std::swap(v->l, v->r);
    top_tree_structure::template flip<Val>(v->subtree_sum);
    v->flip ^= 1;
  }
  static void rotate_right_heavy(node *v){
    node *p = v->p, *pp = p->p;
    if((p->l = v->r)) v->r->p = p;
    v->r = p, p->p = v;
    update_heavy(p), update_heavy(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      update_heavy(pp);
    }
  }
  static void rotate_left_heavy(node *v){
    node *p = v->p, *pp = p->p;
    if((p->r = v->l)) v->l->p = p;
    v->l = p, p->p = v;
    update_heavy(p), update_heavy(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      update_heavy(pp);
    }
  }
  static void splay_heavy(node *v){
    node *u = nullptr;
    push_down(v);
    while(!v->is_root_heavy()){
      node *p = v->p;
      if(p->is_root_heavy()){
        u = p;
        splay_light(u);
        push_down(p), push_down(v);
        if(p->l == v) rotate_right_heavy(v);
        else rotate_left_heavy(v);
      }else{
        node *pp = p->p;
        if(pp->is_root_heavy()) u = pp, splay_light(u);
        push_down(pp), push_down(p), push_down(v);
        if(pp->l == p){
          if(p->l == v) rotate_right_heavy(p);
          else rotate_left_heavy(v);
          rotate_right_heavy(v);
        }else{
          if(p->r == v) rotate_left_heavy(p);
          else rotate_right_heavy(v);
          rotate_left_heavy(v);
        }
      }
    }
    if(u){
      node *l2 = u->l2, *r2 = u->r2, *p2 = u->p2;
      u->l2 = u->r2 = u->p2 = nullptr;
      if((v->l2 = l2)) l2->p2 = v;
      if((v->r2 = r2)) r2->p2 = v;
      update_light(v);
      if((v->p2 = p2)){
        p2->light_top = v;
        update_heavy(p2);
      }
    }
  }
  static void rotate_right_light(node *v){
    node *p = v->p2, *pp = p->p2;
    if((p->l2 = v->r2)) v->r2->p2 = p;
    v->r2 = p, p->p2 = v;
    update_light(p), update_light(v);
    if((v->p2 = pp)){
      if(pp->l2 == p) pp->l2 = v;
      if(pp->r2 == p) pp->r2 = v;
      update_light(pp);
    }
  }
  static void rotate_left_light(node *v){
    node *p = v->p2, *pp = p->p2;
    if((p->r2 = v->l2)) v->l2->p2 = p;
    v->l2 = p, p->p2 = v;
    update_light(p), update_light(v);
    if((v->p2 = pp)){
      if(pp->l2 == p) pp->l2 = v;
      if(pp->r2 == p) pp->r2 = v;
      update_light(pp);
    }
  }
  static void splay_light(node *v){
    push_down(v);
    while(!v->is_root_light()){
      node *p = v->p2;
      if(p->is_root_light()){
        push_down(p), push_down(v);
        if(p->l2 == v) rotate_right_light(v);
        else rotate_left_light(v);
      }else{
        node *pp = p->p2;
        push_down(pp), push_down(p), push_down(v);
        if(pp->l2 == p){
          if(p->l2 == v) rotate_right_light(p);
          else rotate_left_light(v);
          rotate_right_light(v);
        }else{
          if(p->r2 == v) rotate_left_light(p);
          else rotate_right_light(v);
          rotate_left_light(v);
        }
      }
    }
    if(v->p2){
      v->p2->light_top = v;
      update_heavy(v->p2);
    }
  }
  static void insert_light(node *p, node *c){
    push_down(p);
    push_down(c);
    if((c->l2 = p->light_top)) p->light_top->p2 = c;
    update_light(c);
    p->light_top = c;
    c->p2 = p;
    update_heavy(p);
  }
  static void erase_light(node *p, node *c){
    splay_light(c);
    node *l = c->l2, *r = c->r2;
    c->l2 = c->r2 = c->p2 = nullptr;
    if(l && r){
      l->p2 = r->p2 = nullptr;
      while(l->r2) l = l->r2;
      splay_light(l);
      l->r2 = r;
      r->p2 = l;
      update_light(l);
    }else if(r) l = r;
    if(l) l->p2 = p;
    p->light_top = l;
    update_heavy(p);
  }
  static void swap_light(node *p, node *a, node *b){
    push_down(a);
    splay_light(b);
    node *l = b->l2, *r = b->r2;
    b->l2 = b->r2 = b->p2 = nullptr;
    if((a->l2 = l)) l->p2 = a;
    if((a->r2 = r)) r->p2 = a;
    if((a->p2 = p)) p->light_top = a;
    update_light(a);
  }
  static node *expose(node *v){
    node *c = nullptr;
    for(node *u = v; u; u = u->p){
      splay_heavy(u);
      if(c && u->r) swap_light(u, u->r, c);
      else if(c) erase_light(u, c);
      else if(u->r) insert_light(u, u->r);
      u->r = c;
      update_heavy(u);
      c = u;
    }
    splay_heavy(v);
    return c;
  }
public:
  // vを根にする
  static node *evert(node *v){
    expose(v);
    flip(v);
    push_down(v);
    return v;
  }
  // 非連結であることが保証されている場合
  static void _link(node *p, node *c){
    evert(c);
    expose(p);
    c->p = p;
    p->r = c;
    update_heavy(p);
  }
  // 辺(a, b)があることが保証されている場合
  static bool _cut(node *u, node *v){
    evert(u);
    return cut_from_parent(v);
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
  // cの親との辺を切る, false: 辺を切れなかった
  static bool cut_from_parent(node *c){
    expose(c);
    node *p = c->l;
    if(p == nullptr) return false;
    c->l = p->p = nullptr;
    update_heavy(c);
    return true;
  }
  static node *get_root(node *v){
    expose(v);
    while(v->l){
      push_down(v);
      v = v->l;
    }
    splay_heavy(v);
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
    return v->heavy_size - 1;
  }
  static node *lca(node *u, node *v){
    if(!is_same(u, v)) return nullptr;
    return _lca(u, v);
  }
  // depth(v) < kの場合根を返す
  static node *la(node *v, int k){
    expose(v);
    k = std::max(0, v->heavy_size - 1 - k);
    while(true){
      push_down(v);
      int lsz = v->l ? v->l->heavy_size : 0;
      if(lsz == k) break;
      if(lsz > k) v = v->l;
      else v = v->r, k -= lsz + 1;
    }
    splay_heavy(v);
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
  static int size_subtree(node *v){
    expose(v);
    if(!v->l) return v->subtree_size;
    else{
      node *p = v->l;
      v->l = p->p = nullptr;
      update_heavy(v);
      int ret = v->subtree_size;
      v->l = p;
      p->p = v;
      update_heavy(v);
      return ret;
    }
  }
  static void set(node *v, Val x){
    expose(v);
    v->val = x;
    update_heavy(v);
  }
  static Val get(node *v){
    expose(v);
    return v->val;
  }
  // root...vを更新
  static void update_path(node *v, Lazy x){
    expose(v);
    propagate_heavy(v, x);
  }
  static Val query_path(node *v){
    expose(v);
    return v->heavy_sum;
  }
  static void update_subtree(node *v, Lazy x){
    expose(v);
    if(!v->l) propagate_subtree(v, x);
    else{
      node *p = v->l;
      v->l = p->p = nullptr;
      update_heavy(v);
      propagate_subtree(v, x);
      push_down(v);
      v->l = p;
      p->p = v;
      update_heavy(v);
    }
  }
  static Val query_subtree(node *v){
    expose(v);
    if(!v->l) return v->subtree_sum;
    else{
      node *p = v->l;
      v->l = p->p = nullptr;
      update_heavy(v);
      Val ret = v->subtree_sum;
      v->l = p;
      p->p = v;
      update_heavy(v);
      return ret;
    }
  }
};

namespace top_tree_structures{
  /*
    パス更新, 部分木クエリを両立させる場合,
    v->subtree_sum <- apply(v->subtree_sum, x, v->heavy_size)が成り立つ必要がある
  */
  /*
  struct name{
    // 要素の単位元
    template<typename T>
    static constexpr T id1(){

    }
    // 遅延させるものの単位元
    template<typename E>
    static constexpr E id2(){

    }
    // heavy_edgeの子, 親を反転
    template<typename T>
    static constexpr void flip(T &a){

    }
    // vにheavy_edgeの親(左), 子(右), light_edgeをマージしたものをマージ
    // light_edgeを無視するとパスクエリ
    template<typename T>
    static constexpr void merge_subtree(T &v, const T &p, const T &c, const T &l){
      //e.g. v += p + c + l
    }
    // light_edgeをマージ
    template<typename T>
    static constexpr void merge_light(T &v, const T &c){

    }
    // 作用
    template<typename T, typename E>
    static constexpr void apply(T &v, const E &lz, int sz){

    }
    // 遅延させているものをマージ
    template<typename E>
    static constexpr E propagate(E &a, const E &b){

    }
  };
  */
  struct subtree_add_subtree_sum{
    // 要素の単位元
    template<typename T>
    static constexpr T id1(){
      return 0;
    }
    // 遅延させるものの単位元
    template<typename E>
    static constexpr E id2(){
      return 0;
    }
    // heavy_edgeの子, 親を反転
    template<typename T>
    static constexpr void flip(T &a){
      return;
    }
    // vにheavy_edgeの親(左), 子(右), light_edgeをマージしたものをマージ
    // light_edgeを無視するとパスクエリ
    template<typename T>
    static constexpr void merge_subtree(T &v, const T &p, const T &c, const T &l){
      //e.g. v += p + c + l
      v += p + c + l;
    }
    // light_edgeをマージ
    template<typename T>
    static constexpr void merge_light(T &v, const T &c){
      v += c;
    }
    // 作用
    template<typename T, typename E>
    static constexpr void apply(T &v, const E &lz, int sz){
      v += lz * sz;
    }
    // 遅延させているものをマージ
    template<typename E>
    static constexpr void propagate(E &a, const E &b){
      a += b;
    }
  };


  struct max_dist_struct{
    long long key, allsz, plen, clen;
    max_dist_struct(long long key = 0): key(key), allsz(key), plen(key), clen(key){}
  };
  struct max_dist{
    template<typename T>
    static constexpr T id1(){
      return T(0);
    }
    // no use
    template<typename E>
    static constexpr E id2(){
      return E(0);
    }
    // heavy_edgeの子, 親を反転
    template<typename T>
    static constexpr void flip(T &a){
      std::swap(a.plen, a.clen);
    }
    template<typename T>
    static constexpr void merge_subtree(T &v, const T &p, const T &c, const T &l){
      v.allsz = p.allsz + v.key + c.allsz;
      v.plen = std::max(c.plen, c.allsz + v.key + std::max(l.clen, p.plen));
      v.clen = std::max(p.clen, p.allsz + v.key + std::max(l.clen, c.clen));
    }
    // 高さが増えない(light_edge同士をマージ)
    template<typename T>
    static constexpr void merge_light(T &v, const T &c){
      v.clen = std::max(v.clen, c.clen);
    }
    // no use
    template<typename T, typename E>
    static constexpr void apply(T &a, const E &b, int sz){

    }
    // no use
    template<typename E>
    static constexpr void propagate(E &a, const E &b){

    }
  };
};
#endif