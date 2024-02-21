#ifndef _DYNAMIC_ARRAY_H_
#define _DYNAMIC_ARRAY_H_
#include <vector>
#include <iostream>
#include <cassert>

template<typename Val>
struct dynamic_array{
private:
  struct node{
    node *l, *r;
    int sz;
    bool flip;
    Val val;
    node(Val _val): l(nullptr), r(nullptr), sz(1), flip(false), val(_val){}
  };
  node *root;
  int size(node *v){
    return !v ? 0 : v->sz;
  }
  inline void update(node *v){
    v->sz = 1 + (v->l ? v->l->sz : 0) + (v->r ? v->r->sz : 0);
  }
  void push_down(node *v){
    if(!v) return;
    if(v->flip){
      if(v->l) flip(v->l);
      if(v->r) flip(v->r);
      v->flip = false;
    }
  }
  inline void flip(node *v){
    std::swap(v->l, v->r);
    v->flip ^= 1;
  }
  // vの左の子をvの位置に持ってくる
  node *rotate_right(node *v){
    node *l = v->l;
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  // vの右の子をvの位置に持ってくる
  node *rotate_left(node *v){
    node *r = v->r;
    v->r = r->l;
    r->l = v;
    update(v);
    update(r);
    return r;
  }
  // zig-zig, zig-zag
  node *splay_top_down(node *v, int k, node* &u){
    push_down(v);
    int szl = v->l ? v->l->sz : 0;
    if(k == szl) return u = v;
    if(k < szl){
      v->l = splay_top_down(v->l, k, u);
      update(v);
      if(v->l == u) return v;
      if(v->l->l == u) v = rotate_right(v);
      else v->l = rotate_left(v->l);
      return rotate_right(v);
    }else{
      v->r = splay_top_down(v->r, k - szl - 1, u);
      update(v);
      if(v->r == u) return v;
      if(v->r->r == u) v = rotate_left(v);
      else v->r = rotate_right(v->r);
      return rotate_left(v);
    }
    return v;
  }
  // zig
  node *splay_top_down(node *v, int k){
    node *u = nullptr;
    v = splay_top_down(v, k, u);
    if(v->l == u) return rotate_right(v);
    else if(v->r == u) return rotate_left(v);
    return v;
  }
  node *merge(node *l, node *r){
    if(!l || !r) return !l ? r : l;
    r = splay_top_down(r, 0);
    r->l = l;
    update(r);
    return r;
  }
  // 左がkノードになるように分割
  std::pair<node*, node*> split(node *v, int k){
    int n = size(v);
    if(k >= n) return {v, nullptr};
    v = splay_top_down(v, k);
    node *l = v->l;
    v->l = nullptr;
    update(v);
    return {l, v};
  }
  std::tuple<node*, node*, node*> split3(node *v, int l, int r){
    if(l == 0){
      auto [b, c] = split(v, r);
      return {nullptr, b, c};
    }
    v = splay_top_down(v, l - 1); //    (l - 1個)  /  v  / (残り)
    auto [b, c] = split(v->r, r - l); // cがnullptrまたはcの左が空
    v->r = nullptr; // vの右が空
    update(v);
    return {v, b, c};
  }
  // split3によって分割された組でないと壊れる
  node *merge3(node *a, node *b, node *c){
    node *v = merge(b, c); // O(1)
    if(!a) return v;
    a->r = v; // O(1)
    update(a);
    return a;
  }
  node *set_inner(node *v, int k, Val x){
    v = splay_top_down(v, k);
    v->val = x;
    update(v);
    return v;
  }
  node *get_inner(node *v, int k, Val &x){
    v = splay_top_down(v, k);
    x = v->val;
    return v;
  }
  node *flip_inner(node *v, int l, int r){
    if(r == l) return v;
    auto [a, b, c] = split3(v, l, r);
    if(b) flip(b);
    return merge3(a, b, c);
  }
  node *insert_inner(node *v, int k, node *u){
    if(k == size(v)){
      u->l = v;
      update(u);
      return u;
    }
    if(k == 0){
      u->r = v;
      update(u);
      return u;
    }
    v = splay_top_down(v, k);
    u->l = v->l;
    v->l = u;
    update(u);
    update(v);
    return v;
  }
  node *erase_inner(node *v, int k){
    v = splay_top_down(v, k);
    return merge(v->l, v->r);
  }
  node *build(const std::vector<Val> &v, int l, int r){
    int m = (l + r) >> 1;
    node *u = new node(v[m]);
    if(m > l) u->l = build(v, l, m);
    if(r > m + 1) u->r = build(v, m + 1, r);
    update(u);
    return u;
  }
  dynamic_array(node *r): root(r){}
public:
  dynamic_array(): root(nullptr){}
  dynamic_array(const std::vector<Val> &v): root(nullptr){
    if(!v.empty()) root = build(v, 0, v.size());
  }
  int size(){
    return size(root);
  }
  void insert(int k, Val x){
    assert(0 <= k && k <= size());
    root = insert_inner(root, k, new node(x));
  }
  void erase(int k){
    assert(0 <= k && k < size());
    root = erase_inner(root, k);
  }
  void set(int k, Val x){
    assert(0 <= k && k < size());
    root = set_inner(root, k, x);
  }
  Val get(int k){
    assert(0 <= k && k < size());
    Val res;
    root = get_inner(root, k, res);
    return res;
  }
  void flip(int l, int r){
    assert(0 <= l && r <= size());
    root = flip_inner(root, l, r);
  }
  // [l, r)をk左巡回シフト(k = 1のとき, a[l]がr-1に移動する)
  void cyclic_lshift(int l, int r, int k){
    if(l == r) return;
    k %= r - l;
    if(!k) return;
    auto [a, b, c] = split3(root, l, r);
    auto [bl, br] = split(b, k);
    root = merge3(a, merge(br, bl), c);
  }
  // [l, r)をk左巡回シフト(k = 1のとき, a[l]がl+1に移動する)
  void cyclic_rshift(int l, int r, int k){
    if(l == r) return;
    int len = r - l;
    k %= len;
    if(!k) return;
    cyclic_lshift(l, r, len - k);
  }
};
#endif