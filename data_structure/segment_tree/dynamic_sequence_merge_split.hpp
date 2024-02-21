#ifndef _DYNAMIC_SEQUENCE_MERGE_SPLIT_H_
#define _DYNAMIC_SEQUENCE_MERGE_SPLIT_H_
#include <vector>
#include <iostream>
#include <cassert>
#include "../../algebraic_structure/monoid.hpp"

template<typename monoid>
struct dynamic_sequence_merge_split{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge_val = monoid::merge;
  static constexpr auto flip_val = monoid::flip;

private:
  struct node{
    node *l, *r;
    int sz;
    bool flip;
    Val val, sum;
    node(Val _val = id()): l(nullptr), r(nullptr), sz(1),
    flip(false), val(_val), sum(_val){}
  };
  node *root;
  int size(node *v){
    return !v ? 0 : v->sz;
  }
  void update(node *v){
    v->sz = 1;
    v->sum = v->val;
    if(v->l){
      v->sz += v->l->sz;
      v->sum = merge_val(v->l->sum, v->sum);
    }
    if(v->r){
      v->sz += v->r->sz;
      v->sum = merge_val(v->sum, v->r->sum);
    }
  }
  void push_down(node *v){
    if(!v) return;
    if(v->flip){
      if(v->l) flip(v->l);
      if(v->r) flip(v->r);
      v->flip = false;
    }
  }
  void flip(node *v){
    std::swap(v->l, v->r);
    v->sum = flip_val(v->sum);
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
  node *query_inner(node *v, int l, int r, Val &res){
    if(r == l) return v;
    auto [a, b, c] = split3(v, l, r);
    res = b->sum;
    return merge3(a, b, c);
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
  template<typename F>
  int bisect_from_left(node *v, F &f, Val &ok){
    if(!v) return -1;
    push_down(v);
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    Val m = merge_val(ok, v->sum);
    if(!f(m)){
      ok = m;
      return -1;
    }
    int x = bisect_from_left(v->l, f, ok);
    if(x != -1) return x;
    ok = merge_val(ok, v->val);
    if(f(ok)) return szl;
    int res = bisect_from_left(v->r, f, ok);
    return res == -1 ? res : res + szm;
  }
  template<typename F>
  int bisect_from_right(node *v, F &f, Val &ok){
    if(!v) return -1;
    push_down(v);
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    Val m = merge_val(ok, v->sum);
    if(!f(m)){
      ok = m;
      return -1;
    }
    int x = bisect_from_right(v->r, f, ok);
    if(x != -1) return x + szm;
    ok = merge_val(ok, v->val);
    if(f(ok)) return szl;
    return bisect_from_right(v->l, f, ok);
  }
  dynamic_sequence_merge_split(node *r): root(r){}
public:
  dynamic_sequence_merge_split(): root(nullptr){}
  dynamic_sequence_merge_split(const std::vector<Val> &v): root(nullptr){
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
    Val res = id();
    root = get_inner(root, k, res);
    return res;
  }
  Val query(int l, int r){
    assert(0 <= l && r <= size());
    Val res = id();
    root = query_inner(root, l, r, res);
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
  // f(sum[l, r])がtrueになる最左のr. ない場合は-1
  template<typename F>
  int bisect_from_left(int l, F f){
    if(l >= size()) return -1;
    Val ok = id();
    if(!l){
      int i = bisect_from_left(root, f, ok);
      if(i != -1) root = splay_top_down(root, i);
      return i;
    }
    root = splay_top_down(root, l - 1);
    int i = bisect_from_left(root->r, f, ok);
    if(i != -1){
      i += l;
      root = splay_top_down(root, i);
    }
    return i;
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合は-1
  template<typename F>
  int bisect_from_right(int r, F f){
    if(r < 0) return -1;
    Val ok = id();
    if(r + 1 == size()){
      int i = bisect_from_right(root, f, ok);
      if(i != -1) root = splay_top_down(root, i);
      return i;
    }
    root = splay_top_down(root, r + 1);
    int i = bisect_from_right(root->l, f, ok);
    if(i != -1) root = splay_top_down(root, i);
    return i;
  }
};

template<typename monoid>
struct lazy_dynamic_sequence_merge_split{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge_val = monoid::merge;
  static constexpr auto flip_val = monoid::flip;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
private:
  struct node{
    node *l, *r;
    int sz;
    bool flip;
    Val val, sum;
    Lazy lazy;
    node(Val _val = id()): l(nullptr), r(nullptr), sz(1),
    flip(false), val(_val), sum(_val), lazy(id_lazy()){}
  };
  node *root;
  int size(node *v){
    return !v ? 0 : v->sz;
  }
  void update(node *v){
    v->sz = 1;
    v->sum = v->val;
    if(v->l){
      v->sz += v->l->sz;
      v->sum = merge_val(v->l->sum, v->sum);
    }
    if(v->r){
      v->sz += v->r->sz;
      v->sum = merge_val(v->sum, v->r->sum);
    }
  }
  void push_down(node *v){
    if(!v) return;
    if(v->lazy != id_lazy()){
      if(v->l) propagate(v->l, v->lazy);
      if(v->r) propagate(v->r, v->lazy);
      v->lazy = id_lazy();
    }
    if(v->flip){
      if(v->l) flip(v->l);
      if(v->r) flip(v->r);
      v->flip = false;
    }
  }
  void propagate(node *v, Lazy x){
    v->lazy = propagate_lazy(v->lazy, x);
    v->val = apply(v->val, x, 0, 1);
    v->sum = apply(v->sum, x, 0, v->sz);
  }
  void flip(node *v){
    std::swap(v->l, v->r);
    v->sum = flip_val(v->sum);
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
  node *update_inner(node *v, int l, int r, Lazy x){
    if(r == l) return v;
    auto [a, b, c] = split3(v, l, r);
    propagate(b, x);
    return merge3(a, b, c);
  }
  node *query_inner(node *v, int l, int r, Val &res){
    if(r == l) return v;
    auto [a, b, c] = split3(v, l, r);
    res = b->sum;
    return merge3(a, b, c);
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
  template<typename F>
  int bisect_from_left(node *v, F &f, Val &ok){
    if(!v) return -1;
    push_down(v);
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    Val m = merge_val(ok, v->sum);
    if(!f(m)){
      ok = m;
      return -1;
    }
    int x = bisect_from_left(v->l, f, ok);
    if(x != -1) return x;
    ok = merge_val(ok, v->val);
    if(f(ok)) return szl;
    int res = bisect_from_left(v->r, f, ok);
    return res == -1 ? res : res + szm;
  }
  template<typename F>
  int bisect_from_right(node *v, F &f, Val &ok){
    if(!v) return -1;
    push_down(v);
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    Val m = merge_val(ok, v->sum);
    if(!f(m)){
      ok = m;
      return -1;
    }
    int x = bisect_from_right(v->r, f, ok);
    if(x != -1) return x + szm;
    ok = merge_val(ok, v->val);
    if(f(ok)) return szl;
    return bisect_from_right(v->l, f, ok);
  }
  lazy_dynamic_sequence_merge_split(node *r): root(r){}
public:
  lazy_dynamic_sequence_merge_split(): root(nullptr){}
  lazy_dynamic_sequence_merge_split(const std::vector<Val> &v): root(nullptr){
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
    Val res = id();
    root = get_inner(root, k, res);
    return res;
  }
  void update(int l, int r, Lazy x){
    assert(0 <= l && r <= size());
    root = update_inner(root, l, r, x);
  }
  Val query(int l, int r){
    assert(0 <= l && r <= size());
    Val res = id();
    root = query_inner(root, l, r, res);
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
  // f(sum[l, r])がtrueになる最左のr. ない場合は-1
  template<typename F>
  int bisect_from_left(int l, F f){
    if(l >= size()) return -1;
    Val ok = id();
    if(!l){
      int i = bisect_from_left(root, f, ok);
      if(i != -1) root = splay_top_down(root, i);
      return i;
    }
    root = splay_top_down(root, l - 1);
    int i = bisect_from_left(root->r, f, ok);
    if(i != -1){
      i += l;
      root = splay_top_down(root, i);
    }
    return i;
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合は-1
  template<typename F>
  int bisect_from_right(int r, F f){
    if(r < 0) return -1;
    Val ok = id();
    if(r + 1 == size()){
      int i = bisect_from_right(root, f, ok);
      if(i != -1) root = splay_top_down(root, i);
      return i;
    }
    root = splay_top_down(root, r + 1);
    int i = bisect_from_right(root->l, f, ok);
    if(i != -1) root = splay_top_down(root, i);
    return i;
  }
};
#endif
