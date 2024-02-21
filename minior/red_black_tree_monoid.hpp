#ifndef _RED_BLACK_TREE_MONOID_H_
#define _RED_BLACK_TREE_MONOID_H_
#include <vector>
#include <cassert>
#include "../algebraic_structure/monoid.hpp"

template<typename monoid>
struct red_black_tree_monoid{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge_val = monoid::merge;
  struct node{
    node *l, *r, *p;
    bool red;
    int ra, sz;
    Val val;
    // 葉
    node(Val val): l(nullptr), r(nullptr), p(nullptr), red(false), ra(0), sz(1), val(val){}
    // 中間ノード
    node(node *l, node *r, bool red): l(l), r(r), p(nullptr), red(red), ra(l->ra + !(l->red)), sz(l->sz + r->sz){
      l->p = r->p = this;
      val = merge_val(l->val, r->val);
    }
  };
  static std::vector<node*> stock;
  static node *reuse(node *a, node *l, node *r, bool red){
    a->l = l, a->r = r, a->p = nullptr, a->red = red;
    l->p = r->p = a;
    a->ra = l->ra + !(l->red);
    a->sz = l->sz + r->sz;
    a->val = merge_val(a->l->val, a->r->val);
    return a;
  }
  node *root;
  using pn = std::pair<node*, node*>;
  static node *merge_sub(node *a, node *b){
    if(a->ra < b->ra){
      node *c = merge_sub(a, b->l);
      if(!b->red && c->red && c->l->red){
        if(!b->r->red){
          return reuse(c, c->l, reuse(b, c->r, b->r, 1), 0);
        }else{
          b->r->red = 0;
          c->red = 0;
          return reuse(b, c, b->r, 1);
        }
      }else{
        return reuse(b, c, b->r, b->red);
      }
    }else if(a->ra > b->ra){
      node *c = merge_sub(a->r, b);
      if(!a->red && c->red && c->r->red){
        if(!a->l->red){
          return reuse(c, reuse(a, a->l, c->l, 1), c->r, 0);
        }else{
          a->l->red = 0;
          c->red = 0;
          return reuse(a, a->l, c, 1);
        }
      }else{
        return reuse(a, a->l, c, a->red);
      }
    }else{
      if(stock.empty()) return new node(a, b, 1);
      node *u = stock.back();
      stock.pop_back();
      return reuse(u, a, b, 1);
    }
  }
  static node *merge(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    node *c = merge_sub(a, b);
    c->red = 0;
    return c;
  }
  static node *as_root(node *a){
    if(!a) return nullptr;
    a->red = 0;
    return a;
  }
  static pn split(node *a, int k){
    int sz = a->sz, szl = (a->l ? a->l->sz : 0);
    if(k == 0 || k == sz) return (!k ? pn{nullptr, a} : pn{a, nullptr});
    pn res;
    if(k < szl){
      auto [l, r] = split(a->l, k);
      res = pn{l, merge(r, as_root(a->r))};
    }else if(k > szl){
      auto [l, r] = split(a->r, k - szl);
      res = pn{merge(as_root(a->l), l), r};
    }else{
      res = pn{as_root(a->l), as_root(a->r)};
    }
    if(a) stock.push_back(a);
    return res;
  }
  void set(node *a, int k, Val x){
    if(!a->l && !a->r){
      assert(k == 0);
      a->val = x;
      return;
    }
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) set(a->l, k, x);
    else set(a->r, k - szl, x);
    a->val = merge_val(a->l->val, a->r->val);
  }
  Val get(node *a, int k){
    if(!a->l && !a->r){
      assert(k == 0);
      return a->val;
    }
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) return get(a->l, k);
    else return get(a->r, k - szl);
  }
  static Val query(node *a, int l, int r){
    if(!a || l >= r || a->sz <= l || r <= 0) return id();
    if(l <= 0 && a->sz <= r) return a->val;
    if(!a->l && !a->r) return a->val;
    int szl = a->l->sz;
    if(r <= szl) return query(a->l, l, r);
    if(szl <= l) return query(a->r, l - szl, r - szl);
    return merge_val(query(a->l, l, szl), query(a->r, 0, r - szl));
  }
  node *build(const std::vector<Val> &v, int l, int r){
    if(l == r) return nullptr;
    if(r - l == 1) return new node(v[l]);
    int mid = (l + r) / 2;
    node *L = build(v, l, mid);
    node *R = build(v, mid, r);
    return merge(L, R);
  }
  red_black_tree_monoid(node *a): root(a){}
  red_black_tree_monoid(): root(nullptr){}
  red_black_tree_monoid(const std::vector<Val> &v): root(build(v, 0, v.size())){}

  int size()const{
    return root ? root->sz : 0;
  }
  void set(int k, Val x){
    assert(k < size());
    set(root, k, x);
  }
  Val get(int k){
    assert(k < size());
    return get(root, k);
  }
  // k番目にxを挿入
  void insert(int k, Val x){
    auto [a, b] = split(root, k);
    root = merge(a, merge(new node(x), b));
  }
  void erase(int k){
    assert(root->sz > k);
    auto [a, b] = split(root, k);
    assert(b);
    auto [b2, c] = split(b, 1);
    root = merge(a, c);
    if(b2) stock.push_back(b2);
  }
  Val query(int l, int r)const{
    assert(0 <= l && r <= size());
    return query(root, l, r);
  }
  Val query_all(){
    return !root ? id() : root->val;
  }
  using rbtm = red_black_tree_monoid<monoid>;
  std::pair<rbtm, rbtm> split(int k){
    return split(*this, k);
  }
  // 2つに分割. 永続でないためこれ自身のaのrootはnullptrになる
  static std::pair<rbtm, rbtm> split(rbtm &a, int k){
    assert(k <= a.size());
    auto [l, r] = split(a.root, k);
    a.root = nullptr;
    return {rbtm(l), rbtm(r)};
  }
  // a, bをマージ. 永続でないためa, bのrootはnullptrになる
  static rbtm merge(rbtm &a, rbtm &b){
    rbtm res(merge(a.root, b.root));
    a.root = b.root = nullptr;
    return res;
  }
};
template<typename monoid>
std::vector<typename red_black_tree_monoid<monoid>::node*> red_black_tree_monoid<monoid>::stock;



template<typename monoid>
struct lazy_red_black_tree_monoid{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge_val = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
  struct node{
    node *l, *r, *p;
    bool red;
    int ra, sz;
    Val val;
    Lazy lazy;
    // 葉
    node(Val val): l(nullptr), r(nullptr), p(nullptr), red(false), ra(0), sz(1), val(val), lazy(id_lazy()){}
    // 中間ノード
    node(node *l, node *r, bool red): l(l), r(r), p(nullptr), red(red), ra(l->ra + !(l->red)), sz(l->sz + r->sz), lazy(id_lazy()){
      l->p = r->p = this;
      val = merge_val(l->val, r->val);
    }
  };
  static std::vector<node*> stock;
  static node *reuse(node *a, node *l, node *r, bool red){
    a->l = l, a->r = r, a->p = nullptr, a->red = red;
    l->p = r->p = a;
    a->ra = l->ra + !(l->red);
    a->sz = l->sz + r->sz;
    a->val = merge_val(a->l->val, a->r->val);
    a->lazy = id_lazy();
    return a;
  }
  static void propagate(node *a, Lazy x){
    a->val = apply(a->val, x, 0, a->sz);
    a->lazy = propagate_lazy(a->lazy, x);
  }
  static void push_down(node *a){
    if(a->lazy == id_lazy()) return;
    if(a->l) propagate(a->l, a->lazy);
    if(a->r) propagate(a->r, a->lazy);
    a->lazy = id_lazy();
  }
  node *root;
  using pn = std::pair<node*, node*>;
  static node *merge_sub(node *a, node *b){
    if(a->ra < b->ra){
      push_down(b);
      node *c = merge_sub(a, b->l);
      if(!b->red && c->red && c->l->red){
        if(!b->r->red){
          return reuse(c, c->l, reuse(b, c->r, b->r, 1), 0);
        }else{
          b->r->red = 0;
          c->red = 0;
          return reuse(b, c, b->r, 1);
        }
      }else{
        return reuse(b, c, b->r, b->red);
      }
    }else if(a->ra > b->ra){
      push_down(a);
      node *c = merge_sub(a->r, b);
      if(!a->red && c->red && c->r->red){
        if(!a->l->red){
          return reuse(c, reuse(a, a->l, c->l, 1), c->r, 0);
        }else{
          a->l->red = 0;
          c->red = 0;
          return reuse(a, a->l, c, 1);
        }
      }else{
        return reuse(a, a->l, c, a->red);
      }
    }else{
      if(stock.empty()) return new node(a, b, 1);
      node *u = stock.back();
      stock.pop_back();
      return reuse(u, a, b, 1);
    }
  }
  static node *merge(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    node *c = merge_sub(a, b);
    c->red = 0;
    return c;
  }
  static node *as_root(node *a){
    if(!a) return nullptr;
    a->red = 0;
    return a;
  }
  static pn split(node *a, int k){
    int sz = a->sz, szl = (a->l ? a->l->sz : 0);
    if(k == 0 || k == sz) return (!k ? pn{nullptr, a} : pn{a, nullptr});
    pn res;
    push_down(a);
    if(k < szl){
      auto [l, r] = split(a->l, k);
      res = pn{l, merge(r, as_root(a->r))};
    }else if(k > szl){
      auto [l, r] = split(a->r, k - szl);
      res = pn{merge(as_root(a->l), l), r};
    }else{
      res = pn{as_root(a->l), as_root(a->r)};
    }
    if(a) stock.push_back(a);
    return res;
  }
  void set(node *a, int k, Val x){
    if(!a->l && !a->r){
      assert(k == 0);
      a->val = x;
      return;
    }
    push_down(a);
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) set(a->l, k, x);
    else set(a->r, k - szl, x);
    a->val = merge_val(a->l->val, a->r->val);
  }
  Val get(node *a, int k){
    if(!a->l && !a->r){
      assert(k == 0);
      return a->val;
    }
    push_down(a);
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) return get(a->l, k);
    else return get(a->r, k - szl);
  }
  static void update(node *a, int l, int r, Lazy x){
    if(!a || l >= r || a->sz <= l || r <= 0) return;
    if(l <= 0 && a->sz <= r){
      propagate(a, x);
      return;
    }
    if(!a->l && !a->r){
      a->val = apply(a->val, x, 0, 1);
      return;
    }
    push_down(a);
    int szl = a->l->sz;
    update(a->l, l, r, x);
    update(a->r, l - szl, r - szl, x);
    a->val = merge_val(a->l->val, a->r->val);
  }
  static Val query(node *a, int l, int r){
    if(!a || l >= r || a->sz <= l || r <= 0) return id();
    if(l <= 0 && a->sz <= r) return a->val;
    if(!a->l && !a->r) return a->val;
    push_down(a);
    int szl = a->l->sz;
    if(r <= szl) return query(a->l, l, r);
    if(szl <= l) return query(a->r, l - szl, r - szl);
    return merge_val(query(a->l, l, szl), query(a->r, 0, r - szl));
  }
  node *build(const std::vector<Val> &v, int l, int r){
    if(l == r) return nullptr;
    if(r - l == 1) return new node(v[l]);
    int mid = (l + r) / 2;
    node *L = build(v, l, mid);
    node *R = build(v, mid, r);
    return merge(L, R);
  }
  lazy_red_black_tree_monoid(node *a): root(a){}
  lazy_red_black_tree_monoid(): root(nullptr){}
  lazy_red_black_tree_monoid(const std::vector<Val> &v): root(build(v, 0, v.size())){}

  int size()const{
    return root ? root->sz : 0;
  }
  void set(int k, Val x){
    assert(k < size());
    set(root, k, x);
  }
  Val get(int k){
    assert(k < size());
    return get(root, k);
  }
  // k番目にxを挿入
  void insert(int k, Val x){
    auto [a, b] = split(root, k);
    root = merge(a, merge(new node(x), b));
  }
  void erase(int k){
    assert(root->sz > k);
    auto [a, b] = split(root, k);
    assert(b);
    auto [b2, c] = split(b, 1);
    root = merge(a, c);
    if(b2) stock.push_back(b2);
  }
  void update(int l, int r, Lazy x){
    assert(0 <= l && r <= size());
    return update(root, l, r, x);
  }
  Val query(int l, int r)const{
    assert(0 <= l && r <= size());
    return query(root, l, r);
  }
  void update_all(Lazy x){
    if(root) propagate(root, x);
  }
  Val query_all(){
    return !root ? id() : root->val;
  }
  using rbtm = lazy_red_black_tree_monoid<monoid>;
  std::pair<rbtm, rbtm> split(int k){
    return split(*this, k);
  }
  // 2つに分割. 永続でないためこれ自身のaのrootはnullptrになる
  static std::pair<rbtm, rbtm> split(rbtm &a, int k){
    assert(k <= a.size());
    auto [l, r] = split(a.root, k);
    a.root = nullptr;
    return {rbtm(l), rbtm(r)};
  }
  rbtm merge(rbtm &b){
    rbtm res = merge(*this, b);
    root = b.root = nullptr;
    return res;
  }
  // a, bをマージ. 永続でないためa, bのrootはnullptrになる
  static rbtm merge(rbtm &a, rbtm &b){
    rbtm res(merge(a.root, b.root));
    a.root = b.root = nullptr;
    return res;
  }
};
template<typename monoid>
std::vector<typename lazy_red_black_tree_monoid<monoid>::node*> lazy_red_black_tree_monoid<monoid>::stock;
#endif