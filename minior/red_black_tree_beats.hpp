#ifndef _RED_BLACK_TREE_BEATS_H_
#define _RED_BLACK_TREE_BEATS_H_
#include <vector>
#include <cassert>

template<typename beats_struct>
struct red_black_tree_beats{
  using Val = typename beats_struct::Val;
  using Lazy = typename beats_struct::Lazy;
  struct node{
    node *l, *r, *p;
    bool red;
    int ra, sz;
    Val val;
    Lazy lazy;
    bool is_id;
    // 葉
    node(Val val): l(nullptr), r(nullptr), p(nullptr), red(false), ra(0), sz(1), val(val), is_id(true){}
    // 中間ノード
    node(node *l, node *r, bool red): l(l), r(r), p(nullptr), red(red), ra(l->ra + !(l->red)), sz(l->sz + r->sz), is_id(true){
      l->p = r->p = this;
      beats_struct::merge_val(val, l->val, r->val);
    }
  };
  static std::vector<node*> stock;
  static node *reuse(node *a, node *l, node *r, bool red){
    a->l = l, a->r = r, a->p = nullptr, a->red = red;
    l->p = r->p = a;
    a->ra = l->ra + !(l->red);
    a->sz = l->sz + r->sz;
    beats_struct::merge_val(a->val, l->val, r->val);
    a->is_id = true;
    return a;
  }
  static void propagate(node *a, Lazy &x){
    beats_struct::apply(a->val, x, 0, a->sz);
    beats_struct::propagate_lazy(a->lazy, x);
  }
  static void push_down(node *a){
    if(a->is_id) return;
    if(a->l) propagate(a->l, a->lazy);
    if(a->r) propagate(a->r, a->lazy);
    a->is_id = true;
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
  void set(node *a, int k, Val &x){
    if(!a->l && !a->r){
      assert(k == 0);
      a->val = x;
      return;
    }
    push_down(a);
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) set(a->l, k, x);
    else set(a->r, k - szl, x);
    beats_struct::merge_val(a->val, a->l->val, a->r->val);
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
  template<int id>
  static void update(node *a, int l, int r, Lazy &x){
    if(!a || l >= r || a->sz <= l || r <= 0 || beats_struct::template break_check<id>(a->val, x)) return;
    if(l <= 0 && a->sz <= r && beats_struct::template tag_check<id>(a->val, x)){
      propagate(a, x);
      return;
    }
    if(!a->l && !a->r){
      beats_struct::apply(a->val, x, 0, 1);
      return;
    }
    push_down(a);
    int szl = a->l->sz;
    update(a->l, l, r, x);
    update(a->r, l - szl, r - szl, x);
    beats_struct::merge_val(a->val, a->l->val, a->r->val);
  }
  static void query(node *a, int l, int r, Val &ans){
    if(!a || l >= r || a->sz <= l || r <= 0) return;
    if((l <= 0 && a->sz <= r) || (!a->l && !a->r)){
      Val tmp = beats_struct::id_val();
      beats_struct::merge_val(tmp, ans, a->val);
      ans = tmp;
      return;
    }
    push_down(a);
    int szl = a->l->sz;
    query(a->l, l, r, ans);
    query(a->r, l - szl, r - szl, ans);
  }
  template<typename T>
  node *build(const std::vector<T> &v, int l, int r){
    if(l == r) return nullptr;
    if(r - l == 1) return new node(Val(v[l]));
    int mid = (l + r) / 2;
    node *L = build<T>(v, l, mid);
    node *R = build<T>(v, mid, r);
    return merge(L, R);
  }
  red_black_tree_beats(node *a): root(a){}
  red_black_tree_beats(): root(nullptr){}
  template<typename T>
  red_black_tree_beats(const std::vector<T> &v): root(build(v, 0, v.size())){}
  int size()const{
    return root ? root->sz : 0;
  }
  template<typename T>
  void set(int k, T x){
    assert(k < size());
    set(root, k, Val(x));
  }
  Val get(int k){
    assert(k < size());
    return get(root, k);
  }
  // k番目にxを挿入
  template<typename T>
  void insert(int k, T x){
    auto [a, b] = split(root, k);
    root = merge(a, merge(new node(Val(x)), b));
  }
  void erase(int k){
    assert(root->sz > k);
    auto [a, b] = split(root, k);
    assert(b);
    auto [b2, c] = split(b, 1);
    root = merge(a, c);
    if(b2) stock.push_back(b2);
  }
  template<int id>
  void update(int l, int r, Lazy x){
    assert(0 <= l && r <= size());
    return update<id>(root, l, r, x);
  }
  Val query(int l, int r)const{
    assert(0 <= l && r <= size());
    Val res = beats_struct::id_val();
    query(root, l, r, res);
    return res;
  }
  template<typename T>
  void update_all(T x){
    if(root) propagate(root, Lazy(x));
  }
  Val query_all(){
    return !root ? beats_struct::id_val() : root->val;
  }
  using rbt = red_black_tree_beats<beats_struct>;
  std::pair<rbt, rbt> split(int k){
    return split(*this, k);
  }
  // 2つに分割. 永続でないためこれ自身のaのrootはnullptrになる
  static std::pair<rbt, rbt> split(rbt &a, int k){
    assert(k <= a.size());
    auto [l, r] = split(a.root, k);
    a.root = nullptr;
    return {rbt(l), rbt(r)};
  }
  rbt merge(rbt &b){
    rbt res = merge(*this, b);
    root = b.root = nullptr;
    return res;
  }
  // a, bをマージ. 永続でないためa, bのrootはnullptrになる
  static rbt merge(rbt &a, rbt &b){
    rbt res(merge(a.root, b.root));
    a.root = b.root = nullptr;
    return res;
  }
};
template<typename beats_struct>
std::vector<typename red_black_tree_beats<beats_struct>::node*> red_black_tree_beats<beats_struct>::stock;
#endif