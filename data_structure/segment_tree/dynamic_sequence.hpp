#ifndef _DYNAMIC_SEQUENCE_H_
#define _DYNAMIC_SEQUENCE_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include "../../algebraic_structure/monoid.hpp"

template<typename monoid>
struct dynamic_sequence{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
private:
  struct node{
    int h, sz;
    Val val, sum;
    node *l, *r;
    node(Val _val = id()): h(1), sz(1), val(_val), sum(val), l(nullptr), r(nullptr){}
    int balanace_factor(){return (l ? l->h : 0) - (r ? r->h : 0);}
  };
  node *root, *tmp_node;
  static int size(node *v){return v ? v->sz : 0;}
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz = 1;
    v->sum = v->val;
    if(v->l){
      v->sz += v->l->sz;
      v->sum = merge(v->l->sum, v->sum);
    }
    if(v->r){
      v->sz += v->r->sz;
      v->sum = merge(v->sum, v->r->sum);
    }
  }
  
  static node *rotate_right(node *v){
    node *l = v->l;
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  static node *rotate_left(node *v){
    node *r = v->r;
    v->r = r->l;
    r->l = v;
    update(v);
    update(r);
    return r;
  }
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  node *balance(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    if(bf == 2){
      if(v->l->balanace_factor() == -1){
        v->l = rotate_left(v->l);
        update(v);
      }
      return rotate_right(v);
    }else if(bf == -2){
      if(v->r->balanace_factor() == 1){
        v->r = rotate_right(v->r);
        update(v);
      }
      return rotate_left(v);
    }
    return v;
  }
  Val query_inner(node *v, int l, int r){
    if(!v) return id();
    if(l == 0 && r == v->sz) return v->sum;
    int szl = size(v->l), szv = szl + 1;
    Val left_q = id(), right_q = left_q;
    if(l < szl){
      if(r <= szl) return query_inner(v->l, l, r);
      left_q = query_inner(v->l, l, szl);
      l = szl;
    }
    if(szv < r){
      if(szv <= l) return query_inner(v->r, l - szv, r - szv);
      right_q = query_inner(v->r, 0, r - szv);
      r = szv;
    }
    Val res = (l == r ? id() : v->val);
    res = merge(left_q, res);
    res = merge(res, right_q);
    return res;
  }
  node *cut_left_most(node *v){
    if(v->l){
      v->l = cut_left_most(v->l);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  node *cut_right_most(node *v){
    if(v->r){
      v->r = cut_right_most(v->r);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->l;
  }
  void set_inner(node *v, int k, Val x){
    int szl = v->l ? v->l->sz : 0;
    if(k < szl) set_inner(v->l, k, x);
    else if(k > szl) set_inner(v->r, k - szl - 1, x);
    else v->val = x;
    update(v);
  }
  Val get_inner(node *v, int k){
    int szl = v->l ? v->l->sz : 0;
    if(k < szl) return get_inner(v->l, k);
    else if(k > szl) return get_inner(v->r, k - szl - 1);
    else return v->val;
  }
  node *insert_inner(node *v, int k, Val x){
    assert(size() >= k);
    if(!v) return new node(x);
    int szl = v->l ? v->l->sz : 0;
    if(k <= szl) v->l = insert_inner(v->l, k, x);
    else if(k > szl) v->r = insert_inner(v->r, k - szl - 1, x);
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, int k){
    assert(0 <= k && k < size());
    int szl = v->l ? v->l->sz : 0;
    if(k < szl) v->l = erase_inner(v->l, k);
    else if(k > szl) v->r = erase_inner(v->r, k - szl - 1);
    else{
      if(!v->r) return v->l;
      node *u = cut_left_most(v->r);
      tmp_node->l = v->l;
      tmp_node->r = u;
      v = tmp_node;
    }
    update(v);
    return balance(v);
  }
  // https://atcoder.jp/contests/practice2/submissions/42114957
  template<typename F>
  int bisect_from_left(node *v, F &f, int l, Val &ok){
    if(!v || v->sz <= l) return -1;
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    Val m = merge(ok, v->sum);
    if(!l && !f(m)){
      ok = m;
      return -1;
    }
    int x = bisect_from_left(v->l, f, l, ok);
    if(x != -1) return x;
    if(l <= szl){
      ok = merge(ok, v->val);
      if(f(ok)) return szl;
    }
    int res = bisect_from_left(v->r, f, std::max(l - szm, 0), ok);
    return res == -1 ? res : res + szm;
  }
  // https://atcoder.jp/contests/practice2/submissions/42115281
  template<typename F>
  int bisect_from_right(node *v, F &f, int r, Val &ok){
    if(!v || r < 0) return -1;
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    Val m = merge(ok, v->sum);
    if(v->sz <= r && !f(m)){
      ok = m;
      return -1;
    }
    int x = bisect_from_right(v->r, f, r - szm, ok);
    if(x != -1) return x + szm;
    if(szl <= r){
      ok = merge(ok, v->val);
      if(f(ok)) return szl;
    }
    return bisect_from_right(v->l, f, r, ok);
  }
public:
  dynamic_sequence(): root(nullptr){}
  dynamic_sequence(const std::vector<Val> &v){
    if(v.empty()){
      root = nullptr;
      return;
    }
    int n = v.size();
    std::vector<node*> nodes(n);
    for(int i = 0; i < n; i++) nodes[i] = new node(v[i]);
    root = build(nodes, 0, n);
  }
  int size(){
    return size(root);
  }
  void set(int k, Val x){
    assert(0 <= k && k < size());
    set_inner(root, k, x);
  }
  Val get(int k){
    assert(0 <= k && k < size());
    return get_inner(root, k);
  }
  void insert(int k, Val x){
    assert(0 <= k && k <= size());
    root = insert_inner(root, k, x);
  }
  void erase(int k){
    assert(0 <= k && k < size());
    root = erase_inner(root, k);
  }
  Val query(int l, int r){
    if(l >= r) return id();
    assert(0 <= l && r <= size());
    return query_inner(root, l, r);
  }
  Val query_all(){
    if(!root) return id();
    return root->sum;
  }
  // f(sum[l, r])がtrueになる最左のr. ない場合は-1
  template<typename F>
  int bisect_from_left(int l, const F &f){
    Val x = id();
    return bisect_from_left(root, f, l, x);
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合は-1
  template<typename F>
  int bisect_from_right(int r, const F &f){
    Val x = id();
    return bisect_from_right(root, f, r, x);
  }
};

template<typename monoid>
struct lazy_dynamic_sequence{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
private:
  struct node{
    int h, sz;
    Val val, sum;
    Lazy lazy;
    node *l, *r;
    node(Val _val = id()): h(1), sz(1), val(_val), sum(val), lazy(id_lazy()), l(nullptr), r(nullptr){}
    int balanace_factor(){return (l ? l->h : 0) - (r ? r->h : 0);}
  };
  node *root, *tmp_node;
  static int size(node *v){return v ? v->sz : 0;}
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz = 1;
    v->sum = v->val;
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
    v->lazy = propagate_lazy(v->lazy, x);
    v->val = apply(v->val, x, 0, 1);
    v->sum = apply(v->sum, x, 0, v->sz);
  }
  static void push_down(node *v){
    if(v->lazy != id_lazy()){
      if(v->l) propagate(v->l, v->lazy);
      if(v->r) propagate(v->r, v->lazy);
      v->lazy = id_lazy();
    }
  }
  static node *rotate_right(node *v){
    push_down(v);
    node *l = v->l;
    push_down(l);
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  static node *rotate_left(node *v){
    push_down(v);
    node *r = v->r;
    push_down(r);
    v->r = r->l;
    r->l = v;
    update(v);
    update(r);
    return r;
  }
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  node *balance(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    push_down(v);
    if(bf == 2){
      if(v->l->balanace_factor() == -1){
        v->l = rotate_left(v->l);
        update(v);
      }
      return rotate_right(v);
    }else if(bf == -2){
      if(v->r->balanace_factor() == 1){
        v->r = rotate_right(v->r);
        update(v);
      }
      return rotate_left(v);
    }
    return v;
  }
  void set_inner(node *v, int k, Val x){
    push_down(v);
    int szl = v->l ? v->l->sz : 0;
    if(k < szl) set_inner(v->l, k, x);
    else if(k > szl) set_inner(v->r, k - szl - 1, x);
    else v->val = x;
    update(v);
  }
  Val get_inner(node *v, int k){
    push_down(v);
    int szl = v->l ? v->l->sz : 0;
    if(k < szl) return get_inner(v->l, k);
    else if(k > szl) return get_inner(v->r, k - szl - 1);
    else return v->val;
  }
  void update_inner(node *v, int l, int r, Lazy x){
    if(!v) return;
    if(l == 0 && r == v->sz){
      propagate(v, x);
      return;
    }
    push_down(v);
    int szl = size(v->l), szv = szl + 1;
    if(l < szl){
      if(r <= szl){
        update_inner(v->l, l, r, x);
        update(v);
        return;
      }
      update_inner(v->l, l, szl, x);
      l = szl;
    }
    if(szv < r){
      if(szv <= l){
        update_inner(v->r, l - szv, r - szv, x);
        update(v);
        return;
      }
      update_inner(v->r, 0, r - szv, x);
      r = szv;
    }
    if(r == l + 1) v->val = aapply(v->val, x, 0, 1);
    update(v);
  }
  Val query_inner(node *v, int l, int r){
    if(!v) return id();
    if(l == 0 && r == v->sz) return v->sum;
    push_down(v);
    int szl = size(v->l), szv = szl + 1;
    Val left_q = id(), right_q = left_q;
    if(l < szl){
      if(r <= szl) return query_inner(v->l, l, r);
      left_q = query_inner(v->l, l, szl);
      l = szl;
    }
    if(szv < r){
      if(szv <= l) return query_inner(v->r, l - szv, r - szv);
      right_q = query_inner(v->r, 0, r - szv);
      r = szv;
    }
    Val res = (l == r ? id() : v->val);
    res = merge(left_q, res);
    res = merge(res, right_q);
    return res;
  }
  node *cut_left_most(node *v){
    push_down(v);
    if(v->l){
      v->l = cut_left_most(v->l);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  node *cut_right_most(node *v){
    push_down(v);
    if(v->r){
      v->r = cut_right_most(v->r);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->l;
  }
  node *insert_inner(node *v, int k, Val x){
    if(!v) return new node(x);
    push_down(v);
    int szl = v->l ? v->l->sz : 0;
    if(k <= szl) v->l = insert_inner(v->l, k, x);
    else if(k > szl) v->r = insert_inner(v->r, k - szl - 1, x);
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, int k){
    push_down(v);
    int szl = v->l ? v->l->sz : 0;
    if(k < szl) v->l = erase_inner(v->l, k);
    else if(k > szl) v->r = erase_inner(v->r, k - szl - 1);
    else{
      if(!v->r) return v->l;
      node *u = cut_left_most(v->r);
      tmp_node->l = v->l;
      tmp_node->r = u;
      v = tmp_node;
    }
    update(v);
    return balance(v);
  }
  // https://atcoder.jp/contests/practice2/submissions/42114957
  template<typename F>
  int bisect_from_left(node *v, F &f, int l, Val &ok){
    if(!v || v->sz <= l) return -1;
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    Val m = merge(ok, v->sum);
    if(!l && !f(m)){
      ok = m;
      return -1;
    }
    push_down(v);
    int x = bisect_from_left(v->l, f, l, ok);
    if(x != -1) return x;
    if(l <= szl){
      ok = merge(ok, v->val);
      if(f(ok)) return szl;
    }
    int res = bisect_from_left(v->r, f, std::max(l - szm, 0), ok);
    return res == -1 ? res : res + szm;
  }
  // https://atcoder.jp/contests/practice2/submissions/42115281
  template<typename F>
  int bisect_from_right(node *v, F &f, int r, Val &ok){
    if(!v || r < 0) return -1;
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    Val m = merge(ok, v->sum);
    if(v->sz <= r && !f(m)){
      ok = m;
      return -1;
    }
    push_down(v);
    int x = bisect_from_right(v->r, f, r - szm, ok);
    if(x != -1) return x + szm;
    if(szl <= r){
      ok = merge(ok, v->val);
      if(f(ok)) return szl;
    }
    return bisect_from_right(v->l, f, r, ok);
  }
public:
  lazy_dynamic_sequence(): root(nullptr){}
  lazy_dynamic_sequence(const std::vector<Val> &v){
    if(v.empty()){
      root = nullptr;
      return;
    }
    int n = v.size();
    std::vector<node*> nodes(n);
    for(int i = 0; i < n; i++) nodes[i] = new node(v[i]);
    root = build(nodes, 0, n);
  }
  int size(){
    return size(root);
  }
  void set(int k, Val x){
    assert(0 <= k && k < size());
    set_inner(root, k, x);
  }
  Val get(int k){
    assert(0 <= k && k < size());
    return get_inner(root, k);
  }
  void insert(int k, Val x){
    assert(0 <= k && k <= size());
    root = insert_inner(root, k, x);
  }
  void erase(int k){
    assert(0 <= k && k < size());
    root = erase_inner(root, k);
  }
  void update(int l, int r, Lazy x){
    if(l >= r) return;
    assert(0 <= l && r <= size());
    update_inner(root, l, r, x);
  }
  Val query(int l, int r){
    if(l >= r) return id();
    assert(0 <= l && r <= size());
    return query_inner(root, l, r);
  }
  Val query_all(){
    if(!root) return id();
    return root->sum;
  }
  // f(sum[l, r])がtrueになる最左のr. ない場合は-1
  template<typename F>
  int bisect_from_left(int l, const F &f){
    Val x = id();
    return bisect_from_left(root, f, l, x);
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合は-1
  template<typename F>
  int bisect_from_right(int r, const F &f){
    Val x = id();
    return bisect_from_right(root, f, r, x);
  }
};
#endif