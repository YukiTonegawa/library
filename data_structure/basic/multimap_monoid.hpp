#ifndef _MULTIMAP_MONOID_H_
#define _MULTIMAP_MONOID_H_
#include <vector>
#include <cassert>
#include <iostream>
#include "../../algebraic_structure/monoid.hpp"

// キーの昇順, キーが同じものは追加した順番の昇順
template<typename Key, typename monoid>
struct multimap_monoid{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
private:
  struct node{
    int h, sz_unique, sz_sum;
    Key key;
    Val val, sum;
    node *l, *r, *p, *lmost, *rmost;
    node(Key _key, Val _val): h(1), sz_unique(1), sz_sum(1), key(_key), val(_val), sum(_val), l(nullptr), r(nullptr), p(nullptr), lmost(this), rmost(this){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root, *tmp_node;
  int size_unique(node *v){return v ? v->sz_unique : 0;}
  int size_sum(node *v){return v ? v->sz_sum : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0, v->r ? v->r->h : 0) + 1;
    v->sz_unique = v->sz_sum = 1;
    v->lmost = v->rmost = v;
    v->sum = v->val;
    if(v->l){
      v->sz_unique += v->l->sz_unique + (v->l->rmost->key != v->key);
      v->sz_sum += v->l->sz_sum; 
      v->lmost = v->l->lmost;
      v->sum = merge(v->l->sum, v->sum);
    }
    if(v->r){
      v->sz_unique += v->r->sz_unique + (v->r->lmost->key != v->key);
      v->sz_sum += v->r->sz_sum;
      v->rmost = v->r->rmost;
      v->sum = merge(v->sum, v->r->sum);
    }
  }
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l && ((v->l = build(nodes, l, m)))) v->l->p = v;
    if(r > m + 1 && ((v->r = build(nodes, m + 1, r)))) v->r->p = v;
    update(v);
    return v;
  }
  node *rotate_right(node *v){
    node *l = v->l;
    if((v->l = l->r)) v->l->p = v;
    l->p = v->p;
    l->r = v, v->p = l;
    update(v);
    update(l);
    return l;
  }
  node *rotate_left(node *v){
    node *r = v->r;
    if((v->r = r->l)) v->r->p = v;
    r->p = v->p;
    r->l = v, v->p = r;
    update(v);
    update(r);
    return r;
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
  node *leftmost(node *v){
    return v ? v->lmost : nullptr;
  }
  node *rightmost(node *v){
    return v ? v->rmost : nullptr;
  }
  // kの中で最も早く追加された要素
  node *find_key(node *v, Key k){
    node *ret = nullptr;
    while(v){
      if(k <= v->key){
        if(k == v->key) ret = v;
        v = v->l;
      }else v = v->r;
    }
    return ret;
  }
  node *find_next_inner(node *v){
    while(v){
      if(v->r) return v->r->lmost;
      v = v->p;
    }
    return nullptr;
  }
  node *find_prev_inner(node *v){
    while(v){
      if(v->l) return v->l->rmost;
      v = v->p;
    }
    return nullptr;
  }
  int low_count_unique_inner(node *v, Key k){
    int ret = 0;
    while(v){
      int szl = size_unique(v->l);
      int f = (v->l && v->l->rmost->key != v->key);
      if(k <= v->key) v = v->l;
      else ret += szl + f, v = v->r;
    }
    return ret;
  }
  int low_count_sum_inner(node *v, Key k){
    int ret = 0;
    while(v){
      int szl = size_sum(v->l);
      if(k <= v->key) v = v->l;
      else ret += szl + 1, v = v->r;
    }
    return ret;
  }
  // uの左側にいくつ要素があるか
  int low_count_node_inner(node *v){
    int ret = 0;
    ret += size_sum(v->l);
    while(v){
      if(v->p && v->p->r == v){
        v = v->p;
        ret += size_sum(v->l) + 1;
      }else v = v->p;
    }
    return ret;
  }
  // k種類目のノードの中で最も早く追加されたもの
  node *find_kth_unique(node *v, int k){
    while(true){
      int szl = size_unique(v->l);
      if(szl <= k){
        int f = (v->l && v->l->rmost->key != v->key);
        if(szl == k && f) return v;
        k -= szl + f;
        v = v->r;
      }else v = v->l;
    }
  }
  // k番目の要素(複数ある場合は最も早く追加されたもの)
  node *find_kth_sum(node *v, int k){
    while(true){
      int szl = size_sum(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  node *cut_leftmost(node *v){
    if(v->l){
      if((v->l = cut_leftmost(v->l))) v->l->p = v;
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  node *cut_rightmost(node *v){
    if(v->r){
      if((v->r = cut_rightmost(v->r))) v->r->p = v;
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->l;
  }
  // k以上のキーを持つ最小要素(複数ある場合は最も早く追加したもの)
  node *lower_bound_inner(node *v, Key k){
    node *ret = nullptr;
    while(v){
      if(k <= v->key){
        ret = v;
        v = v->l;
      }else v = v->r;
    }
    return ret;
  }
  // k以下のキーを持つ最大要素(複数ある場合は最も遅く追加したもの)
  node *lower_bound_rev_inner(node *v, Key k){
    node *ret = nullptr;
    while(v){
      if(k >= v->key){
        ret = v;
        v = v->r;
      }else v = v->l;
    }
    return ret;
  }
  // kの要素がすでにある場合それらの中で最も右側に追加
  node *emplace_inner(node *v, Key k, Val val){
    if(!v) return new node(k, val);
    if(k < v->key  && ((v->l = emplace_inner(v->l, k, val)))) v->l->p = v;
    if(k >= v->key && ((v->r = emplace_inner(v->r, k, val)))) v->r->p = v;
    update(v);
    return balance(v);
  }
  // kの要素がすでにある場合何もしない
  node *emplace_unique_inner(node *v, Key k, Val val){
    if(!v) return new node(k, val);
    if(k < v->key  && ((v->l = emplace_unique_inner(v->l, k, val)))) v->l->p = v;
    if(k > v->key && ((v->r = emplace_unique_inner(v->r, k, val)))) v->r->p = v;
    update(v);
    return balance(v);
  }
  // kの要素がすでにある場合そのうちのいずれかと置き換える
  node *emplace_replace_inner(node *v, Key k, Val val){
    if(!v) return new node(k, val);
    if(k < v->key  && ((v->l = emplace_replace_inner(v->l, k, val)))) v->l->p = v;
    if(k > v->key && ((v->r = emplace_replace_inner(v->r, k, val)))) v->r->p = v;
    if(k == v->key) v->val = val;
    update(v);
    return balance(v);
  }
  // kを消す(複数ある場合は最も早く追加した要素を消す)
  node *erase_inner(node *v, Key k){
    if(!v) return nullptr;
    if(k < v->key && ((v->l = erase_inner(v->l, k)))) v->l->p = v;
    if(k > v->key && ((v->r = erase_inner(v->r, k)))) v->r->p = v;
    if(k == v->key){
      if(v->l && v->l->rmost->key == k){
        if((v->l = erase_inner(v->l, k))) v->l->p = v;
        update(v);
        return balance(v);
      }
      if(v->r){
        v->r = cut_leftmost(v->r);
        if((tmp_node->l = v->l)) tmp_node->l->p = tmp_node;
        if((tmp_node->r = v->r)) tmp_node->r->p = tmp_node;
        update(tmp_node);
        return balance(tmp_node);
      }
      return v->l;
    }
    update(v);
    return balance(v);
  }
  // vを消す
  node *erase_node_inner(node *v){
    assert(v);
    node *p = v->p, *u = v;
    if(v->r){
      v->r = cut_leftmost(v->r);
      if((tmp_node->l = v->l)) tmp_node->l->p = tmp_node;
      if((tmp_node->r = v->r)) tmp_node->r->p = tmp_node;
      update(tmp_node);
      v = balance(tmp_node);
    }else{
      v = v->l;
    }
    while(p){
      if(p->l == u && ((p->l = v))) v->p = p;
      if(p->r == u && ((p->r = v))) v->p = p;
      update(p);
      u = p;
      p = u->p;
      v = balance(u);
    }
    return v;
  }
  Val query_inner(node *v, Key l, Key r){
    if(!v || v->rmost->key < l || r <= v->lmost->key) return id();
    if(l <= v->lmost->key && v->rmost->key < r) return v->sum;
    if(r <= v->key) return query_inner(v->l, l, r);
    if(v->key < l) return query_inner(v->r, l, r);
    return merge(query_inner(v->l, l, r), merge(v->val, query_inner(v->r, l, r)));
  }
  template<typename F>
  node *bisect_from_left(node *v, F &f, Key l, Val &ok){
    if(!v || v->rmost->key < l) return nullptr;
    Val m = merge(ok, v->sum);
    if(l <= v->lmost->key && !f(m)){
      ok = m;
      return nullptr;
    }
    node *x = bisect_from_left(v->l, f, l, ok);
    if(x) return x;
    if(l <= v->key){
      ok = merge(ok, v->val);
      if(f(ok)) return v;
    }
    return bisect_from_left(v->r, f, l, ok);
  }
  template<typename F>
  node *bisect_from_right(node *v, F &f, Key r, Val &ok){
    if(!v || v->lmost->key > r) return nullptr;
    Val m = merge(v->sum, ok);
    if(r >= v->rmost->key && !f(m)){
      ok = m;
      return nullptr;
    }
    node *x = bisect_from_right(v->r, f, r, ok);
    if(x) return x;
    if(v->key <= r){
      ok = merge(v->val, ok);
      if(f(ok)) return v;
    }
    return bisect_from_right(v->l, f, r, ok);
  }
  void to_list_inner(node *v, std::vector<std::pair<Key, Val>> &res){
    if(v->l) to_list_inner(v->l, res);
    res.push_back({v->key, v->val});
    if(v->r) to_list_inner(v->r, res);
  }
public:
  multimap_monoid(): root(nullptr){}
  multimap_monoid(std::vector<std::pair<Key, Val>> v){
    std::sort(v.begin(), v.end());
    init_sorted(v);
  }
  // 1
  // すでにソート済み
  void init_sorted(const std::vector<std::pair<Key, Val>> &v){
    if(v.empty()){
      root = nullptr;
      return;
    }
    int n = v.size();
    std::vector<node*> nodes(n);
    for(int i = 0; i < n; i++) nodes[i] = (new node(v[i].first, v[i].second));
    root = build(nodes, 0, n);
  }
  // 要素数
  int size(){
    return size_sum(root);
  }
  // 種類数
  int size_unique(){
    return size_unique(root);
  }
  bool empty(){
    return size_sum(root) == 0;
  }
  // 1
  void emplace(Key k, Val val){
    root = emplace_inner(root, k, val);
    root->p = nullptr;
  }
  // 1
  void emplace_unique(Key k, Val val){
    root = emplace_unique_inner(root, k, val);
    root->p = nullptr;
  }
  // 1
  void emplace_replace(Key k, Val val){
    root = emplace_replace_inner(root, k, val);
    root->p = nullptr;
  }
  // 2
  // 同じキーを持つ要素が複数ある場合, 最後に追加された要素を消す
  void erase(Key k){
    if((root = erase_inner(root, k))) root->p = nullptr;
  }
  void erase_node(node *v){
    if((root = erase_node_inner(v))) root->p = nullptr;
  }
  void clear(){
    root = nullptr;
  }
  bool contain(Key k){
    return find_key(root, k);
  }
  // 2
  node *find(Key k){
    return find_key(root, k);
  }
  node *find_next(node *v){
    if(!v) return nullptr;
    return find_next_inner(v);
  }
  node *find_prev(node *v){
    if(!v) return nullptr;
    return find_prev_inner(v);
  }
  node *min(){
    return leftmost(root);
  }
  node *max(){
    return rightmost(root);
  }
  int low_count_unique(Key k){
    return low_count_unique_inner(root, k);
  }
  int low_count_sum(Key k){
    return low_count_sum_inner(root, k);
  }
  // 1
  int low_count_node(node *v){
    return low_count_node_inner(v);
  }
  node *kth_smallest_unique(int k){
    if(size_unique() <= k) return nullptr;
    return find_kth_unique(root, k);
  }
  node *kth_smallest_sum(int k){
    if(size() <= k) return nullptr;
    return find_kth_sum(root, k);
  }
  // 1
  node *lower_bound(Key k){
    return lower_bound_inner(root, k);
  }
  // 1
  // k以下の最大要素
  node *lower_bound_rev(Key k){
    return lower_bound_rev_inner(root, k);
  }
  Val query(Key l, Key r){
    return query_inner(root, l, r);
  }
  Val query_all(){
    if(!root) return id();
    return root->sum;
  }
  // 1
  template<typename F>
  std::pair<node*, Val> bisect_from_left(Key l, F f){
    Val ok = id();
    node* ret = bisect_from_left(root, f, l, ok);
    return {ret, ok};
  }
  template<typename F>
  std::pair<node*, Val> bisect_from_right(Key r, F f){
    Val ok = id();
    node* ret = bisect_from_right(root, f, r, ok);
    return {ret, ok};
  }
  std::vector<std::pair<Key, Val>> to_list(){
    std::vector<std::pair<Key, Val>> res;
    if(root) to_list_inner(root, res);
    return res;
  }
};

// キーの昇順, キーが同じものは追加した順番の昇順
// ノードを取得してから値を使う前に区間更新をすると壊れる
template<typename Key, typename monoid>
struct lazy_multimap_monoid{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
private:
  struct node{
    int h, sz_unique, sz_sum;
    Key key;
    Val val, sum;
    Lazy lazy;
    node *l, *r, *p, *lmost, *rmost;
    node(Key _key, Val _val): h(1), sz_unique(1), sz_sum(1), key(_key), val(_val), sum(_val), lazy(id_lazy()), l(nullptr), r(nullptr), p(nullptr), lmost(this), rmost(this){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root, *tmp_node;
  int size_unique(node *v){return v ? v->sz_unique : 0;}
  int size_sum(node *v){return v ? v->sz_sum : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0, v->r ? v->r->h : 0) + 1;
    v->sz_unique = v->sz_sum = 1;
    v->lmost = v->rmost = v;
    v->sum = v->val;
    if(v->l){
      v->sz_unique += v->l->sz_unique + (v->l->rmost->key != v->key);
      v->sz_sum += v->l->sz_sum; 
      v->lmost = v->l->lmost;
      v->sum = merge(v->l->sum, v->sum);
    }
    if(v->r){
      v->sz_unique += v->r->sz_unique + (v->r->lmost->key != v->key);
      v->sz_sum += v->r->sz_sum;
      v->rmost = v->r->rmost;
      v->sum = merge(v->sum, v->r->sum);
    }
  }
  void propagate(node *v, Lazy x){
    v->lazy = propagate_lazy(v->lazy, x);
    v->val = apply(v->val, x, 0, 1);
    v->sum = apply(v->sum, x, 0, v->sz_sum);
  }
  void push_down(node *v){
    if(v->lazy != id_lazy()){
      if(v->l) propagate(v->l, v->lazy);
      if(v->r) propagate(v->r, v->lazy);
      v->lazy = id_lazy();
    }
  }
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l && ((v->l = build(nodes, l, m)))) v->l->p = v;
    if(r > m + 1 && ((v->r = build(nodes, m + 1, r)))) v->r->p = v;
    update(v);
    return v;
  }
  node *rotate_right(node *v){
    push_down(v);
    node *l = v->l;
    push_down(l);
    if((v->l = l->r)) v->l->p = v;
    l->p = v->p;
    l->r = v, v->p = l;
    update(v);
    update(l);
    return l;
  }
  node *rotate_left(node *v){
    push_down(v);
    node *r = v->r;
    push_down(r);
    if((v->r = r->l)) v->r->p = v;
    r->p = v->p;
    r->l = v, v->p = r;
    update(v);
    update(r);
    return r;
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
  node *leftmost(node *v){
    while(v->l){
      push_down(v);
      v = v->l;
    }
    return v;
  }
  node *rightmost(node *v){
    while(v->r){
      push_down(v);
      v = v->r;
    }
    return v;
  }
  // kの中で最も早く追加された要素
  node *find_key(node *v, Key k){
    node *ret = nullptr;
    while(v){
      push_down(v);
      if(k <= v->key){
        if(k == v->key) ret = v;
        v = v->l;
      }else v = v->r;
    }
    return ret;
  }
  node *find_next_inner(node *v){
    while(v){
      if(v->r) return leftmost(v->r);
      v = v->p;
    }
    return nullptr;
  }
  node *find_prev_inner(node *v){
    while(v){
      if(v->l) return rightmost(v->l);
      v = v->p;
    }
    return nullptr;
  }
  int low_count_unique_inner(node *v, Key k){
    int ret = 0;
    while(v){
      int szl = size_unique(v->l);
      int f = (v->l && v->l->rmost->key != v->key);
      if(k <= v->key) v = v->l;
      else ret += szl + f, v = v->r;
    }
    return ret;
  }
  int low_count_sum_inner(node *v, Key k){
    int ret = 0;
    while(v){
      int szl = size_sum(v->l);
      if(k <= v->key) v = v->l;
      else ret += szl + 1, v = v->r;
    }
    return ret;
  }
  // uの左側にいくつ要素があるか
  int low_count_node_inner(node *v){
    int ret = 0;
    ret += size_sum(v->l);
    while(v){
      if(v->p && v->p->r == v){
        v = v->p;
        ret += size_sum(v->l) + 1;
      }else v = v->p;
    }
    return ret;
  }
  // k種類目のノードの中で最も早く追加されたもの
  node *find_kth_unique(node *v, int k){
    while(true){
      push_down(v);
      int szl = size_unique(v->l);
      if(szl <= k){
        int f = (v->l && v->l->rmost->key != v->key);
        if(szl == k && f) return v;
        k -= szl + f;
        v = v->r;
      }else v = v->l;
    }
  }
  // k番目の要素(複数ある場合は最も早く追加されたもの)
  node *find_kth_sum(node *v, int k){
    while(true){
      push_down(v);
      int szl = size_sum(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  node *cut_leftmost(node *v){
    push_down(v);
    if(v->l){
      if((v->l = cut_leftmost(v->l))) v->l->p = v;
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  node *cut_rightmost(node *v){
    push_down(v);
    if(v->r){
      if((v->r = cut_rightmost(v->r))) v->r->p = v;
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->l;
  }
  // k以上のキーを持つ最小要素(複数ある場合は最も早く追加したもの)
  node *lower_bound_inner(node *v, Key k){
    node *ret = nullptr;
    while(v){
      push_down(v);
      if(k <= v->key){
        ret = v;
        v = v->l;
      }else v = v->r;
    }
    return ret;
  }
  // k以下のキーを持つ最大要素(複数ある場合は最も遅く追加したもの)
  node *lower_bound_rev_inner(node *v, Key k){
    node *ret = nullptr;
    while(v){
      push_down(v);
      if(k >= v->key){
        ret = v;
        v = v->r;
      }else v = v->l;
    }
    return ret;
  }
  // kの要素がすでにある場合それらの中で最も右側に追加
  node *emplace_inner(node *v, Key k, Val val){
    if(!v) return new node(k, val);
    push_down(v);
    if(k < v->key  && ((v->l = emplace_inner(v->l, k, val)))) v->l->p = v;
    if(k >= v->key && ((v->r = emplace_inner(v->r, k, val)))) v->r->p = v;
    update(v);
    return balance(v);
  }
  // kの要素がすでにある場合何もしない
  node *emplace_unique_inner(node *v, Key k, Val val){
    if(!v) return new node(k, val);
    push_down(v);
    if(k < v->key  && ((v->l = emplace_unique_inner(v->l, k, val)))) v->l->p = v;
    if(k > v->key && ((v->r = emplace_unique_inner(v->r, k, val)))) v->r->p = v;
    update(v);
    return balance(v);
  }
  // kの要素がすでにある場合そのうちのいずれかと置き換える
  node *emplace_replace_inner(node *v, Key k, Val val){
    if(!v) return new node(k, val);
    push_down(v);
    if(k < v->key  && ((v->l = emplace_replace_inner(v->l, k, val)))) v->l->p = v;
    if(k > v->key && ((v->r = emplace_replace_inner(v->r, k, val)))) v->r->p = v;
    if(k == v->key) v->val = val;
    update(v);
    return balance(v);
  }
  // kを消す(複数ある場合は最も早く追加した要素を消す)
  node *erase_inner(node *v, Key k){
    if(!v) return nullptr;
    push_down(v);
    if(k < v->key && ((v->l = erase_inner(v->l, k)))) v->l->p = v;
    if(k > v->key && ((v->r = erase_inner(v->r, k)))) v->r->p = v;
    if(k == v->key){
      if(v->l && v->l->rmost->key == k){
        if((v->l = erase_inner(v->l, k))) v->l->p = v;
        update(v);
        return balance(v);
      }
      if(v->r){
        v->r = cut_leftmost(v->r);
        if((tmp_node->l = v->l)) tmp_node->l->p = tmp_node;
        if((tmp_node->r = v->r)) tmp_node->r->p = tmp_node;
        update(tmp_node);
        return balance(tmp_node);
      }
      return v->l;
    }
    update(v);
    return balance(v);
  }
  // vを消す
  node *erase_node_inner(node *v){
    assert(v);
    node *p = v->p, *u = v;
    if(v->r){
      v->r = cut_leftmost(v->r);
      if((tmp_node->l = v->l)) tmp_node->l->p = tmp_node;
      if((tmp_node->r = v->r)) tmp_node->r->p = tmp_node;
      update(tmp_node);
      v = balance(tmp_node);
    }else{
      v = v->l;
    }
    while(p){
      if(p->l == u && ((p->l = v))) v->p = p;
      if(p->r == u && ((p->r = v))) v->p = p;
      update(p);
      u = p;
      p = u->p;
      v = balance(u);
    }
    return v;
  }
  void update_inner(node *v, Key l, Key r, Lazy x){
    if(!v || v->rmost->key < l || r <= v->lmost->key) return;
    if(l <= v->lmost->key && v->rmost->key < r){
      propagate(v, x);
      return;
    }
    push_down(v);
    update_inner(v->l, l, r, x);
    update_inner(v->r, l, r, x);
    if(l <= v->key && v->key < r) v->val = apply(v->val, x, 0, 1);
    update(v);
  }
  Val query_inner(node *v, Key l, Key r){
    if(!v || v->rmost->key < l || r <= v->lmost->key) return id();
    if(l <= v->lmost->key && v->rmost->key < r) return v->sum;
    push_down(v);
    if(r <= v->key) return query_inner(v->l, l, r);
    if(v->key < l) return query_inner(v->r, l, r);
    return merge(query_inner(v->l, l, r), merge(v->val, query_inner(v->r, l, r)));
  }
  template<typename F>
  node *bisect_from_left(node *v, F &f, Key l, Val &ok){
    if(!v || v->rmost->key < l) return nullptr;
    push_down(v);
    Val m = merge(ok, v->sum);
    if(l <= v->lmost->key && !f(m)){
      ok = m;
      return nullptr;
    }
    node *x = bisect_from_left(v->l, f, l, ok);
    if(x) return x;
    if(l <= v->key){
      ok = merge(ok, v->val);
      if(f(ok)) return v;
    }
    return bisect_from_left(v->r, f, l, ok);
  }
  template<typename F>
  node *bisect_from_right(node *v, F &f, Key r, Val &ok){
    if(!v || v->lmost->key > r) return nullptr;
    push_down(v);
    Val m = merge(v->sum, ok);
    if(r >= v->rmost->key && !f(m)){
      ok = m;
      return nullptr;
    }
    node *x = bisect_from_right(v->r, f, r, ok);
    if(x) return x;
    if(v->key <= r){
      ok = merge(v->val, ok);
      if(f(ok)) return v;
    }
    return bisect_from_right(v->l, f, r, ok);
  }
  void to_list_inner(node *v, std::vector<std::pair<Key, Val>> &res){
    push_down(v);
    if(v->l) to_list_inner(v->l, res);
    res.push_back({v->key, v->val});
    if(v->r) to_list_inner(v->r, res);
  }
public:
  lazy_multimap_monoid(): root(nullptr){}
  lazy_multimap_monoid(std::vector<std::pair<Key, Val>> v){
    std::sort(v.begin(), v.end());
    init_sorted(v);
  }
  // 1
  // すでにソート済み
  void init_sorted(const std::vector<std::pair<Key, Val>> &v){
    if(v.empty()){
      root = nullptr;
      return;
    }
    int n = v.size();
    std::vector<node*> nodes(n);
    for(int i = 0; i < n; i++) nodes[i] = (new node(v[i].first, v[i].second));
    root = build(nodes, 0, n);
  }
  // 要素数
  int size(){
    return size_sum(root);
  }
  // 種類数
  int size_unique(){
    return size_unique(root);
  }
  bool empty(){
    return size_sum(root) == 0;
  }
  void emplace(Key k, Val val){
    root = emplace_inner(root, k, val);
    root->p = nullptr;
  }
  void emplace_unique(Key k, Val val){
    root = emplace_unique_inner(root, k, val);
    root->p = nullptr;
  }
  void emplace_replace(Key k, Val val){
    root = emplace_replace_inner(root, k, val);
    root->p = nullptr;
  }
  // 同じキーを持つ要素が複数ある場合, 最後に追加された要素を消す
  void erase(Key k){
    if((root = erase_inner(root, k))) root->p = nullptr;
  }
  void erase_node(node *v){
    if((root = erase_node_inner(v))) root->p = nullptr;
  }
  void clear(){
    root = nullptr;
  }
  bool contain(Key k){
    return find_key(root, k);
  }
  node *find(Key k){
    return find_key(root, k);
  }
  node *find_next(node *v){
    if(!v) return nullptr;
    return find_next_inner(v);
  }
  node *find_prev(node *v){
    if(!v) return nullptr;
    return find_prev_inner(v);
  }
  node *min(){
    return leftmost(root);
  }
  node *max(){
    return rightmost(root);
  }
  int low_count_unique(Key k){
    return low_count_unique_inner(root, k);
  }
  int low_count_sum(Key k){
    return low_count_sum_inner(root, k);
  }
  int low_count_node(node *v){
    return low_count_node_inner(v);
  }
  node *kth_smallest_unique(int k){
    if(size_unique() <= k) return nullptr;
    return find_kth_unique(root, k);
  }
  node *kth_smallest_sum(int k){
    if(size() <= k) return nullptr;
    return find_kth_sum(root, k);
  }
  node *lower_bound(Key k){
    return lower_bound_inner(root, k);
  }
  // k以下の最大要素
  node *lower_bound_rev(Key k){
    return lower_bound_rev_inner(root, k);
  }
  void update(Key l, Key r, Lazy x){
    update_inner(root, l, r, x);
  }
  void update_all(Lazy x){
    if(root) propagate(root, x);
  }
  Val query(Key l, Key r){
    return query_inner(root, l, r);
  }
  Val query_all(){
    if(!root) return id();
    return root->sum;
  }
  template<typename F>
  std::pair<node*, Val> bisect_from_left(Key l, F f){
    Val ok = id();
    node* ret = bisect_from_left(root, f, l, ok);
    return {ret, ok};
  }
  template<typename F>
  std::pair<node*, Val> bisect_from_right(Key r, F f){
    Val ok = id();
    node* ret = bisect_from_right(root, f, r, ok);
    return {ret, ok};
  }
  std::vector<std::pair<Key, Val>> to_list(){
    std::vector<std::pair<Key, Val>> res;
    if(root) to_list_inner(root, res);
    return res;
  }
};
#endif
