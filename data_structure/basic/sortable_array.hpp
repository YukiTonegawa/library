#ifndef _SORTABLE_ARRAY_H_
#define _SORTABLE_ARRAY_H_
#include <vector>
#include <iostream>
#include <cassert>
#include <numeric>
#include <algorithm>

template<typename Key, typename Val>
struct sortable_array{
private:
  struct node{
    Key k, kmin, kmax;
    Val val;
    bool flip;
    int sz_all, sz_inner;
    unsigned int p;
    node *inner_l, *inner_r, *outer_l, *outer_r;
    node(Key k, Val val, unsigned int p): k(k), kmin(k), kmax(k), val(val), flip(false),
    sz_all(1), sz_inner(1), p(p), inner_l(nullptr), inner_r(nullptr), outer_l(nullptr), outer_r(nullptr){}
  };
  using pn = std::pair<node*, node*>;
  using pt = std::tuple<node*, node*, node*>;
  static unsigned int xor128(){
    static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 86675123;
    unsigned int t = (x ^ (x << 11));
    x = y, y = z, z = w;
    return (w = (w ^ (w >> 19)) ^ (t ^ ( t >> 8)));
  }
  void update(node *v){
    v->sz_inner = 1;
    v->kmin = v->kmax = v->k;
    if(v->inner_l){
      v->sz_inner += v->inner_l->sz_inner;
      v->kmin = std::min(v->kmin, v->inner_l->kmin);
      v->kmax = std::max(v->kmax, v->inner_l->kmax);
    }
    if(v->inner_r){
      v->sz_inner += v->inner_r->sz_inner;
      v->kmin = std::min(v->kmin, v->inner_r->kmin);
      v->kmax = std::max(v->kmax, v->inner_r->kmax);
    }
    v->sz_all = v->sz_inner + (v->outer_l ? v->outer_l->sz_all : 0) + (v->outer_r ? v->outer_r->sz_all : 0);
  }
  int size(node *v){
    return v ? v->sz_all : 0;
  }
  node *make_node(Key k, Val val){
    return new node(k, val, xor128());
  }
  void flip(node *v){
    std::swap(v->inner_l, v->inner_r);
    std::swap(v->outer_l, v->outer_r);
    v->flip ^= 1;
  }
  void push_down(node *v){
    if(!v->flip) return;
    if(v->inner_l) flip(v->inner_l);
    if(v->inner_r) flip(v->inner_r);
    if(v->outer_l) flip(v->outer_l);
    if(v->outer_r) flip(v->outer_r);
    v->flip = false;
  }
  node *insert(node *v, node *n, int k){
    if(!v) return n;
    push_down(v);
    int szl = size(v->outer_l);
    if(v->p < n->p){
      pn s = split_outer(v, k);
      n->outer_l = s.first;
      n->outer_r = s.second;
      update(n);
      return n;
    }else if(szl <= k && k < szl + v->sz_inner){
      pn s = split_outer(v, k);
      s.first = insert(s.first, n, k);
      update(s.first);
      return merge_outer(s.first, s.second);
    }else if(k < szl){
      v->outer_l = insert(v->outer_l, n, k);
    }else{
      v->outer_r = insert(v->outer_r, n, k - szl - v->sz_inner);
    }
    update(v);
    return v;
  }
  node *merge_inner(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    push_down(a), push_down(b);
    if(a->p > b->p){
      a->inner_r = merge_inner(a->inner_r, b);
      update(a);
      return a;
    }else{
      b->inner_l = merge_inner(a, b->inner_l);
      update(b);
      return b;
    }
  }
  node *merge_outer(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    push_down(a), push_down(b);
    if(a->p > b->p){
      a->outer_r = merge_outer(a->outer_r, b);
      update(a);
      return a;
    }else{
      b->outer_l = merge_outer(a, b->outer_l);
      update(b);
      return b;
    }
  }
  pt cut_inner(node *v, int k){
    if(!v) return {nullptr, nullptr, nullptr};
    push_down(v);
    int szl = size(v->inner_l);
    if(k <= szl){
      if(k == szl){
        pt res = {v->inner_l, v, v->inner_r};
        v->inner_l = v->inner_r = nullptr;
        update(v);
        return res;
      }
      auto [a, b, c] = cut_inner(v->inner_l, k);
      v->inner_l = c;
      update(v);
      return {a, b, v};
    }else{
      auto [a, b, c] = cut_inner(v->inner_r, k - szl - 1);
      v->inner_r = a;
      update(v);
      return {v, b, c};
    }
  }
  pt cut_outer(node *v, int k){
    if(!v) return {nullptr, nullptr, nullptr};
    push_down(v);
    int szl = size(v->outer_l);
    int szr = szl + v->sz_inner;
    if(k < szl){
      auto [a, b, c] = cut_outer(v->outer_l, k);
      v->outer_l = c;
      update(v);
      return {a, b, v};
    }else if(szr <= k){
      auto [a, b, c] = cut_outer(v->outer_r, k - szr);
      v->outer_r = a;
      update(v);
      return {v, b, c};
    }else{
      node *tmp_l = v->outer_l, *tmp_r = v->outer_r;
      v->outer_l = v->outer_r = nullptr;
      auto [a, b, c] = cut_inner(v, k - szl);
      a = merge_outer(tmp_l, a);
      c = merge_outer(c, tmp_r);
      return {a, b, c};
    }
  }
  pn split_inner(node *v, int k){
    if(!v) return {nullptr, nullptr};
    push_down(v);
    int szl = size(v->inner_l);
    if(k <= szl){
      pn s = split_inner(v->inner_l, k);
      v->inner_l = s.second;
      update(v);
      return {s.first, v};
    }else{
      pn s = split_inner(v->inner_r, k - szl - 1);
      v->inner_r = s.first;
      update(v);
      return {v, s.second};
    }
  }
  pn split_outer(node *v, int k){
    if(!v) return {nullptr, nullptr};
    push_down(v);
    int szl = size(v->outer_l);
    int szr = szl + v->sz_inner;
    if(k < szl){
      pn s = split_outer(v->outer_l, k);
      v->outer_l = s.second;
      update(v);
      return {s.first, v};
    }else if(szr <= k){
      pn s = split_outer(v->outer_r, k - szr);
      v->outer_r = s.first;
      update(v);
      return {v, s.second};
    }else{
      node *tmp_l = v->outer_l, *tmp_r = v->outer_r;
      v->outer_l = v->outer_r = nullptr;
      pn s = split_inner(v, k - szl);
      s.first = merge_outer(tmp_l, s.first);
      s.second = merge_outer(s.second, tmp_r);
      return {s.first, s.second};
    }
  }
  pt split_range_outer(node *v, int l, int r){
    auto [a, b] = split_outer(v, l);
    auto [bb, c] = split_outer(b, r - l);
    return {a, bb, c};
  }
  pn split_key(node *v, Key k){
    if(!v) return {nullptr, nullptr};
    if(k < v->kmin) return {nullptr, v};
    else if(v->kmax <= k) return {v, nullptr};
    push_down(v);
    if(k < v->k){
      pn s = split_key(v->inner_l, k);
      v->inner_l = s.second;
      update(v);
      return {s.first, v};
    }else{
      pn s = split_key(v->inner_r, k);
      v->inner_r = s.first;
      update(v);
      return {v, s.second};
    }
  }
  node *merge_compress(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    push_down(a), push_down(b);
    if(a->p < b->p) std::swap(a, b);
    if(a->k <= b->kmin){
      a->inner_r = merge_compress(a->inner_r, b);
    }else if(b->kmax <= a->k){
      a->inner_l = merge_compress(a->inner_l, b);
    }else{
      auto [bl, br] = split_key(b, a->k);
      a->inner_l = merge_compress(a->inner_l, bl);
      a->inner_r = merge_compress(a->inner_r, br);
    }
    update(a);
    return a;
  }
  node *sort_inner(node *v){
    if(!v) return nullptr;
    node *tmp_l = v->outer_l, *tmp_r = v->outer_r;
    push_down(v);
    v->outer_l = v->outer_r = nullptr;
    node *res = sort_inner(tmp_l);
    if((v->inner_l && v->k < v->inner_l->kmax) || (v->inner_r && v->inner_r->kmin < v->k)) flip(v);
    res = merge_compress(res, v);
    res = merge_compress(res, sort_inner(tmp_r));
    return res;
  }
  void enumerate_inner(node *v, std::vector<std::pair<Key, Val>> &res){
    if(!v) return;
    push_down(v);
    if(v->inner_l) enumerate_inner(v->inner_l, res);
    res.push_back({v->k, v->val});
    if(v->inner_r) enumerate_inner(v->inner_r, res);
  }
  void enumerate_outer(node *v, std::vector<std::pair<Key, Val>> &res){
    if(!v) return;
    push_down(v);
    if(v->outer_l) enumerate_outer(v->outer_l, res);
    enumerate_inner(v, res);
    if(v->outer_r) enumerate_outer(v->outer_r, res);
  }
  void p_satisfy(node *v){
    if(!v->outer_l){
      if(!v->outer_r || v->p > v->outer_r->p) return;
      std::swap(v->p, v->outer_r->p);
      p_satisfy(v->outer_r);
    }else if(!v->outer_r){
      if(v->p > v->outer_l->p) return;
      std::swap(v->p, v->outer_l->p);
      p_satisfy(v->outer_l);
    }else{
      if(v->outer_l->p > v->outer_r->p){
        if(v->p > v->outer_l->p) return;
        std::swap(v->p, v->outer_l->p);
        p_satisfy(v->outer_l);
      }else{
        if(v->p > v->outer_r->p) return;
        std::swap(v->p, v->outer_r->p);
        p_satisfy(v->outer_r);
      }
    }
  }
  node *build(int l, int r, std::vector<std::pair<Key, Val>> &v){
    int mid = (l + r) / 2;
    node *u = make_node(v[mid].first, v[mid].second);
    if(l < mid){
      u->outer_l = build(l, mid, v);
    }else u->outer_l = nullptr;
    if(mid + 1 < r){
      u->outer_r = build(mid + 1, r, v);
    }else u->outer_r = nullptr;
    p_satisfy(u);
    update(u);
    return u;
  }
  node *root;
public:
  sortable_array(): root(nullptr){}
  sortable_array(std::vector<std::pair<Key, Val>> &v): root(v.empty() ? nullptr : build(0, v.size(), v)){}

  // k番目に{key, x}を追加
  void insert(int k, Key key, Val x){
    root = insert(root, make_node(key, x), k);
  }
  // k番目の要素を削除
  void erase(int k){
    auto [a, b, c] = cut_outer(root, k);
    root = merge_outer(a, c);
  }
  // k番目の要素の　{キー, 値}を{key, x}に変更
  void set(int k, Key key, Val x){
    auto [a, b, c] = cut_outer(root, k);
    b->k = key;
    b->val = x;
    update(b);
    root = merge_outer(a, merge_outer(b, c));
  }
  // 値
  Val get(int k){
    auto [a, b, c] = cut_outer(root, k);
    Val ret = b->val;
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // {キー, 値}
  std::pair<Key, Val> get2(int k){
    auto [a, b, c] = cut_outer(root, k);
    std::pair<Key, Val> ret = {b->k, b->val};
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // [l, r)のキーの最小値
  Key min_key(int l, int r){
    auto [a, b, c] = split_range_outer(root, l, r);
    assert(b);
    Key ret = b->kmin;
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // [l, r)のキーの最大値
  Key max_key(int l, int r){
    auto [a, b, c] = split_range_outer(root, l, r);
    assert(b);
    Key ret = b->kmax;
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // [l, r)を反転
  void flip(int l, int r){
    if(r == l) return;
    auto [a, b, c] = split_range_outer(root, l, r);
    flip(b);
    root = merge_outer(a, merge_outer(b, c));
  }
  // [l, r)を昇順ソート
  void sort_ascending(int l, int r){
    if(r == l) return;
    auto [a, b, c] = split_range_outer(root, l, r);
    b = sort_inner(b);
    root = merge_outer(a, merge_outer(b, c));
  }
  // [l, r)を降順ソート
  void sort_descending(int l, int r){
    if(r == l) return;
    auto [a, b, c] = split_range_outer(root, l, r);
    b = sort_inner(b);
    flip(b);
    root = merge_outer(a, merge_outer(b, c));
  }
  std::vector<std::pair<Key, Val>> to_vector(){
    std::vector<std::pair<Key, Val>> res;
    enumerate_outer(root, res);
    return res;
  }
};


// 区間がソート済か判定できるバージョン
template<typename Key, typename Val>
struct sortable_array2{
private:
  struct node{
    Key k, kmin, kmax, kminall, kmaxall;
    Val val;
    bool flip, asc, des;
    int sz_all, sz_inner;
    unsigned int p;
    node *inner_l, *inner_r, *outer_l, *outer_r;
    node(Key k, Val val, unsigned int p): k(k), kmin(k), kmax(k), kminall(k), kmaxall(k), val(val), flip(false), asc(true), des(true), 
    sz_all(1), sz_inner(1), p(p), inner_l(nullptr), inner_r(nullptr), outer_l(nullptr), outer_r(nullptr){}
  };
  using pn = std::pair<node*, node*>;
  using pt = std::tuple<node*, node*, node*>;
  static unsigned int xor128(){
    static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 86675123;
    unsigned int t = (x ^ (x << 11));
    x = y, y = z, z = w;
    return (w = (w ^ (w >> 19)) ^ (t ^ ( t >> 8)));
  }
  void update(node *v){
    v->sz_inner = 1;
    v->kmin = v->kmax = v->k;
    v->asc = v->des = true;
    if(v->inner_l){
      v->sz_inner += v->inner_l->sz_inner;
      if(v->inner_l->kmax > v->k) v->asc = false;
      if(v->inner_l->kmin < v->k) v->des = false;
      v->kmin = std::min(v->kmin, v->inner_l->kmin);
      v->kmax = std::max(v->kmax, v->inner_l->kmax);
    }
    if(v->inner_r){
      v->sz_inner += v->inner_r->sz_inner;
      if(v->inner_r->kmin < v->k) v->asc = false;
      if(v->inner_r->kmax > v->k) v->des = false;
      v->kmin = std::min(v->kmin, v->inner_r->kmin);
      v->kmax = std::max(v->kmax, v->inner_r->kmax);
    }
    v->kminall = v->kmin, v->kmaxall = v->kmax;
    if(v->outer_l){
      if(!v->outer_l->asc || v->outer_l->kmaxall > v->kminall) v->asc = false;
      if(!v->outer_l->des || v->outer_l->kminall < v->kmaxall) v->des = false;
      v->kminall = std::min(v->kminall, v->outer_l->kminall);
      v->kmaxall = std::max(v->kmaxall, v->outer_l->kmaxall);
    }
    if(v->outer_r){
      if(!v->outer_r->asc || v->outer_r->kminall < v->kmaxall) v->asc = false;
      if(!v->outer_r->des || v->outer_r->kmaxall > v->kminall) v->des = false;
      v->kminall = std::min(v->kminall, v->outer_r->kminall);
      v->kmaxall = std::max(v->kmaxall, v->outer_r->kmaxall);
    }
    v->sz_all = v->sz_inner + (v->outer_l ? v->outer_l->sz_all : 0) + (v->outer_r ? v->outer_r->sz_all : 0);
  }
  int size(node *v){
    return v ? v->sz_all : 0;
  }
  node *make_node(Key k, Val val){
    return new node(k, val, xor128());
  }
  void flip(node *v){
    std::swap(v->inner_l, v->inner_r);
    std::swap(v->outer_l, v->outer_r);
    std::swap(v->asc, v->des);
    v->flip ^= 1;
  }
  void push_down(node *v){
    if(!v->flip) return;
    if(v->inner_l) flip(v->inner_l);
    if(v->inner_r) flip(v->inner_r);
    if(v->outer_l) flip(v->outer_l);
    if(v->outer_r) flip(v->outer_r);
    v->flip = false;
  }
  node *insert(node *v, node *n, int k){
    if(!v) return n;
    push_down(v);
    int szl = size(v->outer_l);
    if(v->p < n->p){
      pn s = split_outer(v, k);
      n->outer_l = s.first;
      n->outer_r = s.second;
      update(n);
      return n;
    }else if(szl <= k && k < szl + v->sz_inner){
      pn s = split_outer(v, k);
      s.first = insert(s.first, n, k);
      update(s.first);
      return merge_outer(s.first, s.second);
    }else if(k < szl){
      v->outer_l = insert(v->outer_l, n, k);
    }else{
      v->outer_r = insert(v->outer_r, n, k - szl - v->sz_inner);
    }
    update(v);
    return v;
  }
  node *merge_inner(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    push_down(a), push_down(b);
    if(a->p > b->p){
      a->inner_r = merge_inner(a->inner_r, b);
      update(a);
      return a;
    }else{
      b->inner_l = merge_inner(a, b->inner_l);
      update(b);
      return b;
    }
  }
  node *merge_outer(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    push_down(a), push_down(b);
    if(a->p > b->p){
      a->outer_r = merge_outer(a->outer_r, b);
      update(a);
      return a;
    }else{
      b->outer_l = merge_outer(a, b->outer_l);
      update(b);
      return b;
    }
  }
  pt cut_inner(node *v, int k){
    if(!v) return {nullptr, nullptr, nullptr};
    push_down(v);
    int szl = size(v->inner_l);
    if(k <= szl){
      if(k == szl){
        pt res = {v->inner_l, v, v->inner_r};
        v->inner_l = v->inner_r = nullptr;
        update(v);
        return res;
      }
      auto [a, b, c] = cut_inner(v->inner_l, k);
      v->inner_l = c;
      update(v);
      return {a, b, v};
    }else{
      auto [a, b, c] = cut_inner(v->inner_r, k - szl - 1);
      v->inner_r = a;
      update(v);
      return {v, b, c};
    }
  }
  pt cut_outer(node *v, int k){
    if(!v) return {nullptr, nullptr, nullptr};
    push_down(v);
    int szl = size(v->outer_l);
    int szr = szl + v->sz_inner;
    if(k < szl){
      auto [a, b, c] = cut_outer(v->outer_l, k);
      v->outer_l = c;
      update(v);
      return {a, b, v};
    }else if(szr <= k){
      auto [a, b, c] = cut_outer(v->outer_r, k - szr);
      v->outer_r = a;
      update(v);
      return {v, b, c};
    }else{
      node *tmp_l = v->outer_l, *tmp_r = v->outer_r;
      v->outer_l = v->outer_r = nullptr;
      auto [a, b, c] = cut_inner(v, k - szl);
      a = merge_outer(tmp_l, a);
      c = merge_outer(c, tmp_r);
      return {a, b, c};
    }
  }
  pn split_inner(node *v, int k){
    if(!v) return {nullptr, nullptr};
    push_down(v);
    int szl = size(v->inner_l);
    if(k <= szl){
      pn s = split_inner(v->inner_l, k);
      v->inner_l = s.second;
      update(v);
      return {s.first, v};
    }else{
      pn s = split_inner(v->inner_r, k - szl - 1);
      v->inner_r = s.first;
      update(v);
      return {v, s.second};
    }
  }
  pn split_outer(node *v, int k){
    if(!v) return {nullptr, nullptr};
    push_down(v);
    int szl = size(v->outer_l);
    int szr = szl + v->sz_inner;
    if(k < szl){
      pn s = split_outer(v->outer_l, k);
      v->outer_l = s.second;
      update(v);
      return {s.first, v};
    }else if(szr <= k){
      pn s = split_outer(v->outer_r, k - szr);
      v->outer_r = s.first;
      update(v);
      return {v, s.second};
    }else{
      node *tmp_l = v->outer_l, *tmp_r = v->outer_r;
      v->outer_l = v->outer_r = nullptr;
      pn s = split_inner(v, k - szl);
      s.first = merge_outer(tmp_l, s.first);
      s.second = merge_outer(s.second, tmp_r);
      return {s.first, s.second};
    }
  }
  pt split_range_outer(node *v, int l, int r){
    auto [a, b] = split_outer(v, l);
    auto [bb, c] = split_outer(b, r - l);
    return {a, bb, c};
  }
  pn split_key(node *v, Key k){
    if(!v) return {nullptr, nullptr};
    if(k < v->kmin) return {nullptr, v};
    else if(v->kmax <= k) return {v, nullptr};
    push_down(v);
    if(k < v->k){
      pn s = split_key(v->inner_l, k);
      v->inner_l = s.second;
      update(v);
      return {s.first, v};
    }else{
      pn s = split_key(v->inner_r, k);
      v->inner_r = s.first;
      update(v);
      return {v, s.second};
    }
  }
  node *merge_compress(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    push_down(a), push_down(b);
    if(a->p < b->p) std::swap(a, b);
    if(a->k <= b->kmin){
      a->inner_r = merge_compress(a->inner_r, b);
    }else if(b->kmax <= a->k){
      a->inner_l = merge_compress(a->inner_l, b);
    }else{
      auto [bl, br] = split_key(b, a->k);
      a->inner_l = merge_compress(a->inner_l, bl);
      a->inner_r = merge_compress(a->inner_r, br);
    }
    update(a);
    return a;
  }
  node *sort_inner(node *v){
    if(!v) return nullptr;
    node *tmp_l = v->outer_l, *tmp_r = v->outer_r;
    push_down(v);
    v->outer_l = v->outer_r = nullptr;
    node *res = sort_inner(tmp_l);
    if((v->inner_l && v->k < v->inner_l->kmax) || (v->inner_r && v->inner_r->kmin < v->k)) flip(v);
    res = merge_compress(res, v);
    res = merge_compress(res, sort_inner(tmp_r));
    return res;
  }
  void enumerate_inner(node *v, std::vector<std::pair<Key, Val>> &res){
    if(!v) return;
    push_down(v);
    if(v->inner_l) enumerate_inner(v->inner_l, res);
    res.push_back({v->k, v->val});
    if(v->inner_r) enumerate_inner(v->inner_r, res);
  }
  void enumerate_outer(node *v, std::vector<std::pair<Key, Val>> &res){
    if(!v) return;
    push_down(v);
    if(v->outer_l) enumerate_outer(v->outer_l, res);
    enumerate_inner(v, res);
    if(v->outer_r) enumerate_outer(v->outer_r, res);
  }
  void p_satisfy(node *v){
    if(!v->outer_l){
      if(!v->outer_r || v->p > v->outer_r->p) return;
      std::swap(v->p, v->outer_r->p);
      p_satisfy(v->outer_r);
    }else if(!v->outer_r){
      if(v->p > v->outer_l->p) return;
      std::swap(v->p, v->outer_l->p);
      p_satisfy(v->outer_l);
    }else{
      if(v->outer_l->p > v->outer_r->p){
        if(v->p > v->outer_l->p) return;
        std::swap(v->p, v->outer_l->p);
        p_satisfy(v->outer_l);
      }else{
        if(v->p > v->outer_r->p) return;
        std::swap(v->p, v->outer_r->p);
        p_satisfy(v->outer_r);
      }
    }
  }
  node *build(int l, int r, std::vector<std::pair<Key, Val>> &v){
    int mid = (l + r) / 2;
    node *u = make_node(v[mid].first, v[mid].second);
    if(l < mid){
      u->outer_l = build(l, mid, v);
    }else u->outer_l = nullptr;
    if(mid + 1 < r){
      u->outer_r = build(mid + 1, r, v);
    }else u->outer_r = nullptr;
    p_satisfy(u);
    update(u);
    return u;
  }
  // lkmaxから初めて何連続で昇順に並んでいるか
  int find_right_ascending(node *v, Key lkmax){
    if(v->asc){
      if(lkmax <= v->kminall) return v->sz_all;
      return 0;
    }
    push_down(v);
    int ret = 0;
    if(v->outer_l){
      if(!v->outer_l->asc) return ret + find_right_ascending(v->outer_l, lkmax);
      else if(lkmax > v->outer_l->kminall) return ret;
      else{
        ret += v->outer_l->sz_all;
        lkmax = v->outer_l->kmaxall;
      }
    }
    if(v->inner_l){
      if(!v->inner_l->asc) return ret + find_right_ascending(v->inner_l, lkmax);
      else if(lkmax > v->inner_l->kminall) return ret;
      else{
        ret += v->inner_l->sz_all;
        lkmax = v->inner_l->kmaxall;
      }
    }
    if(lkmax > v->k) return ret;
    ret++;
    lkmax = v->k;
    if(v->inner_r){
      if(!v->inner_r->asc) return ret + find_right_ascending(v->inner_r, lkmax);
      else if(lkmax > v->inner_r->kminall) return ret;
      else{
        ret += v->inner_r->sz_all;
        lkmax = v->inner_r->kmaxall;
      }
    }
    if(v->outer_r){
      if(!v->outer_r->asc) return ret + find_right_ascending(v->outer_r, lkmax);
      else if(lkmax > v->outer_r->kminall) return ret;
      else{
        ret += v->outer_r->sz_all;
        lkmax = v->outer_r->kmaxall;
      }
    }
    return ret;
  }
   // lkminから初めて何連続で降順に並んでいるか
  int find_right_descending(node *v, Key lkmin){
    if(v->des){
      if(lkmin >= v->kmaxall) return v->sz_all;
      return 0;
    }
    push_down(v);
    int ret = 0;
    if(v->outer_l){
      if(!v->outer_l->des) return ret + find_right_descending(v->outer_l, lkmin);
      else if(lkmin < v->outer_l->kmaxall) return ret;
      else{
        ret += v->outer_l->sz_all;
        lkmin = v->outer_l->kminall;
      }
    }
    if(v->inner_l){
      if(!v->inner_l->des) return ret + find_right_descending(v->inner_l, lkmin);
      else if(lkmin < v->inner_l->kmaxall) return ret;
      else{
        ret += v->inner_l->sz_all;
        lkmin = v->inner_l->kminall;
      }
    }
    if(lkmin < v->k) return ret;
    ret++;
    lkmin = v->k;
    if(v->inner_r){
      if(!v->inner_r->des) return ret + find_right_descending(v->inner_r, lkmin);
      else if(lkmin < v->inner_r->kmaxall) return ret;
      else{
        ret += v->inner_r->sz_all;
        lkmin = v->inner_r->kminall;
      }
    }
    if(v->outer_r){
      if(!v->outer_r->des) return ret + find_right_descending(v->outer_r, lkmin);
      else if(lkmin < v->outer_r->kmaxall) return ret;
      else{
        ret += v->outer_r->sz_all;
        lkmin = v->outer_r->kminall;
      }
    }
    return ret;
  }
  // rkminから初めて何連続で昇順に並んでいるか
  int find_left_ascending(node *v, Key rkmin){
    if(v->asc){
      if(rkmin >= v->kmaxall) return v->sz_all;
      return 0;
    }
    push_down(v);
    int ret = 0;
    if(v->outer_r){
      if(!v->outer_r->asc) return ret + find_left_ascending(v->outer_r, rkmin);
      else if(rkmin < v->outer_r->kmaxall) return ret;
      else{
        ret += v->outer_r->sz_all;
        rkmin = v->outer_r->kminall;
      }
    }
    if(v->inner_r){
      if(!v->inner_r->asc) return ret + find_left_ascending(v->inner_r, rkmin);
      else if(rkmin < v->inner_r->kmaxall) return ret;
      else{
        ret += v->inner_r->sz_all;
        rkmin = v->inner_r->kminall;
      }
    }
    if(rkmin < v->k) return ret;
    ret++;
    rkmin = v->k;
    if(v->inner_l){
      if(!v->inner_l->asc) return ret + find_left_ascending(v->inner_l, rkmin);
      else if(rkmin < v->inner_l->kmaxall) return ret;
      else{
        ret += v->inner_l->sz_all;
        rkmin = v->inner_l->kminall;
      }
    }
    if(v->outer_l){
      if(!v->outer_l->asc) return ret + find_left_ascending(v->outer_l, rkmin);
      else if(rkmin < v->outer_l->kmaxall) return 0;
      else{
        ret += v->outer_l->sz_all;
        rkmin = v->outer_l->kminall;
      }
    }
    return ret;
  }
  // rkmaxから初めて何連続で降順に並んでいるか
  int find_left_descending(node *v, Key rkmax){
    if(v->des){
      if(rkmax <= v->kminall) return v->sz_all;
      return 0;
    }
    push_down(v);
    int ret = 0;
    if(v->outer_r){
      if(!v->outer_r->des) return ret + find_left_descending(v->outer_r, rkmax);
      else if(rkmax > v->outer_r->kminall) return ret;
      else{
        ret += v->outer_r->sz_all;
        rkmax = v->outer_r->kmaxall;
      }
    }
    if(v->inner_r){
      if(!v->inner_r->des) return ret + find_left_descending(v->inner_r, rkmax);
      else if(rkmax > v->inner_r->kminall) return ret;
      else{
        ret += v->inner_r->sz_all;
        rkmax = v->inner_r->kmaxall;
      }
    }
    if(rkmax > v->k) return ret;
    ret++;
    rkmax = v->k;
    if(v->inner_l){
      if(!v->inner_l->des) return ret + find_left_descending(v->inner_l, rkmax);
      else if(rkmax > v->inner_l->kminall) return ret;
      else{
        ret += v->inner_l->sz_all;
        rkmax = v->inner_l->kmaxall;
      }
    }
    if(v->outer_l){
      if(!v->outer_l->des) return ret + find_left_descending(v->outer_l, rkmax);
      else if(rkmax > v->outer_l->kminall) return 0;
      else{
        ret += v->outer_l->sz_all;
        rkmax = v->outer_l->kmaxall;
      }
    }
    return ret;
  }
  node *root;
public:
  sortable_array2(): root(nullptr){}
  sortable_array2(std::vector<std::pair<Key, Val>> &v): root(v.empty() ? nullptr : build(0, v.size(), v)){}

  // k番目に{key, x}を追加
  void insert(int k, Key key, Val x){
    root = insert(root, make_node(key, x), k);
  }
  // k番目の要素を削除
  void erase(int k){
    auto [a, b, c] = cut_outer(root, k);
    root = merge_outer(a, c);
  }
  // k番目の要素の　{キー, 値}を{key, x}に変更
  void set(int k, Key key, Val x){
    auto [a, b, c] = cut_outer(root, k);
    b->k = key;
    b->val = x;
    update(b);
    root = merge_outer(a, merge_outer(b, c));
  }
  // 値
  Val get(int k){
    auto [a, b, c] = cut_outer(root, k);
    Val ret = b->val;
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // {キー, 値}
  std::pair<Key, Val> get2(int k){
    auto [a, b, c] = cut_outer(root, k);
    std::pair<Key, Val> ret = {b->k, b->val};
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // [l, r)のキーの最小値
  Key min_key(int l, int r){
    auto [a, b, c] = split_range_outer(root, l, r);
    assert(b);
    Key ret = b->kmin;
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // [l, r)のキーの最大値
  Key max_key(int l, int r){
    auto [a, b, c] = split_range_outer(root, l, r);
    assert(b);
    Key ret = b->kmax;
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // [l, r)を反転
  void flip(int l, int r){
    if(r == l) return;
    auto [a, b, c] = split_range_outer(root, l, r);
    flip(b);
    root = merge_outer(a, merge_outer(b, c));
  }
  // [l, r)を昇順ソート
  void sort_ascending(int l, int r){
    if(r == l) return;
    auto [a, b, c] = split_range_outer(root, l, r);
    b = sort_inner(b);
    root = merge_outer(a, merge_outer(b, c));
  }
  // [l, r)を降順ソート
  void sort_descending(int l, int r){
    if(r == l) return;
    auto [a, b, c] = split_range_outer(root, l, r);
    b = sort_inner(b);
    flip(b);
    root = merge_outer(a, merge_outer(b, c));
  }
  // [l, r)が昇順ソート済か　(長さ1の場合昇順かつ降順)
  bool is_ascending(int l, int r){
    if(r == l) return false;
    auto [a, b, c] = split_range_outer(root, l, r);
    bool ret = b->asc;
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // [l, r)が降順ソート済か(長さ1の場合昇順かつ降順)
  bool is_descending(int l, int r){
    if(r == l) return false;
    auto [a, b, c] = split_range_outer(root, l, r);
    bool ret = b->des;
    root = merge_outer(a, merge_outer(b, c));
    return ret;
  }
  // lから始めて何連続で昇順に並んでいるか
  int find_right_ascending(int l){
    auto [a, b] = split_outer(root, l);
    int ret = find_right_ascending(b, b->kminall);
    root = merge_outer(a, b);
    return ret;
  }
  // lから始めて何連続で降順に並んでいるか
  int find_right_descending(int l){
    auto [a, b] = split_outer(root, l);
    int ret = find_right_descending(b, b->kmaxall);
    root = merge_outer(a, b);
    return ret;
  }
  // rから始めて何連続で昇順(左から見ての昇順)に並んでいるか(rを含む)
  int find_left_ascending(int r){
    auto [a, b] = split_outer(root, r + 1);
    int ret = find_left_ascending(b, b->kmaxall);
    root = merge_outer(a, b);
    return ret;
  }
  // rから始めて何連続で降順(左から見ての降順)に並んでいるか(rを含む)
  int find_left_descending(int r){
    auto [a, b] = split_outer(root, r + 1);
    int ret = find_left_descending(b, b->kminall);
    root = merge_outer(a, b);
    return ret;
  }
  std::vector<std::pair<Key, Val>> to_vector(){
    std::vector<std::pair<Key, Val>> res;
    enumerate_outer(root, res);
    return res;
  }
};
#endif