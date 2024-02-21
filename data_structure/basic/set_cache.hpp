#ifndef _SET_CACHE_H_
#define _SET_CACHE_H_
#include <vector>
#include <cassert>
#include <array>
#include <limits>
#include <algorithm>
#include <iostream>

template<typename Key>
struct set_avl_cache{
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static constexpr int limit_size_per_node = 32;
private:
  struct node{
    int h, sz, sz_sum;
    std::array<Key, limit_size_per_node> keys;
    node *l, *r;
    node(): h(1), l(nullptr), r(nullptr){}
    node(Key _key): h(1), sz(1), sz_sum(1), l(nullptr), r(nullptr){keys[0] = _key;}
    node(const std::vector<Key> &v, int l, int r): h(1), sz(r - l), sz_sum(sz), l(nullptr), r(nullptr){
      assert(sz < limit_size_per_node);
      for(int i = 0; i < sz; i++) keys[i] = v[l + i];
    }
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
    node *split_half(){
      assert(sz == limit_size_per_node);
      node *u = new node();
      sz = limit_size_per_node / 2;
      u->sz_sum = u->sz = limit_size_per_node - sz;
      for(int i = 0; i < u->sz; i++) u->keys[i] = keys[sz + i];
      return u;
    }
  };
  node *root, *tmp_node;
  int size(node *v){return v ? v->sz_sum : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz_sum = (v->l ? v->l->sz_sum : 0) + (v->r ? v->r->sz_sum : 0) + v->sz;
  }
  node *rotate_right(node *v){
    node *l = v->l;
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  node *rotate_left(node *v){
    node *r = v->r;
    v->r = r->l;
    r->l = v;
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
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  node *leftmost(node *v){
    while(v->l) v = v->l;
    return v;
  }
  node *rightmost(node *v){
    while(v->r) v = v->r;
    return v;
  }
  Key find_kth(node *v, int k){
    while(true){
      int szl = size(v->l);
      if(szl <= k){
        if(szl + v->sz > k) return v->keys[k - szl];
        k -= szl + v->sz;
        v = v->r;
      }else v = v->l;
    }
  }
  bool find_key(node *v, Key k){
    while(v){
      if(k < v->keys[0]) v = v->l;
      else if(k > v->keys[v->sz - 1]) v = v->r;
      else{
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        if(idx < v->sz && v->keys[idx] == k) return true;
        return false;
      }
    }
    return false;
  }
  node *cut_leftmost(node *v){
    if(v->l){
      v->l = cut_leftmost(v->l);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  node *cut_rightmost(node *v){
    if(v->r){
      v->r = cut_rightmost(v->r);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->l;
  }
  node *insert_leftmost(node *v, node *u){
    if(!v) return u;
    v->l = insert_leftmost(v->l, u);
    update(v);
    return balance(v);
  }
  node *insert_inner(node *v, Key k){
    if(!v) return new node(k);
    if(v->l && k < v->keys[0]){
      v->l = insert_inner(v->l, k);
    }else if(v->r && v->keys[v->sz - 1] < k){
      v->r = insert_inner(v->r, k);
    }else{
      int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
      if(idx < v->sz && v->keys[idx] == k) return v;
      for(int i = v->sz; i > idx; i--) v->keys[i] = v->keys[i - 1];
      v->keys[idx] = k;
      v->sz++;
      if(v->sz == limit_size_per_node){
        v->r = insert_leftmost(v->r, v->split_half());
      }
    }
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, Key k){
    if(!v) return nullptr;
    if(k < v->keys[0]){
      v->l = erase_inner(v->l, k);
    }else if(v->keys[v->sz - 1] < k){
      v->r = erase_inner(v->r, k);
    }else{
      for(int i = 0; i < v->sz; i++){
        if(v->keys[i] == k){
          for(int j = i + 1; j < v->sz; j++) v->keys[j - 1] = v->keys[j];
          v->sz--;
          break;
        }
      }
      if(!v->sz){
        if(v->r){
          v->r = cut_leftmost(v->r);
          tmp_node->l = v->l;
          tmp_node->r = v->r;
          update(tmp_node);
          return balance(tmp_node);
        }else{
          return v->l;
        }
      }
    }
    update(v);
    return balance(v);
  }
  int low_count_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(k < v->keys[0]) v = v->l;
      else if(k > v->keys[v->sz - 1]) res += szl + v->sz, v = v->r;
      else{
        res += szl;
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        return res + idx;
      }
    }
    return res;
  }
  Key lower_bound_inner(node *v, Key k){
    Key res = inf;
    while(v){
      if(k < v->keys[0]){
        res = v->keys[0];
        v = v->l;
      }else if(k > v->keys[v->sz - 1]){
        v = v->r;
      }else{
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        assert(idx < v->sz);
        return v->keys[idx];
      }
    }
    return res;
  }
  Key lower_bound_rev_inner(node *v, Key k){
    Key res = inf;
    while(v){
      if(k < v->keys[0]){
        v = v->l;
      }else if(k > v->keys[v->sz - 1]){
        res = v->keys[v->sz - 1];
        v = v->r;
      }else{
        int idx = std::upper_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        assert(idx != 0);
        return v->keys[idx - 1];
      }
    }
    return res;
  }
public:
  set_avl_cache(): root(nullptr){}
  set_avl_cache(const std::vector<Key> &v){
    init_sorted(v, false);
  }
  // すでにソート済み
  void init_sorted(std::vector<Key> v, bool sorted = true){
    if(v.empty()){
      root = nullptr;
      return;
    }
    if(!sorted) std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    int n = v.size();
    int size_per_node = limit_size_per_node / 2;
    int m = (n + size_per_node - 1) / size_per_node;
    std::vector<node*> nodes(m);
    for(int i = 0; i < m; i++){
      nodes[i] = new node(v, i * size_per_node, std::min((i + 1) * size_per_node, n));
    }
    root = build(nodes, 0, m);
  }
  int size(){
    return size(root);
  }
  void insert(Key key){
    root = insert_inner(root, key);
  }
  void erase(Key k){
    root = erase_inner(root, k);
  }
  bool find(Key k){
    return find_key(root, k);
  }
  Key min(){
    assert(size());
    return leftmost(root)->keys[0];
  }
  Key max(){
    assert(size());
    node *v = rightmost(root);
    return v->keys[v->sz - 1];
  }
  // k未満の値の数
  int low_count(Key k){
    return low_count_inner(root, k);
  }
  // k番目(0-indexed)に小さいキー
  Key kth_smallest(int k){
    if(size() <= k) return inf;
    return find_kth(root, k);
  }
  // k以上の最小要素
  Key lower_bound(Key k){
    return lower_bound_inner(root, k);
  }
  // k以下の最大要素
  Key lower_bound_rev(Key k){
    return lower_bound_rev_inner(root, k);
  }
};

template<typename Key>
struct multiset_avl_cache{
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static constexpr int limit_size_per_node = 32;
private:
  struct node{
    int h, sz, sz_sum;
    std::array<Key, limit_size_per_node> keys;
    node *l, *r;
    node(): h(1), l(nullptr), r(nullptr){}
    node(Key _key): h(1), sz(1), sz_sum(1), l(nullptr), r(nullptr){keys[0] = _key;}
    node(const std::vector<Key> &v, int l, int r): h(1), sz(r - l), sz_sum(sz), l(nullptr), r(nullptr){
      assert(sz < limit_size_per_node);
      for(int i = 0; i < sz; i++) keys[i] = v[l + i];
    }
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
    node *split_half(){
      assert(sz == limit_size_per_node);
      node *u = new node();
      sz = limit_size_per_node / 2;
      u->sz_sum = u->sz = limit_size_per_node - sz;
      for(int i = 0; i < u->sz; i++) u->keys[i] = keys[sz + i];
      return u;
    }
  };
  node *root, *tmp_node;
  int size(node *v){return v ? v->sz_sum : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz_sum = (v->l ? v->l->sz_sum : 0) + (v->r ? v->r->sz_sum : 0) + v->sz;
  }
  node *rotate_right(node *v){
    node *l = v->l;
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  node *rotate_left(node *v){
    node *r = v->r;
    v->r = r->l;
    r->l = v;
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
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  node *leftmost(node *v){
    while(v->l) v = v->l;
    return v;
  }
  node *rightmost(node *v){
    while(v->r) v = v->r;
    return v;
  }
  Key find_kth(node *v, int k){
    while(true){
      int szl = size(v->l);
      if(szl <= k){
        if(szl + v->sz > k) return v->keys[k - szl];
        k -= szl + v->sz;
        v = v->r;
      }else v = v->l;
    }
  }
  bool find_key(node *v, Key k){
    while(v){
      if(k < v->keys[0]) v = v->l;
      else if(k > v->keys[v->sz - 1]) v = v->r;
      else{
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        if(idx < v->sz && v->keys[idx] == k) return true;
        return false;
      }
    }
    return false;
  }
  node *cut_leftmost(node *v){
    if(v->l){
      v->l = cut_leftmost(v->l);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  node *cut_rightmost(node *v){
    if(v->r){
      v->r = cut_rightmost(v->r);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->l;
  }
  node *insert_leftmost(node *v, node *u){
    if(!v) return u;
    v->l = insert_leftmost(v->l, u);
    update(v);
    return balance(v);
  }
  node *insert_inner(node *v, Key k){
    if(!v) return new node(k);
    if(v->l && k < v->keys[0]){
      v->l = insert_inner(v->l, k);
    }else if(v->r && v->keys[v->sz - 1] < k){
      v->r = insert_inner(v->r, k);
    }else{
      int idx = v->sz;
      while(idx){
        if(v->keys[idx - 1] <= k) break;
        v->keys[idx] = v->keys[idx - 1];
        idx--;
      }
      v->keys[idx] = k;
      v->sz++;
      if(v->sz == limit_size_per_node){
        v->r = insert_leftmost(v->r, v->split_half());
      }
    }
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, Key k){
    if(!v) return nullptr;
    if(k < v->keys[0]){
      v->l = erase_inner(v->l, k);
    }else if(v->keys[v->sz - 1] < k){
      v->r = erase_inner(v->r, k);
    }else{
      for(int i = 0; i < v->sz; i++){
        if(v->keys[i] == k){
          for(int j = i + 1; j < v->sz; j++) v->keys[j - 1] = v->keys[j];
          v->sz--;
          break;
        }
      }
      if(!v->sz){
        if(v->r){
          v->r = cut_leftmost(v->r);
          tmp_node->l = v->l;
          tmp_node->r = v->r;
          update(tmp_node);
          return balance(tmp_node);
        }else{
          return v->l;
        }
      }
    }
    update(v);
    return balance(v);
  }
  int low_count_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(k <= v->keys[0]) v = v->l;
      else if(k > v->keys[v->sz - 1]) res += szl + v->sz, v = v->r;
      else{
        res += szl;
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        return res + idx;
      }
    }
    return res;
  }
  int low_count_eq_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(k < v->keys[0]) v = v->l;
      else if(k >= v->keys[v->sz - 1]) res += szl + v->sz, v = v->r;
      else{
        res += szl;
        int idx = std::upper_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        return res + idx;
      }
    }
    return res;
  }
  Key lower_bound_inner(node *v, Key k){
    Key res = inf;
    while(v){
      if(k < v->keys[0]){
        res = v->keys[0];
        v = v->l;
      }else if(k > v->keys[v->sz - 1]){
        v = v->r;
      }else{
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        assert(idx < v->sz);
        return v->keys[idx];
      }
    }
    return res;
  }
  Key lower_bound_rev_inner(node *v, Key k){
    Key res = inf;
    while(v){
      if(k < v->keys[0]){
        v = v->l;
      }else if(k > v->keys[v->sz - 1]){
        res = v->keys[v->sz - 1];
        v = v->r;
      }else{
        int idx = std::upper_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        assert(idx != 0);
        return v->keys[idx - 1];
      }
    }
    return res;
  }
public:
  multiset_avl_cache(): root(nullptr){}
  multiset_avl_cache(std::vector<Key> v){
    std::sort(v.begin(), v.end());
    init_sorted(v);
  }
  // すでにソート済み
  void init_sorted(const std::vector<Key> &v){
    if(v.empty()){
      root = nullptr;
      return;
    }
    int n = v.size();
    int size_per_node = limit_size_per_node / 2;
    int m = (n + size_per_node - 1) / size_per_node;
    std::vector<node*> nodes(m);
    for(int i = 0; i < m; i++){
      nodes[i] = new node(v, i * size_per_node, std::min((i + 1) * size_per_node, n));
    }
    root = build(nodes, 0, m);
  }
  int size(){
    return size(root);
  }
  void insert(Key key){
    root = insert_inner(root, key);
  }
  void erase(Key k){
    root = erase_inner(root, k);
  }
  bool find(Key k){
    return find_key(root, k);
  }
  // kの数
  int count(Key k){
    return low_count_eq_inner(root, k) - low_count_inner(root, k);
  }
  Key min(){
    assert(size());
    return leftmost(root)->keys[0];
  }
  Key max(){
    assert(size());
    node *v = rightmost(root);
    return v->keys[v->sz - 1];
  }
  // k未満の値の数
  int low_count(Key k){
    return low_count_inner(root, k);
  }
  // k番目(0-indexed)に小さいキー
  Key kth_smallest(int k){
    if(size() <= k) return inf;
    return find_kth(root, k);
  }
  // k以上の最小要素
  Key lower_bound(Key k){
    return lower_bound_inner(root, k);
  }
  // k以下の最大要素
  Key lower_bound_rev(Key k){
    return lower_bound_rev_inner(root, k);
  }
};

template<typename Key, typename Val>
struct map_avl_cache{
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static Val inf_val;
  static constexpr int limit_size_per_node = 32;
private:
  struct node{
    int h, sz, sz_sum;
    std::array<Key, limit_size_per_node> keys;
    std::array<Val, limit_size_per_node> vals;
    node *l, *r;
    node(): h(1), l(nullptr), r(nullptr){}
    node(Key _key, Val _val): h(1), sz(1), sz_sum(1), l(nullptr), r(nullptr){keys[0] = _key; vals[0] = _val;}
    node(const std::vector<std::pair<Key, Val>> &v, int l, int r): h(1), sz(r - l), sz_sum(sz), l(nullptr), r(nullptr){
      assert(sz < limit_size_per_node);
      for(int i = 0; i < sz; i++){
        keys[i] = v[l + i].first;
        vals[i] = v[l + i].second;
      }
    }
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
    node *split_half(){
      assert(sz == limit_size_per_node);
      node *u = new node();
      sz = limit_size_per_node / 2;
      u->sz_sum = u->sz = limit_size_per_node - sz;
      for(int i = 0; i < u->sz; i++){
        u->keys[i] = keys[sz + i];
        u->vals[i] = vals[sz + i];
      }
      return u;
    }
  };
  node *root, *tmp_node;
  int size(node *v){return v ? v->sz_sum : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz_sum = (v->l ? v->l->sz_sum : 0) + (v->r ? v->r->sz_sum : 0) + v->sz;
  }
  node *rotate_right(node *v){
    node *l = v->l;
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  node *rotate_left(node *v){
    node *r = v->r;
    v->r = r->l;
    r->l = v;
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
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  node *leftmost(node *v){
    while(v->l) v = v->l;
    return v;
  }
  node *rightmost(node *v){
    while(v->r) v = v->r;
    return v;
  }
  std::pair<Key, Val> find_kth(node *v, int k){
    while(true){
      int szl = size(v->l);
      if(szl <= k){
        if(szl + v->sz > k) return std::make_pair(v->keys[k - szl], v->vals[k - szl]);
        k -= szl + v->sz;
        v = v->r;
      }else v = v->l;
    }
  }
  Val find_key(node *v, Key k){
    while(v){
      if(k < v->keys[0]) v = v->l;
      else if(k > v->keys[v->sz - 1]) v = v->r;
      else{
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        if(idx < v->sz && v->keys[idx] == k) return v->vals[idx];
        return inf_val;
      }
    }
    return inf_val;
  }
  node *cut_leftmost(node *v){
    if(v->l){
      v->l = cut_leftmost(v->l);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  node *cut_rightmost(node *v){
    if(v->r){
      v->r = cut_rightmost(v->r);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->l;
  }
  node *insert_leftmost(node *v, node *u){
    if(!v) return u;
    v->l = insert_leftmost(v->l, u);
    update(v);
    return balance(v);
  }
  node *emplace_inner(node *v, Key k, Val val, bool replace = false){
    if(!v) return new node(k, val);
    if(v->l && k < v->keys[0]){
      v->l = emplace_inner(v->l, k, val, replace);
    }else if(v->r && v->keys[v->sz - 1] < k){
      v->r = emplace_inner(v->r, k, val, replace);
    }else{
      int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
      if(idx < v->sz && v->keys[idx] == k){
        if(replace) v->vals[idx] = val;
        return v;
      }
      for(int i = v->sz; i > idx; i--){
        v->keys[i] = v->keys[i - 1];
        v->vals[i] = v->vals[i - 1];
      }
      v->keys[idx] = k;
      v->vals[idx] = val;
      v->sz++;
      if(v->sz == limit_size_per_node){
        v->r = insert_leftmost(v->r, v->split_half());
      }
    }
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, Key k){
    if(!v) return nullptr;
    if(k < v->keys[0]){
      v->l = erase_inner(v->l, k);
    }else if(v->keys[v->sz - 1] < k){
      v->r = erase_inner(v->r, k);
    }else{
      for(int i = 0; i < v->sz; i++){
        if(v->keys[i] == k){
          for(int j = i + 1; j < v->sz; j++){
            v->keys[j - 1] = v->keys[j];
            v->vals[j - 1] = v->vals[i];
          }
          v->sz--;
          break;
        }
      }
      if(!v->sz){
        if(v->r){
          v->r = cut_leftmost(v->r);
          tmp_node->l = v->l;
          tmp_node->r = v->r;
          update(tmp_node);
          return balance(tmp_node);
        }else{
          return v->l;
        }
      }
    }
    update(v);
    return balance(v);
  }
  int low_count_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(k < v->keys[0]) v = v->l;
      else if(k > v->keys[v->sz - 1]) res += szl + v->sz, v = v->r;
      else{
        res += szl;
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        return res + idx;
      }
    }
    return res;
  }
  std::pair<Key, Val> lower_bound_inner(node *v, Key k){
    Key res = inf;
    Val res_val = inf_val;
    while(v){
      if(k < v->keys[0]){
        res = v->keys[0];
        res_val = v->vals[0];
        v = v->l;
      }else if(k > v->keys[v->sz - 1]){
        v = v->r;
      }else{
        int idx = std::lower_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        assert(idx < v->sz);
        return std::make_pair(v->keys[idx], v->vals[idx]);
      }
    }
    return std::make_pair(res, res_val);
  }
  std::pair<Key, Val> lower_bound_rev_inner(node *v, Key k){
    Key res = inf;
    Val res_val = inf_val;
    while(v){
      if(k < v->keys[0]){
        v = v->l;
      }else if(k > v->keys[v->sz - 1]){
        res = v->keys[v->sz - 1];
        res_val = v->vals[v->sz - 1];
        v = v->r;
      }else{
        int idx = std::upper_bound(v->keys.begin(), v->keys.begin() + v->sz, k) - v->keys.begin();
        assert(idx != 0);
        return std::make_pair(v->keys[idx - 1], v->vals[idx - 1]);
      }
    }
    return std::make_pair(res, res_val);
  }
public:
  map_avl_cache(): root(nullptr){}
  map_avl_cache(std::vector<std::pair<Key, Val>> v){
    std::sort(v.begin(), v.end());
    init_sorted(v);
  }
  // すでにソート済み
  void init_sorted(const std::vector<std::pair<Key, Val>> &_v){
    if(_v.empty()){
      root = nullptr;
      return;
    }
    std::vector<std::pair<Key, Val>> v;
    // キーがユニークでない場合前にある要素を優先
    for(int i = 0; i < _v.size(); i++) if(v.empty() || v.back().first != _v[i].first) v.push_back(_v[i]);
    int n = v.size();
    int size_per_node = limit_size_per_node / 2;
    int m = (n + size_per_node - 1) / size_per_node;
    std::vector<node*> nodes(m);
    for(int i = 0; i < m; i++){
      nodes[i] = new node(v, i * size_per_node, std::min((i + 1) * size_per_node, n));
    }
    root = build(nodes, 0, m);
  }
  int size(){
    return size(root);
  }
  void emplace(Key key, Val val){
    root = emplace_inner(root, key, val);
  }
  void emplace_replace(Key key, Val val){
    root = emplace_inner(root, key, val, true);
  }
  void erase(Key k){
    root = erase_inner(root, k);
  }
  bool find(Key k){
    return find_key(root, k) != inf_val;
  }
  Val at(Key k){
    return find_key(root, k);
  }
  std::pair<Key, Val> min(){
    assert(size());
    node *v = leftmost(root);
    return std::make_pair(v->keys[0], v->vals[0]);
  }
  std::pair<Key, Val> max(){
    assert(size());
    node *v = rightmost(root);
    return std::make_pair(v->key[v->sz - 1], v->vals[v->sz - 1]);
  }
  // k未満の値の数
  int low_count(Key k){
    return low_count_inner(root, k);
  }
  // k番目(0-indexed)に小さいキー
  std::pair<Key, Val> kth_smallest(int k){
    if(size() <= k) return std::make_pair(inf, inf_val);
    return find_kth(root, k);
  }
  // k以上の最小要素
  std::pair<Key, Val> lower_bound(Key k){
    return lower_bound_inner(root, k);
  }
  // k以下の最大要素
  std::pair<Key, Val> lower_bound_rev(Key k){
    return lower_bound_rev_inner(root, k);
  }
};
template<typename Key, typename Val>
Val map_avl_cache<Key, Val>::inf_val = std::numeric_limits<Val>::max();
#endif