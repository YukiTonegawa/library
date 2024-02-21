#ifndef _MULTISET_SUM_H
#define _MULTISET_SUM_H

#include <cassert>
#include <limits>
#include <vector>

template<typename Key, typename Val>
struct map_sum_cache{
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static constexpr Val inf_val = std::numeric_limits<Val>::max();
  static constexpr int limit_size_per_node = 16;
private:
  struct node{
    int h, sz, sz_sum;
    std::array<Key, limit_size_per_node> keys;
    std::array<Val, limit_size_per_node> vals;
    Val sum, sumsub;
    node *l, *r;
    node(): h(1), l(nullptr), r(nullptr), sum(0){}
    node(Key _key, Val _val): h(1), sz(1), sz_sum(1), l(nullptr), r(nullptr){keys[0] = _key; vals[0] = sum = sumsub = _val;}
    node(const std::vector<std::pair<Key, Val>> &v, int l, int r): h(1), sz(r - l), sz_sum(sz), sum(0), l(nullptr), r(nullptr){
      assert(sz < limit_size_per_node);
      for(int i = 0; i < sz; i++){
        keys[i] = v[l + i].first;
        vals[i] = v[l + i].second;
        sum += vals[i];
      }
      sumsub = sum;
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
        u->sum += u->vals[i];
      }
      u->sumsub = u->sum;
      sum -= u->sum;
      return u;
    }
  };
  node *root;
  int size(node *v){return v ? v->sz_sum : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz_sum = (v->l ? v->l->sz_sum : 0) + (v->r ? v->r->sz_sum : 0) + v->sz;
    v->sumsub = (v->l ? v->l->sumsub : 0) + (v->r ? v->r->sumsub : 0) + v->sum;
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
        v->vals[idx] += val;
        v->sum += val;
        v->sumsub += val;
        return v;
      }
      for(int i = v->sz; i > idx; i--){
        v->keys[i] = v->keys[i - 1];
        v->vals[i] = v->vals[i - 1];
      }
      v->keys[idx] = k;
      v->vals[idx] = val;
      v->sum += val;
      v->sz++;
      if(v->sz == limit_size_per_node){
        v->r = insert_leftmost(v->r, v->split_half());
      }
    }
    update(v);
    return balance(v);
  }
  Val query_inner(node *v, Key k){
    Val ret = 0;
    while(v){
      if(k < v->keys[0]){
        v = v->l;
      }else if(k > v->keys[v->sz - 1]){
        ret += (v->l ? v->l->sumsub : 0) + v->sum;
        v = v->r;
      }else{
        ret += (v->l ? v->l->sumsub : 0);
        for(int i = 0; i < v->sz; i++){
          if(v->keys[i] >= k) return ret;
          ret += v->vals[i];
        }
      }
    }
    return ret;
  }
public:
  map_sum_cache(): root(nullptr){}
  map_sum_cache(std::vector<std::pair<Key, Val>> v){
    std::sort(v.begin(), v.end());
    init_sorted(v);
  }
  // すでにソート済みかつキーがユニーク
  void init_sorted(const std::vector<std::pair<Key, Val>> &v){
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
  void update(Key key, Val val){
    root = emplace_inner(root, key, val);
  }
  Val query(Key l, Key r){
    return query_inner(root, r) - query_inner(root, l);
  }
};
#endif