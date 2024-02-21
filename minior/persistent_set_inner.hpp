#ifndef _PERSISTENT_SET_INNER_H_
#define _PERSISTENT_SET_INNER_H_
#include <vector>
#include <cassert>
#include <numeric>

template<typename Key>
struct persistent_set{
public:
  static constexpr Key inf = std::numeric_limits<Key>::max();
  struct node{
    char h;
    int sz;
    Key key;
    node *l, *r;
    node(Key _key): h(1), sz(1), key(_key), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  static int size(node *v){return v ? v->sz : 0;}
private:
  static node *tmp_node;
  static node *copy_node(node *v){
    if(!v) return nullptr;
    return new node(*v);
  }
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz = (v->l ? v->l->sz : 0) + (v->r ? v->r->sz : 0) + 1;
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
  static node *balance_insert(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    if(bf == 2){
      if(v->l->balanace_factor() == -1) v->l = rotate_left(v->l);
      return rotate_right(v);
    }else if(bf == -2){
      if(v->r->balanace_factor() == 1) v->r = rotate_right(v->r);
      return rotate_left(v);
    }
    return v;
  }
  static node *balance_erase(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    if(bf == 2){
      v->l = copy_node(v->l);
      if(v->l->balanace_factor() == -1){
        v->l->r = copy_node(v->l->r);
        v->l = rotate_left(v->l);
      }
      return rotate_right(v);
    }else if(bf == -2){
      v->r = copy_node(v->r);
      if(v->r->balanace_factor() == 1){
        v->r->l = copy_node(v->r->l);
        v->r = rotate_right(v->r);
      }
      return rotate_left(v);
    }
    return v;
  }
  static node *leftmost(node *v){
    while(v->l) v = v->l;
    return v;
  }
  static node *rightmost(node *v){
    while(v->r) v = v->r;
    return v;
  }
  static node *cut_leftmost(node *v){
    v = copy_node(v);
    if(v->l){
      v->l = cut_leftmost(v->l);
      update(v);
      return balance_erase(v);
    }
    tmp_node = v;
    return v->r;
  }
  static node *cut_rightmost(node *v){
    v = copy_node(v);
    if(v->r){
      v->r = cut_rightmost(v->r);
      update(v);
      return balance_erase(v);
    }
    tmp_node = v;
    return v->l;
  }
  static node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  static node *find_key(node *v, Key k){
    while(v){
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r;
      else return v;
    }
    return v;
  }
  static node *find_kth(node *v, int k){
    while(true){
      int szl = size(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
public:
  persistent_set(){}
  // すでにソート済み
  static node *build_sorted(const std::vector<Key> &v){
    if(v.empty()) return nullptr;
    int n = v.size();
    std::vector<node*> nodes;
    for(int i = 0; i < n; i++) if(nodes.empty() || nodes.back()->key != v[i]) nodes.push_back(new node(v[i]));
    return build(nodes, 0, nodes.size());
  }
  static node *insert(node *v, Key x){
    if(!v) return new node(x);
    v = copy_node(v);
    if(x < v->key) v->l = insert(v->l, x);
    else if(x > v->key) v->r = insert(v->r, x);
    else return v;
    update(v);
    return balance_insert(v);
  }
  static node *erase(node *v, Key x){
    if(!v) return v;
    if(x < v->key){
      v = copy_node(v);
      v->l = erase(v->l, x);
    }
    else if(x > v->key){
      v = copy_node(v);
      v->r = erase(v->r, x);
    }else{
      if(!v->r) return copy_node(v->l);
      node *u = cut_leftmost(v->r);
      tmp_node->l = v->l;
      tmp_node->r = u;
      v = tmp_node;
    }
    update(v);
    return balance_erase(v);
  }
  static bool empty(node *v){
    return size(v) == 0;
  }
  static bool find(node *v, Key k){
    return find_key(v, k);
  }
  static Key min(node *v){
    assert(size(v));
    return leftmost(v)->key;
  }
  static Key max(node *v){
    assert(size(v));
    return rightmost(v)->key;
  }
  static int low_count(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r, res += szl + 1;
      else return res + szl;
    }
    return res;
  }
  static Key lower_bound(node *v, Key k){
    Key res = inf;
    while(v){
      if(k < v->key){
        res = v->key;
        v = v->l;
      }else if(k > v->key){
        v = v->r;
      }else return k;
    }
    return res;
  }
  static Key lower_bound_rev(node *v, Key k){
    Key res = inf;
    while(v){
      if(k < v->key){
        v = v->l;
      }else if(k > v->key){
        res = v->key;
        v = v->r;
      }else return k;
    }
    return res;
  }
  static Key kth_smallest(node *v, int k){
    if(size(v) <= k) return inf;
    return find_kth(v, k)->key;
  }
};
template<typename T>
typename persistent_set<T>::node *persistent_set<T>::tmp_node = nullptr;


template<typename Key>
struct persistent_multiset{
public:
  static constexpr Key inf = std::numeric_limits<Key>::max();
  using Count = int;
  struct node{
    char h;
    int sz_unique;
    Count sz_sum, cnt;
    Key key;
    node *l, *r;
    node(Key _key, Count _cnt): h(1), sz_unique(1), sz_sum(_cnt), cnt(_cnt), key(_key), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  static int size_unique(node *v){return v ? v->sz_unique : 0;}
  static Count size_sum(node *v){return v ? v->sz_sum : 0;}
private:
  static node *tmp_node;
  static node *copy_node(node *v){
    if(!v) return nullptr;
    return new node(*v);
  }
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz_unique = (v->l ? v->l->sz_unique : 0) + (v->r ? v->r->sz_unique : 0) + 1;
    v->sz_sum = (v->l ? v->l->sz_sum : 0) + (v->r ? v->r->sz_sum : 0) + v->cnt;
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
  static node *balance_insert(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    if(bf == 2){
      if(v->l->balanace_factor() == -1) v->l = rotate_left(v->l);
      return rotate_right(v);
    }else if(bf == -2){
      if(v->r->balanace_factor() == 1) v->r = rotate_right(v->r);
      return rotate_left(v);
    }
    return v;
  }
  static node *balance_erase(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    if(bf == 2){
      v->l = copy_node(v->l);
      if(v->l->balanace_factor() == -1){
        v->l->r = copy_node(v->l->r);
        v->l = rotate_left(v->l);
      }
      return rotate_right(v);
    }else if(bf == -2){
      v->r = copy_node(v->r);
      if(v->r->balanace_factor() == 1){
        v->r->l = copy_node(v->r->l);
        v->r = rotate_right(v->r);
      }
      return rotate_left(v);
    }
    return v;
  }
  static node *leftmost(node *v){
    while(v->l) v = v->l;
    return v;
  }
  static node *rightmost(node *v){
    while(v->r) v = v->r;
    return v;
  }
  static node *cut_leftmost(node *v){
    v = copy_node(v);
    if(v->l){
      v->l = cut_leftmost(v->l);
      update(v);
      return balance_erase(v);
    }
    tmp_node = v;
    return v->r;
  }
  static node *cut_rightmost(node *v){
    v = copy_node(v);
    if(v->r){
      v->r = cut_rightmost(v->r);
      update(v);
      return balance_erase(v);
    }
    tmp_node = v;
    return v->l;
  }
  static node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  static node *find_key(node *v, Key k){
    while(v){
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r;
      else return v;
    }
    return v;
  }
  static node *find_kth_unique(node *v, int k){
    while(true){
      int szl = size_unique(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  static node *find_kth_sum(node *v, Count k){
    while(true){
      Count szl = size_sum(v->l);
      if(szl <= k){
        if(k < szl + v->cnt) return v;
        k -= szl + v->cnt;
        v = v->r;
      }else v = v->l;
    }
  }
public:
  persistent_multiset(){}
  // すでにソート済み
  static node *build_sorted(const std::vector<Key> &v){
    if(v.empty()) return nullptr;
    int n = v.size();
    std::vector<node*> nodes;
    for(int i = 0; i < n; i++){
      if(nodes.empty() || nodes.back()->key != v[i]){
        nodes.push_back(new node(v[i], 1));
      }else{
        nodes.back()->cnt++;
        nodes.back()->sz_sum++;
      }
    }
    return build(nodes, 0, nodes.size());
  }
  static node *insert(node *v, Key x, int cnt = 1){
    if(!v) return new node(x, cnt);
    v = copy_node(v);
    if(x < v->key) v->l = insert(v->l, x, cnt);
    else if(x > v->key) v->r = insert(v->r, x, cnt);
    else{
      v->cnt += cnt;
    }
    update(v);
    return balance_insert(v);
  }
  static node *erase(node *v, Key x, int cnt = 1){
    if(!v) return v;
    if(x < v->key){
      v = copy_node(v);
      v->l = erase(v->l, x, cnt);
    }
    else if(x > v->key){
      v = copy_node(v);
      v->r = erase(v->r, x, cnt);
    }else{
      if(v->cnt <= cnt){
        if(!v->r) return copy_node(v->l);
        node *u = cut_leftmost(v->r);
        tmp_node->l = v->l;
        tmp_node->r = u;
        v = tmp_node;
      }else{
        v = copy_node(v);
        v->cnt -= cnt;
      }
    }
    update(v);
    return balance_erase(v);
  }
  static bool empty(node *v){
    return size(v) == 0;
  }
  static Count find(node *v, Key k){
    auto u = find_key(v, k);
    return u ? u->cnt : 0;
  }
  static Key min(node *v){
    assert(size(v));
    return leftmost(v)->key;
  }
  static Key max(node *v){
    assert(size(v));
    return rightmost(v)->key;
  }
  static int low_count_unique(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size_unique(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) res += szl + 1, v = v->r;
      else return res + szl;
    }
    return res;
  }
  static Count low_count_sum(node *v, Key k){
    Count res = 0;
    while(v){
      int szl = size_sum(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) res += szl + v->cnt, v = v->r;
      else return res + szl;
    }
    return res;
  }
  static Key lower_bound(node *v, Key k){
    Key res = inf;
    while(v){
      if(k < v->key){
        res = v->key;
        v = v->l;
      }else if(k > v->key){
        v = v->r;
      }else return k;
    }
    return res;
  }
  static Key lower_bound_rev(node *v, Key k){
    Key res = inf;
    while(v){
      if(k < v->key){
        v = v->l;
      }else if(k > v->key){
        res = v->key;
        v = v->r;
      }else return k;
    }
    return res;
  }
  // k番目(0-indexed)に小さいキーの種類, ない場合はinf
  static Key kth_smallest_unique(node *v, int k){
    if(size_unique(v) <= k) return inf;
    return find_kth_unique(v, k)->key;
  }
  // k番目(0-indexed)に小さいキー(同じ要素を重複して数える), ない場合はinf
  static Key kth_smallest_sum(node *v, Count k){
    if(size(v) <= k) return inf;
    return find_kth_sum(v, k)->key;
  }
};
template<typename T>
typename persistent_multiset<T>::node *persistent_multiset<T>::tmp_node = nullptr;

template<typename Key, typename Val>
struct persistent_map{
public:
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static constexpr Val inf_val = std::numeric_limits<Val>::max();
  struct node{
    char h;
    int sz;
    Key key;
    Val val;
    node *l, *r;
    node(Key _key, Val _val): h(1), sz(1), key(_key), val(_val), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  static int size(node *v){return v ? v->sz : 0;}
private:
  static node *tmp_node;
  static node *copy_node(node *v){
    if(!v) return nullptr;
    return new node(*v);
  }
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz = (v->l ? v->l->sz : 0) + (v->r ? v->r->sz : 0) + 1;
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
  static node *balance_insert(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    if(bf == 2){
      if(v->l->balanace_factor() == -1) v->l = rotate_left(v->l);
      return rotate_right(v);
    }else if(bf == -2){
      if(v->r->balanace_factor() == 1) v->r = rotate_right(v->r);
      return rotate_left(v);
    }
    return v;
  }
  static node *balance_erase(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    if(bf == 2){
      v->l = copy_node(v->l);
      if(v->l->balanace_factor() == -1){
        v->l->r = copy_node(v->l->r);
        v->l = rotate_left(v->l);
      }
      return rotate_right(v);
    }else if(bf == -2){
      v->r = copy_node(v->r);
      if(v->r->balanace_factor() == 1){
        v->r->l = copy_node(v->r->l);
        v->r = rotate_right(v->r);
      }
      return rotate_left(v);
    }
    return v;
  }
  static node *leftmost(node *v){
    while(v->l) v = v->l;
    return v;
  }
  static node *rightmost(node *v){
    while(v->r) v = v->r;
    return v;
  }
  static node *cut_leftmost(node *v){
    v = copy_node(v);
    if(v->l){
      v->l = cut_leftmost(v->l);
      update(v);
      return balance_erase(v);
    }
    tmp_node = v;
    return v->r;
  }
  static node *cut_rightmost(node *v){
    v = copy_node(v);
    if(v->r){
      v->r = cut_rightmost(v->r);
      update(v);
      return balance_erase(v);
    }
    tmp_node = v;
    return v->l;
  }
  static node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  static node *find_key(node *v, Key k){
    while(v){
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r;
      else return v;
    }
    return v;
  }
  static node *find_kth(node *v, int k){
    while(true){
      int szl = size(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  static node *lower_bound_inner(node *v, Key k){
    node *res = nullptr;
    while(v){
      if(k < v->key){
        res = v;
        v = v->l;
      }else if(k > v->key){
        v = v->r;
      }else return v;
    }
    return res;
  }
  static node *lower_bound_rev_inner(node *v, Key k){
    node *res = nullptr;
    while(v){
      if(k < v->key){
        v = v->l;
      }else if(k > v->key){
        res = v;
        v = v->r;
      }else return v;
    }
    return res;
  }
public:
  persistent_map(){}
  // すでにソート済み
  // キーがユニークでない場合前にあるものを優先
  static node *build_sorted_unsafe(const std::vector<std::pair<Key, Val>> &v){
    if(v.empty()) return nullptr;
    int n = v.size();
    std::vector<node*> nodes;
    for(int i = 0; i < n; i++){
      if(nodes.empty() || nodes.back()->key != v[i].first){
        nodes.push_back(new node(v[i].first, v[i].second));
      }
    }
    return build(nodes, 0, nodes.size());
  }
  static node *emplace(node *v, Key x, Val val){
    if(!v) return new node(x, val);
    v = copy_node(v);
    if(x < v->key) v->l = emplace(v->l, x, val);
    else if(x > v->key) v->r = emplace(v->r, x, val);
    else return v;
    update(v);
    return balance_insert(v);
  }
  static node *emplace_replace(node *v, Key x, Val val){
    if(!v) return new node(x, val);
    v = copy_node(v);
    if(x < v->key) v->l = emplace_replace(v->l, x, val);
    else if(x > v->key) v->r = emplace_replace(v->r, x, val);
    else v->val = val;
    update(v);
    return balance_insert(v);
  }
  static node *erase(node *v, Key x){
    if(!v) return v;
    if(x < v->key){
      v = copy_node(v);
      v->l = erase(v->l, x);
    }
    else if(x > v->key){
      v = copy_node(v);
      v->r = erase(v->r, x);
    }else{
      if(!v->r) return copy_node(v->l);
      node *u = cut_leftmost(v->r);
      tmp_node->l = v->l;
      tmp_node->r = u;
      v = tmp_node;
    }
    update(v);
    return balance_erase(v);
  }
  static bool empty(node *v){
    return size(v) == 0;
  }
  static bool find(node *v, Key k){
    return find_key(v, k);
  }
  static std::pair<Key, Val> min(node *v){
    assert(size(v));
    auto u = leftmost(v);
    return {u->key, u->val};
  }
  static std::pair<Key, Val> max(node *v){
    assert(size(v));
    auto u = rightmost(v);
    return {u->key, u->val};
  }
  static int low_count(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r, res += szl + 1;
      else return res + szl;
    }
    return res;
  }
  static Val at(node *v, Key k){
    node *u = find_key(v, k);
    return u ? u->val : inf_val;
  }
  static std::pair<Key, Val> lower_bound(node *v, Key k){
    auto u = lower_bound_inner(v, k);
    if(!u) return {inf, inf_val};
    return {u->key, u->val};
  }
  static std::pair<Key, Val> lower_bound_rev(node *v, Key k){
    auto u = lower_bound_rev_inner(v, k);
    if(!u) return {inf, inf_val};
    return {u->key, u->val};
  }
  static std::pair<Key, Val> kth_smallest(node *v, int k){
    if(size(v) <= k) return {inf, inf_val};
    auto u = find_kth(v, k);
    return {u->key, u->val};
  }
};
template<typename T, typename U>
typename persistent_map<T, U>::node *persistent_map<T, U>::tmp_node = nullptr;
#endif