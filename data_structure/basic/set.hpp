#ifndef _SET_H_
#define _SET_H_
#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>

template<typename Key>
struct set_avl{
  static constexpr Key inf = std::numeric_limits<Key>::max();
private:
  struct node{
    int h, sz;
    Key key;
    node *l, *r;
    node(Key _key): h(1), sz(1), key(_key), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root, *tmp_node;
  int size(node *v){return v ? v->sz : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz = (v->l ? v->l->sz : 0) + (v->r ? v->r->sz : 0) + 1;
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
  node *find_kth(node *v, int k){
    while(true){
      int szl = size(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  node *find_key(node *v, Key k){
    while(v){
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r;
      else return v;
    }
    return v;
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
  node *insert_inner(node *v, Key k){
    if(!v) return new node(k);
    if(k < v->key){
      v->l = insert_inner(v->l, k);
    }else if(k > v->key){
      v->r = insert_inner(v->r, k);
    }else{
      return v; // set
    }
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, Key k){
    if(!v) return nullptr;
    if(k < v->key){
      v->l = erase_inner(v->l, k);
    }else if(k > v->key){
      v->r = erase_inner(v->r, k);
    }else{
      if(v->r){
        v->r = cut_leftmost(v->r);
        tmp_node->l = v->l;
        tmp_node->r = v->r;
        update(tmp_node);
        return balance(tmp_node);
      }
      return v->l;
    }
    update(v);
    return balance(v);
  }
  int low_count_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r, res += szl + 1;
      else return res + szl;
    }
    return res;
  }
  Key lower_bound_inner(node *v, Key k){
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
  Key lower_bound_rev_inner(node *v, Key k){
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
  void to_list_inner(node *v, std::vector<Key> &res){
    if(v->l) to_list_inner(v->l, res);
    res.push_back(v->key);
    if(v->r) to_list_inner(v->r, res);
  }
public:
  set_avl(): root(nullptr){}
  set_avl(std::vector<Key> v){
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
    std::vector<node*> nodes;
    for(int i = 0; i < n; i++) if(nodes.empty() || nodes.back()->key != v[i]) nodes.push_back(new node(v[i]));
    root = build(nodes, 0, nodes.size());
  }
  int size(){
    return size(root);
  }
  bool empty(){
    return size(root) == 0;
  }
  void insert(Key k){
    root = insert_inner(root, k);
  }
  void erase(Key k){
    root = erase_inner(root, k);
  }
  void clear(){
    root = nullptr;
  }
  bool find(Key k){
    return find_key(root, k);
  }
  Key min(){
    assert(size());
    return leftmost(root)->key;
  }
  Key max(){
    assert(size());
    return rightmost(root)->key;
  }
  // k未満の値の数
  int low_count(Key k){
    return low_count_inner(root, k);
  }
  // k番目(0-indexed)に小さいキー
  Key kth_smallest(int k){
    if(size() <= k) return inf;
    return find_kth(root, k)->key;
  }
  // k以上の最小要素
  Key lower_bound(Key k){
    return lower_bound_inner(root, k);
  }
  // k以下の最大要素
  Key lower_bound_rev(Key k){
    return lower_bound_rev_inner(root, k);
  }
  std::vector<Key> to_list(){
    std::vector<Key> res;
    if(root) to_list_inner(root, res);
    return res;
  }
};

template<typename Key>
struct multiset_avl{
  static constexpr Key inf = std::numeric_limits<Key>::max();
  using Count = int;
private:
  struct node{
    int h, sz_unique;
    Count sz_sum, cnt;
    Key key;
    node *l, *r;
    node(Key _key, Count _cnt): h(1), sz_unique(1), sz_sum(_cnt), cnt(_cnt), key(_key), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root, *tmp_node;
  int size_unique(node *v){return v ? v->sz_unique : 0;}
  Count size_sum(node *v){return v ? v->sz_sum : 0;}

  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz_unique = (v->l ? v->l->sz_unique : 0) + (v->r ? v->r->sz_unique : 0) + 1;
    v->sz_sum = (v->l ? v->l->sz_sum : 0) + (v->r ? v->r->sz_sum : 0) + v->cnt;
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
  node *find_kth_unique(node *v, int k){
    while(true){
      int szl = size_unique(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  node *find_kth_sum(node *v, Count k){
    while(true){
      Count szl = size_sum(v->l);
      if(szl <= k){
        if(k < szl + v->cnt) return v;
        k -= szl + v->cnt;
        v = v->r;
      }else v = v->l;
    }
  }
  node *find_key(node *v, Key k){
    while(v){
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r;
      else return v;
    }
    return v;
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
  node *insert_inner(node *v, Key k, Count cnt){
    if(!v) return new node(k, cnt);
    if(k < v->key){
      v->l = insert_inner(v->l, k, cnt);
    }else if(k > v->key){
      v->r = insert_inner(v->r, k, cnt);
    }else{
      v->cnt += cnt; // multiset
      update(v);
      return v;
    }
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, Key k, Count cnt){
    if(!v) return nullptr;
    if(k < v->key){
      v->l = erase_inner(v->l, k, cnt);
    }else if(k > v->key){
      v->r = erase_inner(v->r, k, cnt);
    }else{
      v->cnt -= cnt;
      if(v->cnt <= 0){
        if(v->r){
          v->r = cut_leftmost(v->r);
          tmp_node->l = v->l;
          tmp_node->r = v->r;
          update(tmp_node);
          return balance(tmp_node);
        }
        return v->l;
      }
    }
    update(v);
    return balance(v);
  }
  int low_count_unique_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size_unique(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) res += szl + 1, v = v->r;
      else return res + szl;
    }
    return res;
  }
  Count low_count_sum_inner(node *v, Key k){
    Count res = 0;
    while(v){
      int szl = size_sum(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) res += szl + v->cnt, v = v->r;
      else return res + szl;
    }
    return res;
  }
  Key lower_bound_inner(node *v, Key k){
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
  Key lower_bound_rev_inner(node *v, Key k){
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
  void to_list_inner(node *v, std::vector<std::pair<Key, Count>> &res){
    if(v->l) to_list_inner(v->l, res);
    res.push_back({v->key, v->cnt});
    if(v->r) to_list_inner(v->r, res);
  }
public:
  multiset_avl(): root(nullptr){}
  multiset_avl(std::vector<Key> v){
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
    std::vector<node*> nodes;
    for(int i = 0; i < n; i++){
      if(nodes.empty() || nodes.back()->key != v[i]){
        nodes.push_back(new node(v[i], 1));
      }else{
        nodes.back()->cnt++;
        nodes.back()->sz_sum++;
      }
    }
    root = build(nodes, 0, nodes.size());
  }
  // 要素数
  Count size(){
    return size_sum(root);
  }
  bool empty(){
    return size_sum(root) == 0;
  }
  // 種類数
  int size_unique(){
    return size_unique(root);
  }
  // kをcnt個追加
  void insert(Key k, Count cnt = 1){
    root = insert_inner(root, k, cnt);
  }
  // kをcnt個削除(足りない分は無視)
  void erase(Key k, Count cnt = 1){
    root = erase_inner(root, k, cnt);
  }
  void clear(){
    root = nullptr;
  }
  bool find(Key k){
    return find_key(root, k);
  }
  // kの個数
  Count count(Key k){
    node *v = find_key(root, k);
    return v ? v->cnt : 0;
  }
  Key min(){
    assert(size());
    return leftmost(root)->key;
  }
  Key max(){
    assert(size());
    return rightmost(root)->key;
  }

  std::pair<Key, Count> min2(){
    assert(size());
    auto v = leftmost(root);
    return {v->key, v->cnt};
  }
  std::pair<Key, Count> max2(){
    assert(size());
    auto v = rightmost(root);
    return {v->key, v->cnt};
  }
  // k未満の値の種類
  int low_count_unique(Key k){
    return low_count_unique_inner(root, k);
  }
  // k未満の値の数(同じ要素を重複して数える)
  Count low_count_sum(Key k){
    return low_count_sum_inner(root, k);
  }
  // k番目(0-indexed)に小さいキーの種類, ない場合はinf
  Key kth_smallest_unique(int k){
    if(size_unique() <= k) return inf;
    return find_kth_unique(root, k)->key;
  }
  // k番目(0-indexed)に小さいキー(同じ要素を重複して数える), ない場合はinf
  Key kth_smallest_sum(Count k){
    if(size() <= k) return inf;
    return find_kth_sum(root, k)->key;
  }
  // k以上の最小要素, ない場合はinf
  Key lower_bound(Key k){
    return lower_bound_inner(root, k);
  }
  // k以下の最大要素, ない場合はinf
  Key lower_bound_rev(Key k){
    return lower_bound_rev_inner(root, k);
  }
  std::vector<std::pair<Key, Count>> to_list(){
    std::vector<std::pair<Key, Count>> res;
    if(root) to_list_inner(root, res);
    return res;
  }
};

template<typename Key, typename Val>
struct map_avl{
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static Val inf_val;
private:
  struct node{
    int h, sz;
    Key key;
    Val val;
    node *l, *r;
    node(Key _key, Val _val): h(1), sz(1), key(_key), val(_val), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root, *tmp_node;
  int size(node *v){return v ? v->sz : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz = (v->l ? v->l->sz : 0) + (v->r ? v->r->sz : 0) + 1;
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
  node *find_kth(node *v, int k){
    while(true){
      int szl = size(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  node *find_key(node *v, Key k){
    while(v){
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r;
      else return v;
    }
    return v;
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
  node *emplace_inner(node *v, Key k, Val val, bool replace = false){
    if(!v) return new node(k, val);
    if(k < v->key){
      v->l = emplace_inner(v->l, k, val, replace);
    }else if(k > v->key){
      v->r = emplace_inner(v->r, k, val, replace);
    }else{
      if(replace) v->val = val;
      return v;
    }
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, Key k){
    if(!v) return nullptr;
    if(k < v->key){
      v->l = erase_inner(v->l, k);
    }else if(k > v->key){
      v->r = erase_inner(v->r, k);
    }else{
      if(v->r){
        v->r = cut_leftmost(v->r);
        tmp_node->l = v->l;
        tmp_node->r = v->r;
        update(tmp_node);
        return balance(tmp_node);
      }
      return v->l;
    }
    update(v);
    return balance(v);
  }
  int low_count_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) res += szl + 1, v = v->r;
      else return res + szl;
    }
    return res;
  }
  node *lower_bound_inner(node *v, Key k){
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
  node *lower_bound_rev_inner(node *v, Key k){
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
  void to_list_inner(node *v, std::vector<std::pair<Key, Val>> &res){
    if(v->l) to_list_inner(v->l, res);
    res.push_back({v->key, v->val});
    if(v->r) to_list_inner(v->r, res);
  }
public:
  map_avl(): root(nullptr){}
  map_avl(std::vector<std::pair<Key, Val>> v){
    std::sort(v.begin(), v.end());
    init_sorted(v);
  }
  // すでにソート済み
  // キーがユニークでない場合前にあるものを優先
  void init_unsafe(const std::vector<std::pair<Key, Val>> &v){
    if(v.empty()){
      root = nullptr;
      return;
    }
    int n = v.size();
    std::vector<node*> nodes;
    for(int i = 0; i < n; i++){
      if(nodes.empty() || nodes.back()->key != v[i].first){
        nodes.push_back(new node(v[i].first, v[i].second));
      }
    }
    root = build(nodes, 0, nodes.size());
  }
  int size(){
    return size(root);
  }
  bool empty(){
    return size(root) == 0;
  }
  void emplace(Key k, Val val){
    root = emplace_inner(root, k, val);
  }
  void emplace_replace(Key k, Val val){
    root = emplace_inner(root, k, val, true);
  }
  void erase(Key k){
    root = erase_inner(root, k);
  }
  void clear(){
    root = nullptr;
  }
  bool find(Key k){
    return find_key(root, k);
  }
  Val at(Key k){
    node *v = find_key(root, k);
    return v ? v->val : inf_val;
  }
  std::pair<Key, Val> min(){
    assert(size());
    node *v = leftmost(root);
    return std::make_pair(v->key, v->val);
  }
  std::pair<Key, Val> max(){
    assert(size());
    node *v = rightmost(root);
    return std::make_pair(v->key, v->val);
  }
  // k未満の値の数
  int low_count(Key k){
    return low_count_inner(root, k);
  }
  // k番目(0-indexed)に小さいキー
  std::pair<Key, Val> kth_smallest(int k){
    if(size() <= k) return std::make_pair(inf, inf_val);
    node *v = find_kth(root, k);
    return std::make_pair(v->key, v->val);
  }
  // k以上の最小要素
  std::pair<Key, Val> lower_bound(Key k){
    node *v = lower_bound_inner(root, k);
    return v ? std::make_pair(v->key, v->val) : std::make_pair(inf, inf_val);
  }
  // k以下の最大要素
  std::pair<Key, Val> lower_bound_rev(Key k){
    node *v = lower_bound_rev_inner(root, k);
    return v ? std::make_pair(v->key, v->val) : std::make_pair(inf, inf_val);
  }
  std::vector<std::pair<Key, Val>> to_list(){
    std::vector<std::pair<Key, Val>> res;
    if(root) to_list_inner(root, res);
    return res;
  }
};
template<typename Key, typename Val>
Val map_avl<Key, Val>::inf_val = std::numeric_limits<Val>::max();

// キーの昇順, キーが同じものは追加した順番の昇順
template<typename Key, typename Val>
struct multimap_avl{
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static Val inf_val;
private:
  struct node{
    int h, sz_unique, sz_sum;
    Key key;
    std::vector<Val> val;
    node *l, *r;
    node(Key _key, Val _val): h(1), sz_unique(1), sz_sum(1), key(_key), val{_val}, l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root, *tmp_node;
  int size_unique(node *v){return v ? v->sz_unique : 0;}
  int size_sum(node *v){return v ? v->sz_sum : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz_unique = (v->l ? v->l->sz_unique : 0) + (v->r ? v->r->sz_unique : 0) + 1;
    v->sz_sum = (v->l ? v->l->sz_sum : 0) + (v->r ? v->r->sz_sum : 0) + v->val.size();
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
  node *find_kth_unique(node *v, int k){
    while(true){
      int szl = size_unique(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  node *find_kth_sum(node *v, int &k){
    while(true){
      int szl = size_sum(v->l);
      if(szl <= k){
        if(k < szl + v->cnt){
          k -= szl;
          return v;
        }
        k -= szl + v->cnt;
        v = v->r;
      }else v = v->l;
    }
  }
  node *find_key(node *v, Key k){
    while(v){
      if(k < v->key) v = v->l;
      else if(k > v->key) v = v->r;
      else return v;
    }
    return v;
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
  node *emplace_inner(node *v, Key k, Val val){
    if(!v) return new node(k, val);
    if(k < v->key){
      v->l = emplace_inner(v->l, k, val);
    }else if(k > v->key){
      v->r = emplace_inner(v->r, k, val);
    }else{
      v->val.push_back(val);
      update(v);
      return v;
    }
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, Key k, bool erase_all = false){
    if(!v) return nullptr;
    if(k < v->key){
      v->l = erase_inner(v->l, k, erase_all);
    }else if(k > v->key){
      v->r = erase_inner(v->r, k, erase_all);
    }else{
      if(erase_all) v->val.clear();
      else v->val.pop_back();
      if(v->val.empty()){
        if(v->r){
          v->r = cut_leftmost(v->r);
          tmp_node->l = v->l;
          tmp_node->r = v->r;
          update(tmp_node);
          return balance(tmp_node);
        }
        return v->l;
      }
    }
    update(v);
    return balance(v);
  }
  int low_count_unique_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size_unique(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) res += szl + 1, v = v->r;
      else return res + szl;
    }
    return res;
  }
  int low_count_sum_inner(node *v, Key k){
    int res = 0;
    while(v){
      int szl = size_sum(v->l);
      if(k < v->key) v = v->l;
      else if(k > v->key) res += szl + v->val.size(), v = v->r;
      else return res + szl;
    }
    return res;
  }
  node *lower_bound_inner(node *v, Key k){
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
  node *lower_bound_rev_inner(node *v, Key k){
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
  multimap_avl(): root(nullptr){}
  multimap_avl(std::vector<std::pair<Key, Val>> v){
    std::sort(v.begin(), v.end());
    init_sorted(v);
  }
  // すでにソート済み
  void init_sorted(const std::vector<std::pair<Key, Val>> &v){
    if(v.empty()){
      root = nullptr;
      return;
    }
    std::vector<node*> nodes;
    for(int i = 0; i < v.size(); i++){
      if(nodes.empty() || nodes.back().key != v[i].first){
        nodes.push_back(new node(v[i].first, v[i].second));
      }else{
        nodes.back()->val.push_back(v[i].second);
      }
    }
    root = build(nodes, 0, nodes.size());
  }
  // 要素数
  int size(){
    return size_sum(root);
  }
  bool empty(){
    return size_sum(root) == 0;
  }
  // 種類数
  int size_unique(){
    return size_unique(root);
  }
  void emplace(Key k, Val val){
    root = emplace_inner(root, k, val);
  }
  // 同じキーを持つ要素が複数ある場合, 最後に追加された要素を消す
  void erase(Key k){
    root = erase_inner(root, k);
  }
  // kをキーとする要素を全て消す
  void erase_all(Key k){
    root = erase_inner(root, k, true);
  }
  void clear(){
    root = nullptr;
  }
  bool find(Key k){
    return find_key(root, k);
  }
  // 複数ある場合は最後に追加した要素
  Val at(Key k){
    node *v = find_key(root, k);
    return v ? v->val.back() : inf_val;
  }
  std::pair<Key, Val> min(){
    assert(size());
    node *v = leftmost(root);
    return std::make_pair(v->key, v->val[0]);
  }
  // 複数ある場合は最後に追加した要素
  std::pair<Key, Val> max(){
    assert(size());
    node *v = rightmost(root);
    return std::make_pair(v->key, v->val.back());
  }
  // k未満の値の種類
  int low_count_unique(Key k){
    return low_count_unique_inner(root, k);
  }
  // k未満の値の数(同じ要素を重複して数える)
  int low_count_sum(Key k){
    return low_count_sum_inner(root, k);
  }
  // k種類目(0-indexed)に小さいキーの最初に追加された要素
  std::pair<Key, Val> kth_smallest_unique(int k){
    if(size_unique() <= k) return std::make_pair(inf, inf_val);
    node *v = find_kth_unique(root, k);
    return std::make_pair(v->key, v->val[0]);
  }
  // k番目(0-indexed)に小さいキー(同じ要素を重複して数える)
  std::pair<Key, Val> kth_smallest_sum(int k){
    if(size() <= k) return std::make_pair(inf, inf_val);
    node *v = find_kth_sum(root, k);
    return std::make_pair(v->key, v->val[k]);
  }
  // k以上の最小要素
  std::pair<Key, Val> lower_bound(Key k){
    node *v = lower_bound_inner(root, k);
    return v ? std::make_pair(v->key, v->val[0]) : std::make_pair(inf, inf_val);
  }
  // k以下の最大要素
  std::pair<Key, Val> lower_bound_rev(Key k){
    node *v = lower_bound_rev_inner(root, k);
    return v ? std::make_pair(v->key, v->val.back()) : std::make_pair(inf, inf_val);
  }
};
template<typename Key, typename Val>
Val multimap_avl<Key, Val>::inf_val = std::numeric_limits<Val>::max();
#endif
