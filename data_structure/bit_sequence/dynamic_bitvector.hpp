#ifndef _DYNAMIC_BITVECTOR_H_
#define _DYNAMIC_BITVECTOR_H_
#include <vector>
#include <iostream>
#include <bitset>
#include <cassert>
#include <algorithm>
#include "bit_operation.hpp"

struct dynamic_bit128{
  static constexpr __uint128_t mask_0_64 = ((__uint128_t)1 << 64) - 1;
  __uint128_t x;
  dynamic_bit128(__uint128_t _x): x(_x){}
  std::pair<__uint128_t, __uint128_t> split(){
    return {x >> 64, x & mask_0_64};
  }
  // 先頭k-bitと残りの間にfを挿入
  void insert(int k, bool f){
    __uint128_t y = ((__uint128_t)1 << k) - 1;
    __uint128_t left = x & ~y, right = x & y;
    x = (left << 1) | ((__uint128_t)f << k) | right;
  }
  // k番目を削除
  bool erase(int k){
    bool ret = (x >> k) & 1;
    if(k == 127){
      x ^= (__uint128_t)ret << k;
    }else{
      __uint128_t y = ((__uint128_t)1 << (k + 1)) - 1;
      __uint128_t left = (x & ~y) >> 1, right = (x & (y >> 1));
      x = left | right;
    }
    return ret;
  }
  // k番目にfをセット
  void set(int k, bool f){
    if(f) x |= ((__uint128_t)1 << k);
    else x &= ~((__uint128_t)1 << k);
  }
  // k番目の値を取得
  bool get(int k){
    return (x >> k) & 1;
  }
  int pop_count(){
    return __builtin_popcountll(x >> 64) + __builtin_popcountll(x & mask_0_64);
  }
  bool access(int k){
    return (x >> k) & 1;
  }
  // rank [0, k) 下位k-bitのpop
  int rank(int k){
    if(k == 128) return pop_count();
    if(k < 64) return __builtin_popcountll((x & mask_0_64) & ((1ULL << k) - 1));
    k -= 64;
    return __builtin_popcountll(x & mask_0_64) + __builtin_popcountll((x >> 64)  & ((1ULL << k) - 1));
  }
  // k番目の1, 無い場合壊れる
  int select1(int k){
    int left_pop = __builtin_popcountll(x & mask_0_64);
    if(left_pop > k) return select_64bit(x & mask_0_64, k);
    return 64 + select_64bit(x >> 64, k - left_pop);
  }
  // k番目の0, 無い場合壊れる
  int select0(int k){
    __uint128_t y = ~x;
    int left_unpop = __builtin_popcountll(y & mask_0_64);
    if(left_unpop > k) return select_64bit(y & mask_0_64, k);
    return 64 + select_64bit(y >> 64, k - left_unpop);
  }
};

struct dynamic_bitvector{
  struct node{
    int h;
    uint8_t sz, pop;
    uint32_t szsum, popsum;
    dynamic_bit128 bit;
    node *l, *r, *p;
    node(__uint128_t x, int sz): h(1), sz(sz), szsum(sz), bit(x), l(nullptr), r(nullptr), p(nullptr){
      pop = popsum = __builtin_popcountll(x);
    }
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  static int size(node *v){return v ? v->szsum : 0;}
  static int pop(node *v){return v ? v->popsum : 0;}
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0, v->r ? v->r->h : 0) + 1;
    v->szsum = v->sz;
    v->popsum = v->pop;
    if(v->l) v->szsum += v->l->szsum, v->popsum += v->l->popsum;
    if(v->r) v->szsum += v->r->szsum, v->popsum += v->r->popsum;
  }
  static void update_rec(node *v){
    while(v){
      update(v);
      v = v->p;
    }
  }
  static node *rotate_right(node *v){
    node *p = v->p, *pp = p->p;
    if((p->l = v->r)) v->r->p = p;
    v->r = p, p->p = v;
    update(p), update(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      update(pp);
    }
    return v;
  }
  static node *rotate_left(node *v){
    node *p = v->p, *pp = p->p;
    if((p->r = v->l)) v->l->p = p;
    v->l = p, p->p = v;
    update(p), update(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      update(pp);
    }
    return v;
  }
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l){
      v->l = build(nodes, l, m);
      v->l->p = v;
    }
    if(r > m + 1){
      v->r = build(nodes, m + 1, r);
      v->r->p = v;
    }
    update(v);
    return v;
  }
  static node *left_most(node *v){
    while(v->l) v = v->l;
    return v;
  }
  static node *right_most(node *v){
    while(v->r) v = v->r;
    return v;
  }
  // 左側のサイズの和が始めてkを超えるようなノード, 余ったkを返す
  static node *find_kth(node *v, int &k){
    while(true){
      int lsz = size(v->l);
      if(lsz <= k){
        int msz = lsz + v->sz;
        if(k < msz){
          k -= lsz;
          return v;
        }
        k -= msz;
        v = v->r;
      }else v = v->l;
    }
  }
  // vから親まで辿って更新しながら平衡化
  static node *balancing(node *v){
    while(true){
      update(v);
      int bf = v->balanace_factor();
      assert(-2 <= bf && bf <= 2);
      if(bf == 2){
        if(v->l->balanace_factor() == -1){
          rotate_left(v->l->r);
          v = rotate_right(v->l);
        }else{
          v = rotate_right(v->l);
        }
      }else if(bf == -2){
        if(v->r->balanace_factor() == 1){
          rotate_right(v->r->l);
          v = rotate_left(v->r);
        }else{
          v = rotate_left(v->r);
        }
      }
      if(!v->p) return v;
      v = v->p;
    }
  }
  node *root;
  static constexpr int INIT_SIZE = 64; // 初期化時の1ブロック当たりのビット
  dynamic_bitvector(): root(nullptr){}
  dynamic_bitvector(const std::vector<bool> &bits){
    if(bits.empty()){
      root = nullptr;
      return;
    }
    std::vector<node*> nodes;
    int n = bits.size();
    __uint128_t b = 0;
    for(int i = 0; i < n; i++){
      b += (__uint128_t)bits[i] << (i % INIT_SIZE);
      if(i % INIT_SIZE == (INIT_SIZE - 1) || i == n - 1){
        node *r = new node(b, INIT_SIZE);
        r->sz = r->szsum = (i % INIT_SIZE) + 1;
        r->pop = r->popsum = r->bit.pop_count();
        nodes.push_back(r);
        b = 0;
      }
    }
    root = build(nodes, 0, nodes.size());
  }
  int size(){
    return size(root);
  }
  int pop(){
    return pop(root);
  }
  int rank1(int k){
    if(!root) return 0;
    node *v = root;
    int res = 0;
    while(v){
      int lsz = size(v->l);
      if(lsz <= k){
        int msz = lsz + v->sz;
        res += pop(v->l);
        if(k < msz) return res + v->bit.rank(k - lsz);
        res += v->pop;
        k -= msz;
        v = v->r;
      }else v = v->l;
    }
    return res;
  }
  int rank0(int k){
    return k - rank1(k);
  }
  int select1(int k){
    assert(k < pop());
    node *v = root;
    int lsz = 0;
    while(true){
      int lpop = pop(v->l);
      if(lpop <= k){
        int mpop = lpop + v->pop;
        if(k < mpop) return lsz + size(v->l) + v->bit.select1(k - lpop);
        k -= mpop;
        lsz += size(v->l) + v->sz;
        v = v->r;
      }else v = v->l;
    }
  }
  int select0(int k){
    assert(k < size() - pop());
    node *v = root;
    int lsz = 0;
    while(true){
      int lunpop = size(v->l) - pop(v->l);
      if(lunpop <= k){
        int munpop = lunpop + (v->sz - v->pop);
        if(k < munpop) return lsz + size(v->l) + v->bit.select0(k - lunpop);
        k -= munpop;
        lsz += size(v->l) + v->sz;
        v = v->r;
      }else v = v->l;
    }
  }
  void set(int k, bool f){
    assert(0 <= k && k < size());
    node *v = find_kth(root, k);
    v->bit.set(k, f);
    v->pop = v->bit.pop_count();
    update_rec(v);
  }
  bool get(int k){
    assert(0 <= k && k < size());
    node *v = find_kth(root, k);
    return v->bit.get(k);
  }
  void insert(int k, bool f){
    assert(k <= size());
    if(!root){
      root = new node(f, 1);
      return;
    }
    node *v = (k == size() ? right_most(root) : find_kth(root, k));
    if(k == size()) k = v->sz;
    if(v->sz < 127){
      v->bit.insert(k, f);
      v->sz++;
      if(f) v->pop++;
      update_rec(v);
    }else{
      v->bit.insert(k, f);
      std::pair<__uint128_t, __uint128_t> p = v->bit.split();
      node *u = new node(p.first, 64);
      v->bit.x = p.second;
      v->sz = 64;
      v->pop = v->bit.pop_count();
      if(!v->r){
        v->r = u;
        u->p = v;
      }else{
        node *q = left_most(v->r);
        q->l = u;
        u->p = q;
        update(q);
      }
      root = balancing(u);
    }
  }
  bool erase(int k){
    assert(0 <= k && k < size());
    node *v = find_kth(root, k);
    bool ret;
    ret = v->bit.erase(k);
    v->sz--;
    v->pop = v->bit.pop_count();
    update_rec(v);
    return ret;
  }
  bool erase_safe(int k){
    assert(0 <= k && k < size());
    node *v = find_kth(root, k);
    bool ret = v->bit.erase(k);
    v->sz--;
    v->pop = v->bit.pop_count();

    if(!v->sz){
      node *u = nullptr;
      if(v->l && v->r){
        u = left_most(v->r);
        node *up = u->p;
        u->l = v->l, u->l->p = u;
        if((u->p = v->p)){
          if(u->p->l == v) u->p->l = u;
          else u->p->r = u;
        }
        if(up != v){
          if((up->l = u->r)) up->l->p = up;
          u->r = v->r, u->r->p = u;
          root = balancing(up);
        }else root = balancing(u);
      }else if(v->l || v->r){
        u = v->l ? v->l : v->r;
        assert(!u->l && !u->r);
        if((u->p = v->p)){
          if(u->p->l == v) u->p->l = u;
          else u->p->r = u;
        }
        root = balancing(u);
      }else{ // vが葉
        assert(v->h = 1);
        if(!v->p) root = nullptr;
        else{
          u = v->p;
          if(u->l == v) u->l = nullptr;
          else u->r = nullptr;
          root = balancing(u);
        }
      }
    }else update_rec(v);
    return ret;
  }
};


#endif
