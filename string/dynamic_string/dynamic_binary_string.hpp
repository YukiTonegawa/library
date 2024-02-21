#ifndef _DYNAMIC_BINARY_STRING_H_
#define _DYNAMIC_BINARY_STRING_H_

#include "../../data_structure/bit_sequence/dynamic_bit.hpp"
#include "../rolling_hash.hpp"

struct dynamic_binary_string{
  using ull = unsigned long long;
  using rh = montgomery_rolling_hash;
  
  static std::vector<std::vector<ull>> htable; // 2^16未満のハッシュ
  
  static ull calc_hash(ull x, int sz){
    if(htable.empty()){
      htable.resize(1 << 16, std::vector<ull>(16));
      for(int i = 0; i < (1 << 16); i++){
        ull tmp = 0;
        for(int j = 0; j < 16; j++){
          ull y = rh::__mod('0' + ((i >> j) & 1));
          tmp = rh::add(tmp, rh::mul(y, rh::rpow[j]));
          htable[i][j] = tmp;
        }
      }
    }
    if(sz-- == 0) return 0;
    if(sz < 32){
      if(sz < 16){
        return htable[x & mask_0_16][sz];
      }else{
        return rh::add(htable[x & mask_0_16][15], rh::mul(htable[(x >> 16) & mask_0_16][sz - 16], rh::rpow[16]));
      }
    }else{
      ull a = rh::add(htable[x & mask_0_16][15], rh::mul(htable[(x >> 16) & mask_0_16][15], rh::rpow[16]));
      if(sz < 48){
        return rh::add(a, rh::mul(htable[(x >> 32) & mask_0_16][sz - 32], rh::rpow[32]));
      }else{
        ull b = rh::add(a, rh::mul(htable[(x >> 32) & mask_0_16][15], rh::rpow[32]));
        return rh::add(b, rh::mul(htable[(x >> 48) & mask_0_16][sz - 48], rh::rpow[48]));
      }
    }
  }
  static ull calc_hash(ull x, int l, int r){
    return calc_hash(x >> l, r - l);
  }
  struct node{
    int h, sz, pop, szsum, popsum;
    dynamic_bit64 bit;
    ull hval, hsum;
    node *l, *r;
    node(ull x, int sz): h(1), sz(sz), pop(__builtin_popcount(x)), szsum(sz), popsum(pop), bit(x), hval(calc_hash(x, sz)), hsum(hval), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  // updateはしない
  static std::pair<node*, node*> split_half(node *v){
    assert(v->sz == 64);
    auto [H, L] = v->bit.split_half();
    node *u = new node(H, 32);
    v->bit.x = L;
    v->sz = 32;
    v->pop = v->bit.popcount();
    v->hval = calc_hash(L, 32);
    return {v, u};
  }
  node *root;
  static node *tmp_node;
  static int size(node *v){return v ? v->szsum : 0;}
  static int pop(node *v){return v ? v->popsum : 0;}
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0, v->r ? v->r->h : 0) + 1;
    v->szsum = v->sz;
    v->popsum = v->pop;
    v->hsum = v->hval;
    if(v->l){
      v->szsum += v->l->szsum;
      v->popsum += v->l->popsum;
      v->hsum = rh::add(v->l->hsum, rh::mul(v->hval, rh::rpow[v->l->szsum]));
    }
    if(v->r){
      v->hsum = rh::add(v->hsum, rh::mul(v->r->hsum, rh::rpow[v->szsum]));
      v->szsum += v->r->szsum;
      v->popsum += v->r->popsum;
    }
  }
  static node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) / 2;
    node *v = nodes[m];
    if(l < m) v->l = build(nodes, l, m);
    if(m + 1 < r) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
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
  static node *balance(node *v){
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
  static int __rank1(node *v, int r){
    if(!v) return 0;
    int res = 0;
    while(v){
      int lsz = size(v->l);
      if(lsz <= r){
        int msz = lsz + v->sz;
        res += pop(v->l);
        if(r < msz) return res + v->bit.rank(r - lsz);
        res += v->pop;
        r -= msz;
        v = v->r;
      }else v = v->l;
    }
    return res;
  }
  static int __select1(node *v, int k){
    if(pop(v) <= k) return -1;
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
  static int __select0(node *v, int k){
    if(size(v) - pop(v) <= k) return -1;
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
  static void __set(node *v, int k, bool x){
    int szl = v->l ? v->l->sz : 0;
    int szv = szl + v->sz;
    if(k < szl) __set(v->l, k, x);
    else if(k >= szv) __set(v->r, k - szv, x);
    else{
      if(v->bit.get(x) == x) return;
      v->bit.set(k, x);
      if(x){
        v->pop++, v->popsum++;
        v->hval = rh::add(v->hval, rh::rpow[k - szl]);
        v->hsum = rh::add(v->hval, rh::rpow[k]);
      }else{
        v->pop--, v->popsum--;
        v->hval = rh::sub(v->hval, rh::rpow[k - szl]);
        v->hsum = rh::sub(v->hval, rh::rpow[k]);
      }
      return;
    }
    update(v);
  }
  static bool __get(node *v, int k){
    v = find_kth(v, k);
    return v->bit.get(k);
  }
  static node *__insert_leftmost(node *v, node *u){
    if(!v) return u;
    v->l = __insert_leftmost(v->l, u);
    update(v);
    return balance(v);
  }
  static node *__insert(node *v, int k, bool f){
    if(!v) return new node(f, 1);
    int szl = (v->l ? v->l->szsum : 0), szv = szl + v->sz;
    if(k < szl){
      v->l = __insert(v->l, k, f);
    }else if(k > szv){
      v->r = __insert(v->r, k - szv, f);
    }else{
      k -= szl;
      assert(0 <= k && k <= v->sz);
      v->bit.insert(k, f);
      v->sz++;
      v->hval = calc_hash(v->bit.x, v->sz);
      if(f) v->pop++;
      if(v->sz == 64){
        node *u;
        std::tie(v, u) = split_half(v);
        v->r = __insert_leftmost(v->r, u);
      }
    }
    update(v);
    return balance(v);
  }
  static node *__cut_leftmost(node *v){
    if(v->l){
      v->l = __cut_leftmost(v->l);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  static node *__erase(node *v, int k){
    if(!v) return v;
    int szl = (v->l ? v->l->szsum : 0);
    int szv = szl + v->sz;
    if(k < szl){
      v->l = __erase(v->l, k);
    }else if(k >= szv){
      v->r = __erase(v->r, k - szv);
    }else{
      k -= szl;
      v->bit.erase(k);
      v->sz--;
      v->hval = calc_hash(v->bit.x, v->sz);
      v->pop = v->bit.popcount();
      if(v->sz == 0){
        if(!v->r || !v->l) return !v->l ? v->r : v->l;
        node *u = __cut_leftmost(v->r);
        tmp_node->l = v->l;
        tmp_node->r = u;
        v = tmp_node;
      }
    }
    update(v);
    return balance(v);
  }
  static ull __hash_range(node *v, int r){
    if(!v || r <= 0) return 0;
    if(r == v->szsum) return v->hsum;
    int szl = (v->l ? v->l->szsum : 0), szv = szl + v->sz;
    if(r <= szl) return __hash_range(v->l, r);
    ull ql = (v->l ? v->l->hsum : 0);
    if(szv < r){
      ull res = rh::add(ql, rh::mul(v->hval, rh::rpow[szl]));
      return rh::add(res, rh::mul(__hash_range(v->r, r - szv), rh::rpow[szv]));
    }
    return rh::add(ql, rh::mul(calc_hash(v->bit.x, 0, r - szl), rh::rpow[szl]));
  }
  static ull __hash_range(node *v, int l, int r){
    if(!v || l >= r) return 0;
    if(l == 0 && r == v->szsum) return v->hsum;
    int szl = (v->l ? v->l->szsum : 0), szv = szl + v->sz;
    if(r <= szl) return __hash_range(v->l, l, r);
    if(szv <= l) return __hash_range(v->r, l - szv, r - szv);
    ull ql = __hash_range(v->l, l, szl), qr = __hash_range(v->r, 0, r - szv);
    
    int llen = std::max(szl - l, 0), rlen = szv - l;
    l = std::max(l, szl) - szl, r = std::min(r, szv) - szl; // vの[l, r)
    ull res = rh::add(ql, rh::mul(calc_hash(v->bit.x, l, r), rh::rpow[llen]));
    return rh::add(res, rh::mul(qr, rh::rpow[rlen]));
  }
  // [l, r)を列挙(r - l <= 64)
  static ull __enumerate64(node *v, int l, int r){
    if(!v || l >= r) return 0;
    if(l == 0 && r == v->szsum) return v->bit.x;
    int szl = (v->l ? v->l->szsum : 0), szv = szl + v->sz;
    if(r <= szl) return __enumerate64(v->l, l, r);
    if(szv <= l) return __enumerate64(v->r, l - szv, r - szv);
    ull L = __enumerate64(v->l, l, szl), R = __enumerate64(v->r, 0, r - szv);
    int llen = std::max(szl - l, 0), rlen = szv - l;
    l = std::max(l, szl) - szl, r = std::min(r, szv) - szl; // vの[l, r)
    ull M = (r == 64 ? v->bit.x : v->bit.x & ((1ULL << r) - 1)) >> l;
    return L | (M << llen) | (R << rlen);
  }
  static int __lcp(node *v, dynamic_binary_string &b, int l1, int l2, int &oksz, ull &ok){
    if(!v || l1 >= v->szsum || l2 >= b.size()) return 0;
    if(l1 == 0 && l2 + v->szsum <= b.size()){
      ull h = rh::mul(rh::sub(b.hash_range(l2 + v->szsum), ok), rh::rinvpow[l2]);
      if(v->hsum == h){
        ok = rh::add(ok, rh::mul(v->hsum, rh::rpow[oksz]));
        oksz += v->szsum;
        return v->szsum;
      }
    }
    int szl = (v->l ? v->l->szsum : 0), szv = szl + v->sz;
    if(szv <= l1) return __lcp(v->r, b, l1 - szv, l2, oksz, ok);
    int left_cross = std::max(0, szl - l1), L = __lcp(v->l, b, l1, l2, oksz, ok);
    l2 += left_cross;
    if(L != left_cross || l2 == b.size()) return L;
    l1 = std::max(l1, szl) - szl;
    int mid_cross = v->sz - l1;
    ull vx = (l1 == 0 ? v->hval : calc_hash(v->bit.x, l1, v->sz));
    if(l2 + mid_cross > b.size() || vx != b.hash_range(l2, l2 + mid_cross)){
      int len = std::min(mid_cross, b.size() - l2);
      ull x = v->bit.x >> l1;
      ull y = dynamic_binary_string::__enumerate64(b.root, l2, l2 + len);
      return L + std::min(len, x == y ? 64 : select_64bit(x ^ y, 0));
    }
    ok = rh::add(ok, rh::mul(vx, rh::rpow[oksz]));
    oksz += mid_cross;
    return L + mid_cross + __lcp(v->r, b, 0, l2 + mid_cross, oksz, ok);
  }
public:
  static constexpr int INIT_SIZE = 32; // 初期化時の1ブロック当たりのビット
  dynamic_binary_string(): root(nullptr){}
  dynamic_binary_string(const std::vector<bool> &bits){
    if(bits.empty()){
      root = nullptr;
      return;
    }
    std::vector<node*> nodes;
    int n = bits.size();
    ull b = 0;
    for(int i = 0; i < n; i++){
      int j = i % INIT_SIZE;
      b += (ull)bits[i] << j;
      if(j == (INIT_SIZE - 1) || i == n - 1){
        nodes.push_back(new node(b, j + 1));
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
  // [0, r)の1の数
  int rank1(int r){
    assert(0 <= r && r <= size());
    return __rank1(root, r);
  }
  // [0, r)の0の数
  int rank0(int r){
    return r - rank1(r);
  }
  // k番目の1の位置, ない場合は-1
  int select1(int k){
    return __select1(root, k);
  }
  // k番目の0の位置, ない場合は-1
  int select0(int k){
    return __select0(root, k);
  }
  // k以降の最も左の1(k含む), ない場合は-1
  int find_next1(int k){
    return select1(rank1(k));
  }
  // k以降の最も左の0(k含む), ない場合は-1
  int find_next0(int k){
    return select0(rank0(k));
  }
  // k以前の最も右の1(k含む), ない場合は-1
  int find_prev1(int k){
    int cnt = rank1(k + 1);
    if(cnt == 0) return -1;
    return select1(cnt - 1);
  }
  // k以前の最も右の0(k含む), ない場合は-1
  int find_prev0(int k){
    int cnt = rank0(k + 1);
    if(cnt == 0) return -1;
    return select0(cnt - 1);
  }
  void set(int k, bool x){
    assert(k < size());
    __set(root, k, x);
  }
  bool get(int k){
    assert(k < size());
    return __get(root, k);
  }
  void insert(int k, bool x){
    assert(0 <= k && k <= size());
    root = __insert(root, k, x);
  }
  void erase(int k){
    assert(0 <= k && k < size());
    root = __erase(root, k);
  }
  ull hash_range(int r){
    assert(0 <= r && r <= size());
    return __hash_range(root, r);
  }
  ull hash_range(int l, int r){
    assert(0 <= l && l <= r && r <= size());
    return __hash_range(root, l, r);
  }
  ull hash_all(){
    return root ? root->hsum : 0;
  }
  static int lcp(dynamic_binary_string &a, dynamic_binary_string &b, int l1, int l2){
    int oksz = l2;
    ull ok = b.hash_range(l2);
    return __lcp(a.root, b, l1, l2, oksz, ok);
  }
};
std::vector<std::vector<unsigned long long>> dynamic_binary_string::htable;
dynamic_binary_string::node *dynamic_binary_string::tmp_node = nullptr;
#endif