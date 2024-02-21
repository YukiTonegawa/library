#ifndef _DYNAMIC_STRING_MERGE_SPLIT_PERSISTENT_H_
#define _DYNAMIC_STRING_MERGE_SPLIT_PERSISTENT_H_

#include "../rolling_hash.hpp"
struct dynamic_string_merge_split_persistent{
  using ull = unsigned long long;
  using rh = montgomery_rolling_hash;
  struct node{
    node *l, *r;
    bool red;
    int ra, sz;
    ull val;
    // 葉
    node(ull val): l(nullptr), r(nullptr), red(false), ra(0), sz(1), val(val){}
    // 中間ノード
    node(node *l, node *r, bool red): l(l), r(r), red(red), ra(l->ra + !(l->red)), sz(l->sz + r->sz){
      val = rh::add(l->val, rh::mul(r->val, rh::rpow[l->sz]));
    }
  };
private:
  using pn = std::pair<node*, node*>;
  static node *merge_sub(node *a, node *b){
    if(a->ra < b->ra){
      node *c = merge_sub(a, b->l);
      if(!b->red && c->red && c->l->red){
        if(!b->r->red){
          return new node(c->l, new node(c->r, b->r, 1), 0);
        }else{
          b->r->red = 0;
          c->red = 0;
          return new node(c, b->r, 1);
        }
      }else{
        return new node(c, b->r, b->red);
      }
    }else if(a->ra > b->ra){
      node *c = merge_sub(a->r, b);
      if(!a->red && c->red && c->r->red){
        if(!a->l->red){
          return new node(new node(a->l, c->l, 1), c->r, 0);
        }else{
          a->l->red = 0;
          c->red = 0;
          return new node(a->l, c, 1);
        }
      }else{
        return new node(a->l, c, a->red);
      }
    }else{
      return new node(a, b, 1);
    }
  }
  static node *__merge(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    node *c = merge_sub(a, b);
    c->red = 0;
    return c;
  }
  static node *as_root(node *a){
    if(!a) return nullptr;
    if(a->red) a = new node(a->l, a->r, 0);
    return a;
  }
  static pn __split(node *a, int k){
    int sz = a->sz, szl = (a->l ? a->l->sz : 0);
    if(k == 0 || k == sz) return (!k ? pn{nullptr, a} : pn{a, nullptr});
    pn res;
    if(k < szl){
      auto [l, r] = __split(a->l, k);
      res = pn{l, __merge(r, as_root(a->r))};
    }else if(k > szl){
      auto [l, r] = __split(a->r, k - szl);
      res = pn{__merge(as_root(a->l), l), r};
    }else{
      res = pn{as_root(a->l), as_root(a->r)};
    }
    return res;
  }
  static node *__set(node *a, int k, ull x){
    a = new node(*a);
    if(!a->l && !a->r){
      assert(k == 0);
      a->val = x;
      return a;
    }
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) a->l = __set(a->l, k, x);
    else a->r = __set(a->r, k - szl, x);
    a->val = rh::add(a->l->val, rh::mul(a->r->val, rh::rpow[szl]));
    return a;
  }
  static ull __get(node *a, int k){
    if(!a->l && !a->r){
      assert(k == 0);
      return a->val;
    }
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) return __get(a->l, k);
    else return __get(a->r, k - szl);
  }
  static ull query(node *a, int l, int r){
    if(!a || l >= r || a->sz <= l || r <= 0) return 0;
    if(l <= 0 && a->sz <= r) return a->val;
    if(!a->l && !a->r) return a->val;
    int szl = a->l->sz;
    if(r <= szl) return query(a->l, l, r);
    if(szl <= l) return query(a->r, l - szl, r - szl);
    return rh::add(query(a->l, l, szl), rh::mul(query(a->r, 0, r - szl), rh::rpow[szl]));
  }
  template<typename T>
  static node *build(const std::vector<T> &v, int l, int r){
    if(l == r) return nullptr;
    if(r - l == 1) return new node(rh::__mod(v[l]));
    int mid = (l + r) / 2;
    node *L = build(v, l, mid);
    node *R = build(v, mid, r);
    return merge(L, R);
  }
  static int __lcp(node *a, const rolling_hash_string &b, int l1, int l2){
    if(!a || l2 >= b.size()) return 0;
    if(l1 == 0 && l2 + a->sz <= b.size() && a->val == b.hash_range(l2, l2 + a->sz)) return a->sz;
    int szl = a->l ? a->l->sz : 0;
    if(szl <= l1) return __lcp(a->r, b, l1 - szl, l2);
    int left_cross = szl - l1, L = __lcp(a->l, b, l1, l2);
    if(L != left_cross) return L;
    l2 += left_cross;
    return L + __lcp(a->r, b, 0, l2);
  }
  static int __lcp(node *a, const dynamic_rolling_hash_string &b, int l1, int l2, int &oksz, ull &ok){
    if(!a || l2 >= b.size()) return 0;
    if(l1 == 0 && l2 + a->sz <= b.size()){
      ull rh = rh::mul(rh::sub(b.hash_range(l2 + a->sz), ok), rh::rinvpow[l2]);
      if(a->val == rh){
        ok = rh::add(ok, rh::mul(a->val, rh::rpow[oksz]));
        oksz += a->sz;
        return a->sz;
      }
    }
    int szl = a->l ? a->l->sz : 0;
    if(szl <= l1) return __lcp(a->r, b, l1 - szl, l2, oksz, ok);
    int left_cross = szl - l1, L = __lcp(a->l, b, l1, l2, oksz, ok);
    if(L != left_cross) return L;
    l2 += left_cross;
    return L + __lcp(a->r, b, 0, l2, oksz, ok);
  }
  static int __lcp(node *v, node *b, int l1, int l2, int &oksz, ull &ok){
    if(!v || l2 >= size(b)) return 0;
    if(l1 == 0 && l2 + v->sz <= size(b)){
      ull rh = rh::mul(rh::sub(hash_range(b, 0, l2 + v->sz), ok), rh::rinvpow[l2]);
      if(v->val == rh){
        ok = rh::add(ok, rh::mul(v->val, rh::rpow[oksz]));
        oksz += v->sz;
        return v->sz;
      }
    }
    int szl = v->l ? v->l->sz : 0;
    if(szl <= l1) return __lcp(v->r, b, l1 - szl, l2, oksz, ok);
    int left_cross = szl - l1, L = __lcp(v->l, b, l1, l2, oksz, ok);
    if(L != left_cross) return L;
    l2 += left_cross;
    return L + __lcp(v->r, b, 0, l2, oksz, ok);
  }
public:
  dynamic_string_merge_split_persistent(){}

  template<typename T>
  static node *build(const std::vector<T> &v){
    return build(v, 0, v.size());
  }
  static int size(node *a){
    return a ? a->sz : 0;
  }
  template<typename T>
  static node *set(node *a, int k, T x){
    assert(k < size(a));
    return __set(a, k, rh::__mod(x));
  }
  static ull get(node *a, int k){
    assert(k < size(a));
    return __get(a, k);
  }
  // k番目にxを挿入
  template<typename T>
  static node *insert(node *a, int k, T x){
    auto [A, B] = __split(a, k);
    return __merge(A, __merge(new node(rh::__mod(x)), B));
  }
  static node *erase(node *a, int k){
    assert(size(a) > k);
    auto [A, B] = __split(a, k);
    assert(B);
    auto [B2, C] = __split(B, 1);
    return __merge(A, C);
  }
  static ull hash_range(node *a, int l, int r){
    assert(0 <= l && r <= size(a));
    return query(a, l, r);
  }
  // 2つに分割
  static std::pair<node*, node*> split(node *a, int k){
    assert(k <= size(a));
    return __split(a, k);
  }
  // a, bをマージ
  static node *merge(node *a, node *b){
    return __merge(a, b);
  }
  static int lcp(node *a, const rolling_hash_string &b, int l1, int l2){
    return __lcp(a, b, l1, l2);
  }
  static int lcp(const rolling_hash_string &a, node *b, int l1, int l2){
    return lcp(b, a, l2, l1);
  }
  static int lcp(node *a, const dynamic_rolling_hash_string &b, int l1, int l2){
    int oksz = l2;
    ull ok = b.hash_range(l2);
    return __lcp(a, b, l1, l2, oksz, ok);
  }
  static int lcp(const dynamic_rolling_hash_string &a, node *b, int l1, int l2){
    return lcp(b, a, l2, l1);
  }
  static int lcp(node *a, node *b, int l1, int l2){
    int oksz = l2;
    ull ok = hash_range(b, 0, l2);
    return __lcp(a, b, l1, l2, oksz, ok);
  }
};
#endif