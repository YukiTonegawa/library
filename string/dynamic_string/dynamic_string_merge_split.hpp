#ifndef _DYNAMIC_STRING_MERGE_SPLIT_H_
#define _DYNAMIC_STRING_MERGE_SPLIT_H_
#include "../rolling_hash.hpp"

struct dynamic_string_merge_split{
  using ull = unsigned long long;
  using rh = montgomery_rolling_hash;
private:
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
  static std::vector<node*> stock;
  static node *reuse(node *a, node *l, node *r, bool red){
    a->l = l, a->r = r, a->red = red;
    a->ra = l->ra + !(l->red);
    a->sz = l->sz + r->sz;
    a->val = rh::add(a->l->val, rh::mul(a->r->val, rh::rpow[a->l->sz]));
    return a;
  }
  node *root;
  using pn = std::pair<node*, node*>;
  static node *merge_sub(node *a, node *b){
    if(a->ra < b->ra){
      node *c = merge_sub(a, b->l);
      if(!b->red && c->red && c->l->red){
        if(!b->r->red){
          return reuse(c, c->l, reuse(b, c->r, b->r, 1), 0);
        }else{
          b->r->red = 0;
          c->red = 0;
          return reuse(b, c, b->r, 1);
        }
      }else{
        return reuse(b, c, b->r, b->red);
      }
    }else if(a->ra > b->ra){
      node *c = merge_sub(a->r, b);
      if(!a->red && c->red && c->r->red){
        if(!a->l->red){
          return reuse(c, reuse(a, a->l, c->l, 1), c->r, 0);
        }else{
          a->l->red = 0;
          c->red = 0;
          return reuse(a, a->l, c, 1);
        }
      }else{
        return reuse(a, a->l, c, a->red);
      }
    }else{
      if(stock.empty()) return new node(a, b, 1);
      node *u = stock.back();
      stock.pop_back();
      return reuse(u, a, b, 1);
    }
  }
  static node *merge(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    node *c = merge_sub(a, b);
    c->red = 0;
    return c;
  }
  static node *as_root(node *a){
    if(!a) return nullptr;
    a->red = 0;
    return a;
  }
  static pn split(node *a, int k){
    int sz = a->sz, szl = (a->l ? a->l->sz : 0);
    if(k == 0 || k == sz) return (!k ? pn{nullptr, a} : pn{a, nullptr});
    pn res;
    if(k < szl){
      auto [l, r] = split(a->l, k);
      res = pn{l, merge(r, as_root(a->r))};
    }else if(k > szl){
      auto [l, r] = split(a->r, k - szl);
      res = pn{merge(as_root(a->l), l), r};
    }else{
      res = pn{as_root(a->l), as_root(a->r)};
    }
    if(a) stock.push_back(a);
    return res;
  }
  void set(node *a, int k, ull x){
    if(!a->l && !a->r){
      assert(k == 0);
      a->val = x;
      return;
    }
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) set(a->l, k, x);
    else set(a->r, k - szl, x);
    a->val = rh::add(a->l->val, rh::mul(a->r->val, rh::rpow[szl]));
  }
  ull get(node *a, int k){
    if(!a->l && !a->r){
      assert(k == 0);
      return a->val;
    }
    int szl = a->l ? a->l->sz : 0;
    if(k < szl) return get(a->l, k);
    else return get(a->r, k - szl);
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
  node *build(const std::vector<T> &v, int l, int r){
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
  static int __lcp(node *v, const dynamic_string_merge_split &b, int l1, int l2, int &oksz, ull &ok){
    if(!v || l2 >= b.size()) return 0;
    if(l1 == 0 && l2 + v->sz <= b.size()){
      ull rh = rh::mul(rh::sub(b.hash_range(0, l2 + v->sz), ok), rh::rinvpow[l2]);
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
  dynamic_string_merge_split(node *a): root(a){}
public:
  dynamic_string_merge_split(): root(nullptr){}
  template<typename T>
  dynamic_string_merge_split(const std::vector<T> &v): root(build(v, 0, v.size())){}

  int size()const{
    return root ? root->sz : 0;
  }
  template<typename T>
  void set(int k, T x){
    assert(k < size());
    set(root, k, rh::__mod(x));
  }
  ull get(int k){
    assert(k < size());
    return get(root, k);
  }
  // k番目にxを挿入
  template<typename T>
  void insert(int k, T x){
    auto [a, b] = split(root, k);
    root = merge(a, merge(new node(rh::__mod(x)), b));
  }
  void erase(int k){
    assert(root->sz > k);
    auto [a, b] = split(root, k);
    assert(b);
    auto [b2, c] = split(b, 1);
    root = merge(a, c);
    if(b2) stock.push_back(b2);
  }
  ull hash_range(int l, int r)const{
    assert(0 <= l && r <= size());
    return query(root, l, r);
  }
  std::pair<dynamic_string_merge_split, dynamic_string_merge_split> split(int k){
    return split(*this, k);
  }
  // 2つに分割. 永続でないためこれ自身のaのrootはnullptrになる
  static std::pair<dynamic_string_merge_split, dynamic_string_merge_split> split(dynamic_string_merge_split &a, int k){
    assert(k <= a.size());
    auto [l, r] = split(a.root, k);
    a.root = nullptr;
    return {dynamic_string_merge_split(l), dynamic_string_merge_split(r)};
  }
  // a, bをマージ. 永続でないためa, bのrootはnullptrになる
  static dynamic_string_merge_split merge(dynamic_string_merge_split &a, dynamic_string_merge_split &b){
    dynamic_string_merge_split res(merge(a.root, b.root));
    a.root = b.root = nullptr;
    return res;
  }
  static int lcp(const dynamic_string_merge_split &a, const rolling_hash_string &b, int l1, int l2){
    return __lcp(a.root, b, l1, l2);
  }
  static int lcp(const rolling_hash_string &a, const dynamic_string_merge_split &b, int l1, int l2){
    return lcp(b, a, l2, l1);
  }
  static int lcp(const dynamic_string_merge_split &a, const dynamic_rolling_hash_string &b, int l1, int l2){
    int oksz = l2;
    ull ok = b.hash_range(l2);
    return __lcp(a.root, b, l1, l2, oksz, ok);
  }
  static int lcp(const dynamic_rolling_hash_string &a, const dynamic_string_merge_split &b, int l1, int l2){
    return lcp(b, a, l2, l1);
  }
  static int lcp(const dynamic_string_merge_split &a, const dynamic_string_merge_split &b, int l1, int l2){
    int oksz = l2;
    ull ok = b.hash_range(0, l2);
    return __lcp(a.root, b, l1, l2, oksz, ok);
  }
};
std::vector<dynamic_string_merge_split::node*> dynamic_string_merge_split::stock;
#endif