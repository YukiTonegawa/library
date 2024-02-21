#ifndef _DYNAMIC_STRING_PERSISTENT_H_
#define _DYNAMIC_STRING_PERSISTENT_H_
#include "../rolling_hash.hpp"

struct dynamic_string_persistent{
  using ull = unsigned long long;
  using rh = montgomery_rolling_hash;
  struct node{
    int h, sz;
    ull val, sum;
    node *l, *r;
    template<typename T>
    node(T _val): h(1), sz(1), val(rh::__mod(_val)), sum(val), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
private:
  static node *tmp_node;
  static node *copy_node(node *v){
    if(!v) return nullptr;
    return new node(*v);
  }
  static int __size(node *v){return v ? v->sz : 0;}
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz = 1;
    v->sum = v->val;
    if(v->l){
      v->sum = rh::add(v->l->sum, rh::mul(v->sum, rh::rpow[v->l->sz]));
      v->sz += v->l->sz;
    }
    if(v->r){
      v->sum = rh::add(v->sum, rh::mul(v->r->sum, rh::rpow[v->sz]));
      v->sz += v->r->sz;
    }
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
  static node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
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
  static node *set_inner(node *v, int k, ull x){
    int szl = v->l ? v->l->sz : 0;
    v = copy_node(v);
    if(k < szl) v->l = set_inner(v->l, k, x);
    else if(k > szl) v->r = set_inner(v->r, k - szl - 1, x);
    else v->val = x;
    update(v);
    return v;
  }
  static ull get_inner(node *v, int k){
    int szl = v->l ? v->l->sz : 0;
    if(k < szl) return get_inner(v->l, k);
    else if(k > szl) return get_inner(v->r, k - szl - 1);
    else return v->val;
  }
  static ull hash_inner(node *v, int l, int r){
    if(!v) return 0;
    if(l == 0 && r == v->sz) return v->sum;
    int szl = __size(v->l), szv = szl + 1;
    if(r <= szl) return hash_inner(v->l, l, r);
    if(szv <= l) return hash_inner(v->r, l - szv, r - szv);
    ull ql = hash_inner(v->l, l, szl), qr = hash_inner(v->r, 0, r - szv);
    ull sum = rh::add(ql, rh::mul(v->val, rh::rpow[szl - l]));
    sum = rh::add(sum, rh::mul(qr, rh::rpow[szl - l + 1]));
    return sum;
  }
  static node *insert_inner(node *v, int k, ull x){
    assert(__size(v) >= k);
    if(!v) return new node(x);
    int szl = v->l ? v->l->sz : 0;
    v = copy_node(v);
    if(k <= szl) v->l = insert_inner(v->l, k, x);
    else if(k > szl) v->r = insert_inner(v->r, k - szl - 1, x);
    update(v);
    return balance_insert(v);
  }
  static node *erase_inner(node *v, int k){
    assert(0 <= k && k < __size(v));
    int szl = v->l ? v->l->sz : 0;
    if(k < szl){
      v = copy_node(v);
      v->l = erase_inner(v->l, k);
    }
    else if(k > szl){
      v = copy_node(v);
      v->r = erase_inner(v->r, k - szl - 1);
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
  static int __lcp(node *v, const rolling_hash_string &b, int l1, int l2){
    if(!v || l2 >= b.size()) return 0;
    if(l1 == 0 && l2 + v->sz <= b.size() && v->sum == b.hash_range(l2, l2 + v->sz)) return v->sz;
    int szl = v->l ? v->l->sz : 0;
    if(szl + 1 <= l1) return __lcp(v->r, b, l1 - (szl + 1), l2);
    int left_cross = szl - l1, L = __lcp(v->l, b, l1, l2);
    if(L != left_cross) return L;
    l2 += left_cross;
    if(l2 >= b.size() || v->val != b.get(l2)) return L;
    return L + 1 + __lcp(v->r, b, 0, l2 + 1);
  }
  static int __lcp(node *v, const dynamic_rolling_hash_string &b, int l1, int l2, int &oksz, ull &ok){
    if(!v || l2 >= b.size()) return 0;
    if(l1 == 0 && l2 + v->sz <= b.size()){
      ull rh = rh::mul(rh::sub(b.hash_range(l2 + v->sz), ok), rh::rinvpow[l2]);
      if(v->sum == rh){
        ok = rh::add(ok, rh::mul(v->sum, rh::rpow[oksz]));
        oksz += v->sz;
        return v->sz;
      }
    }
    int szl = v->l ? v->l->sz : 0;
    if(szl + 1 <= l1) return __lcp(v->r, b, l1 - (szl + 1), l2, oksz, ok);
    int left_cross = szl - l1, L = __lcp(v->l, b, l1, l2, oksz, ok);
    if(L != left_cross) return L;
    l2 += left_cross;
    if(l2 >= b.size() || v->val != b.get(l2)) return L;
    ok = rh::add(ok, rh::mul(v->val, rh::rpow[oksz]));
    oksz++;
    return L + 1 + __lcp(v->r, b, 0, l2 + 1, oksz, ok);
  }
  static int __lcp_naive(node *v, node *u, int l1, int l2, int &oksz, ull &ok){
    if(!v || l2 >= __size(u)) return 0;
    if(l1 == 0 && l2 + v->sz <= __size(u)){
      ull rh = rh::mul(rh::sub(hash_inner(u, 0, l2 + v->sz), ok), rh::rinvpow[l2]);
      if(v->sum == rh){
        ok = rh::add(ok, rh::mul(v->sum, rh::rpow[oksz]));
        oksz += v->sz;
        return v->sz;
      }
    }
    int szl = v->l ? v->l->sz : 0;
    if(szl + 1 <= l1) return __lcp_naive(v->r, u, l1 - (szl + 1), l2, oksz, ok);
    int left_cross = szl - l1, L = __lcp_naive(v->l, u, l1, l2, oksz, ok);
    if(L != left_cross) return L;
    l2 += left_cross;
    if(l2 >= __size(u) || v->val != get_inner(u, l2)) return L;
    ok = rh::add(ok, rh::mul(v->val, rh::rpow[oksz]));
    oksz++;
    return L + 1 + __lcp_naive(v->r, u, 0, l2 + 1, oksz, ok);
  }
public:
  dynamic_string_persistent(){}

  template<typename T>
  static node *build(const std::vector<T> &v){
    if(v.empty()) return nullptr;
    int n = v.size();
    std::vector<node*> nodes(n);
    for(int i = 0; i < n; i++) nodes[i] = new node(v[i]);
    return build(nodes, 0, n);
  }
  static int size(node *v){
    return __size(v);
  }
  template<typename T>
  static node *set(node *v, int k, T x){
    assert(0 <= k && k < size(v));
    return set_inner(v, k, rh::__mod(x));
  }
  static ull get(node *v, int k){
    assert(0 <= k && k < size(v));
    return get_inner(v, k);
  }
  static node *insert(node *v, int k, ull x){
    assert(0 <= k && k <= size(v));
    return insert_inner(v, k, x);
  }
  static node *erase(node *v, int k){
    assert(0 <= k && k < size(v));
    return erase_inner(v, k);
  }
  static ull hash_range(node *v, int l, int r){
    assert(0 <= l && r <= size(v));
    if(l >= r) return 0;
    return hash_inner(v, l, r);
  }
  static ull hash_all(node *v){
    if(!v) return 0;
    return v->sum;
  }
  static int lcp(node *v, const rolling_hash_string &b, int l1, int l2){
    return __lcp(v, b, l1, l2);
  }
  static int lcp(const rolling_hash_string &a, node *v, int l1, int l2){
    return lcp(v, a, l2, l1);
  }
  static int lcp(node *v, const dynamic_rolling_hash_string &b, int l1, int l2){
    int oksz = l2;
    ull ok = b.hash_range(l2);
    return __lcp(v, b, l1, l2, oksz, ok);
  }
  static int lcp(const dynamic_rolling_hash_string &a, node *v, int l1, int l2){
    return lcp(v, a, l2, l1);
  }
  static int lcp(node *v, node *b, int l1, int l2){
    int oksz = l2;
    ull ok = hash_range(b, 0, l2);
    return __lcp_naive(v, b, l1, l2, oksz, ok);
  }
};
dynamic_string_persistent::node *dynamic_string_persistent::tmp_node = nullptr;
#endif