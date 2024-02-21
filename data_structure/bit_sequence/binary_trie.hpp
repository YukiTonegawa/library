#ifndef _BINARY_TRIE_H_
#define _BINARY_TRIE_H_
#include <vector>
#include <cassert>
#include "../../algebraic_structure/monoid.hpp"

// @param 0 < bitlen < (Idxのbit長)
template<typename Idx, int bitlen>
struct binary_trie{
  static constexpr Idx inf = Idx(1) << bitlen;
  static constexpr int __lcp(Idx a, Idx b){
    return (a == b ? 0 : __builtin_clzll(a ^ b) - 64) + bitlen;
  }
private:
  struct node{
    Idx idx;
    node *l, *r;
    int sz, szsum;
    node(Idx idx): idx(idx), l(nullptr), r(nullptr), sz(1), szsum(1){}
  };
  node *make_node(Idx idx){
    return new node(idx);
  }
  void eval(node *v){
    v->szsum = v->sz + (v->l ? v->l->szsum : 0) + (v->r ? v->r->szsum : 0);
  }
  void insert(node* &v, Idx k, int x, Idx l, Idx r, bool unique){
    if(!v){
      v = make_node(k);
      return;
    }
    if(v->idx == k){
      if(unique) return;
      v->sz += x;
      eval(v);
      return;
    }
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < mid){
      if(v->idx < k) std::swap(v->idx, k), std::swap(v->sz, x);
      insert(v->l, k, x, l, mid, unique);
    }else{
      if(v->idx > k) std::swap(v->idx, k), std::swap(v->sz, x);
      insert(v->r, k, x, mid, r, unique);
    }
    eval(v);
  }
   void erase(node* &v, Idx k, Idx l, Idx r, bool all){
    if(!v) return;
    if(v->idx == k){
      v->sz = (all ? 0 : std::max(0, v->sz - 1));
      eval(v);
      return;
    }
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < mid) erase(v->l, k, l, mid, all);
    else erase(v->r, k, mid, r, all);
    eval(v);
  }
  int find(node *v, Idx k, Idx l, Idx r){
    if(!v) return 0;
    if(v->idx == k) return v->sz;
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < mid) return find(v->l, k, l, mid);
    else return find(v->r, k, mid, r);
  }
  int count_range(node *v, Idx a, Idx b, Idx l, Idx r){
    if(!v || b <= l || r <= a) return 0;
    if(a <= l && r <= b) return v->szsum;
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    int ret = count_range(v->l, a, b, l, mid);
    if(a <= v->idx && v->idx < b) ret += v->sz;
    return ret + count_range(v->r, a, b, mid, r);
  }
  void count_lcp2(node *v, Idx k, Idx l, Idx r, std::vector<int> &ret){
    if(!v) return;
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < mid){
      count_lcp2(v->l, k, l, mid, ret);
      ret[__lcp(v->idx, k)] += v->sz;
      if(v->r) ret[__lcp(v->r->idx, k)] += v->r->szsum;
    }else{
      if(v->l) ret[__lcp(v->l->idx, k)] += v->l->szsum;
      ret[__lcp(v->idx, k)] += v->sz;
      count_lcp2(v->r, k, mid, r, ret);
    }
  }
  Idx kth_smallest(node *v, int k, Idx l, Idx r){
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    int lsz = (v->l ? v->l->szsum : 0);
    if(k < lsz){
      return kth_smallest(v->l, k, l, mid);
    }else if(k < lsz + v->sz){
      return v->idx;
    }else{
      return kth_smallest(v->r, k - lsz - v->sz, mid, r);
    }
  }
  node *root;
public:
  binary_trie(): root(nullptr){static_assert(bitlen < sizeof(Idx) * 8);}
  int size(){
    return root ? root->szsum : 0;
  }
  int find(Idx k){
    return find(root, k, 0, inf);
  }
  //　すでに同じ要素があり, unique = trueの場合何もしない
  void insert(Idx k, bool unique = false){
    insert(root, k, 1, 0, inf, unique);
  }
  //　all = trueの場合全て, そうでない場合1個だけ消す, 要素kがない場合は何もしない
  void erase(Idx k, bool all = false){
    erase(root, k, 0, inf, all);
  }
  int count_range(Idx l, Idx r){
    return count_range(root, l, r, 0, inf);
  }
  // kと接頭辞(bitlen-1, bitlen-2...)がt個以上一致する要素の数
  int count_lcp_ge(Idx k, int t){
    if(t > bitlen) return 0;
    t = std::max(t, 0);
    int low = bitlen - t;
    Idx p = (~((Idx(1) << low) - 1)) & k;
    return count_range(p, p + (Idx(1) << low));
  }
  // kと接頭辞(bitlen-1, bitlen-2...)がちょうどt個一致する要素の数
  int count_lcp(Idx k, int t){
    if(t > bitlen || t < 0) return 0;
    t = std::max(t, 0);
    int low = bitlen - t;
    Idx p = (~((Idx(1) << low) - 1)) & k;
    Idx l = p, r = p + (Idx(1) << low);
    if(t == bitlen) return count_range(l, r);
    if((k >> (low - 1)) & 1) r -= Idx(1) << (low - 1);
    else l += Idx(1) << (low - 1);
    return count_range(l, r);
  }
  // ans[t] := kと接頭辞(bitlen-1, bitlen-2...)がちょうどt個一致する要素のsum
  std::vector<int> count_lcp2(Idx k){
    std::vector<int> ret(bitlen + 1, 0); // [0, bitlen]
    cpunt_lcp2(root, k, 0, inf, ret);
    return ret;
  }
  // 最小要素, ない場合はinf
  Idx min(){
    return size() ? kth_smallest(0) : inf;
  }
  // 最大要素, ない場合はinf
  Idx max(){
    return size() ? kth_smallest(size() - 1) : inf;
  }
  // k番目に小さい要素(同じ要素を重複して数える), ない場合はinf
  Idx kth_smallest(int k){
    return k < size() ? kth_smallest(root, k, 0, inf) : inf;
  }
  // x以上の最小要素, ない場合はinf
  Idx lower_bound(Idx x){
    assert(0 <= x && x < inf);
    return kth_smallest(count_range(0, x));
  }
  // x以下の最大要素, ない場合はinf
  Idx lower_bound_rev(Idx x){
    assert(0 <= x && x < inf);
    int cnt = count_range(0, x + 1);
    return cnt ? kth_smallest(cnt - 1) : inf;
  }
};


// @param 0 < bitlen < (Idxのbit長)
template<typename Idx, int bitlen, typename monoid>
struct binary_trie_monoid{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  static constexpr Idx inf = Idx(1) << bitlen;
  static constexpr int __lcp(Idx a, Idx b){
    return (a == b ? 0 : __builtin_clzll(a ^ b) - 64) + bitlen;
  }
private:
  struct node{
    Idx idx;
    node *l, *r;
    Val val, sum;
    node(Idx idx, Val val): idx(idx), l(nullptr), r(nullptr), val(val), sum(val){}
  };
  node *make_node(Idx idx,  Val val = id()){
    return new node(idx, val);
  }
  void eval(node *v){
    v->sum = v->val;
    if(v->l) v->sum = merge(v->l->sum, v->sum);
    if(v->r) v->sum = merge(v->sum, v->r->sum);
  }
  void set(node* &v, Idx k, Val x, Idx l, Idx r){
    if(!v){
      v = make_node(k, x);
      return;
    }
    if(v->idx == k){
      v->val = x;
      eval(v);
      return;
    }
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < mid){
      if(v->idx < k) std::swap(v->idx, k), std::swap(v->val, x);
      set(v->l, k, x, l, mid);
    }else{
      if(v->idx > k) std::swap(v->idx, k), std::swap(v->val, x);
      set(v->r, k, x, mid, r);
    }
    eval(v);
  }
  Val get(node *v, Idx k, Idx l, Idx r){
    if(!v) return id();
    if(v->idx == k) return v->val;
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < mid) return get(v->l, k, l, mid);
    else return get(v->r, k, mid, r);
  }
  void update(node* &v, Idx k, Val x, Idx l, Idx r){
    if(!v){
      v = make_node(k, x);
      return;
    }
    if(v->idx == k){
      v->val = merge(v->val, x);
      eval(v);
      return;
    }
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < mid){
      if(v->idx < k) std::swap(v->idx, k), std::swap(v->val, x);
      update(v->l, k, x, l, mid);
    }else{
      if(v->idx > k) std::swap(v->idx, k), std::swap(v->val, x);
      update(v->r, k, x, mid, r);
    }
    eval(v);
  }
  Val query(node *v, Idx a, Idx b, Idx l, Idx r){
    if(!v || b <= l || r <= a) return id();
    if(a <= l && r <= b) return v->sum;
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    Val ret = query(v->l, a, b, l, mid);
    if(a <= v->idx && v->idx < b) ret = merge(ret, v->val);
    return merge(ret, query(v->r, a, b, mid, r));
  }
  void query_lcp2(node *v, Idx k, Idx l, Idx r, std::vector<Val> &ret){
    if(!v) return;
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < mid){
      query_lcp2(v->l, k, l, mid, ret);
      int tmp = __lcp(v->idx, k);
      ret[tmp] = merge(ret[tmp], v->val);
      if(v->r){
        tmp = __lcp(v->r->idx, k);
        ret[tmp] = merge(ret[tmp], v->r->sum);
      }
    }else{
      int tmp;
      if(v->l){
        tmp = __lcp(v->l->idx, k);
        ret[tmp] = merge(ret[tmp], v->l->sum);
      }
      tmp = __lcp(v->idx, k);
      ret[tmp] = merge(ret[tmp], v->val);
      query_lcp2(v->r, k, mid, r, ret);
    }
  }
  template<typename F>
  Idx bisect_from_left(node *v, const Idx k, F &f, Idx l, Idx r, Val &ok){
    if(!v || r <= k) return inf;
    Val m = merge(ok, v->sum);
    if(k <= l && !f(m)){
      ok = m;
      return inf;
    }
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    Idx x = bisect_from_left(v->l, k, f, l, mid, ok);
    if(x != inf) return x;
    if(k <= v->idx){
      ok = merge(ok, v->val);
      if(f(ok)) return v->idx;
    }
    return bisect_from_left(v->r, k, f, mid, r, ok);
  }
  template<typename F>
  Idx bisect_from_right(node *v, const Idx k, F &f, Idx l, Idx r, Val &ok){
    if(!v || k < l) return inf;
    Val m = merge(v->sum, ok);
    if(r <= k + 1 && !f(m)){
      ok = m;
      return inf;
    }
    Idx mid = (l >> 1) + (r >> 1) + ((l & r) & 1);
    Idx x = bisect_from_right(v->r, k, f, mid, r, ok);
    if(x != inf) return x;
    if(v->idx <= k){
      ok = merge(v->val, ok);
      if(f(ok)) return v->idx;
    }
    return bisect_from_right(v->l, k, f, l, mid, ok);
  }
  node *root;
public:
  binary_trie_monoid(): root(nullptr){static_assert(bitlen < sizeof(Idx) * 8);}
  void set(Idx k, Val x){
    set(root, k, x, 0, inf);
  }
  Val get(Idx k){
    return get(root, k, 0, inf);
  }
  void update(Idx k, Val x){
    update(root, k, x, 0, inf);
  }
  Val query(Idx l, Idx r){
    return query(root, l, r, 0, inf);
  }
  Val query_all(){
    return root ? root->sum : id();
  }
  // kと接頭辞(bitlen-1, bitlen-2...)がt個以上一致する要素のsum
  Val query_lcp_ge(Idx k, int t){
    if(t > bitlen) return id();
    t = std::max(t, 0);
    int low = bitlen - t;
    Idx p = (~((Idx(1) << low) - 1)) & k;
    return query(p, p + (Idx(1) << low));
  }
  // kと接頭辞(bitlen-1, bitlen-2...)がちょうどt個一致する要素のsum
  Val query_lcp(Idx k, int t){
    if(t > bitlen || t < 0) return id();
    t = std::max(t, 0);
    int low = bitlen - t;
    Idx p = (~((Idx(1) << low) - 1)) & k;
    Idx l = p, r = p + (Idx(1) << low);
    if(t == bitlen) return query(l, r);
    if((k >> (low - 1)) & 1) r -= Idx(1) << (low - 1);
    else l += Idx(1) << (low - 1);
    return query(l, r);
  }
  // ans[t] := kと接頭辞(bitlen-1, bitlen-2...)がちょうどt個一致する要素のsum
  std::vector<Val> query_lcp2(Idx k){
    std::vector<Val> ret(bitlen + 1,id()); // [0, bitlen]
    query_lcp2(root, k, 0, inf, ret);
    return ret;
  }
  // f(sum[l, r])がtrueになる最左のr. ない場合はinf
  // firstがinfでない場合, その時の値
  template<typename F>
  std::pair<Idx, Val> bisect_from_left(Idx l, const F &f){
    Val x = id();
    Idx ret = bisect_from_left(root, l, f, 0, inf, x);
    return {ret, x};
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合はinf
  // firstがinfでない場合, その時の値
  template<typename F>
  std::pair<Idx, Val> bisect_from_right(Idx r, const F &f){
    Val x = id();
    Idx ret = bisect_from_right(root, r, f, 0, inf, x);
    return {ret, x};
  }
};
#endif