#ifndef _BINARY_TRIE_RANGE_BIT_OPERATION_H_
#define _BINARY_TRIE_RANGE_BIT_OPERATION_H_
#include <vector>
#include <cassert>
#include "../../algebraic_structure/monoid.hpp"

// @param 0 < bitlen < (Idxのbit長)
// @param idxがunsigned
template<typename Idx = uint32_t, int bitlen = 30>
struct binary_trie_range_bit_opetarion{
  static constexpr Idx inf = Idx(1) << bitlen;
private:
  struct Lazy{
    Idx __or, __and, __xor;
    Lazy(): __or(0), __and(~(Idx)0), __xor(0){}
    Lazy(Idx a, Idx b, Idx c): __or(a), __and(b), __xor(c){}
    bool is_id(){return __or == 0 && __and == ~(Idx)0 && __xor == 0;}
    void reset(){__or = 0, __and = ~(Idx)0, __xor = 0;}
  };
  static constexpr Idx or_and_xor(Idx i, const Lazy &lz){
    return ((i | lz.__or) & lz.__and) ^ lz.__xor;
  }
  static constexpr void propagate_lazy(Lazy &l, const Lazy &r){
    // 例外の部分だけ, {l.or, l.and, l.xor ^ r.xor}にしたい
    Idx __exception_flag = (~r.__or) & r.__and;
    Idx a = __exception_flag & l.__or;
    Idx b = __exception_flag & l.__and;
    Idx c = __exception_flag & (l.__xor ^ r.__xor);
    l.__or  = ((l.__or & r.__or) ^ l.__xor) | r.__or;
    l.__and = ((l.__and ^ l.__xor) | r.__or) & r.__and;
    l.__xor = r.__xor;
    __exception_flag = ~__exception_flag;
    l.__or  = (l.__or  & __exception_flag) | a;
    l.__and = (l.__and & __exception_flag) | b;
    l.__xor = (l.__xor & __exception_flag) | c;
  }
  struct node{
    int sz;
    Lazy lazy;
    node *l, *r;
    node(int _num = 0): sz(_num), lazy(), l(nullptr), r(nullptr){}
  };
  node *root;
  void eval(node *v, int h){
    if(h == -1) return;
    v->sz = (v->l ? v->l->sz : 0) + (v->r ? v->r->sz : 0);
  }
  // 高さhのノードに作用
  void propagate(node *v, int h, const Lazy &lz){
    if(!v || h == -1) return;
    propagate_lazy(v->lazy, lz);
    Idx L = (or_and_xor(0, lz) >> h) & 1; // 左のノードが左右どちらになるか
    Idx R = (or_and_xor((Idx)1 << h, lz) >> h) & 1; // 右のノードが左右どちらになるか
    if(L == R){
      push_down(v, h);
      if(L){
        v->r = merge(v->l, v->r, h - 1);
        v->l = nullptr;
      }else{
        v->l = merge(v->l, v->r, h - 1);
        v->r = nullptr;
      }
    }else if(L){
      assert(!R);
      std::swap(v->l, v->r);
    }
    eval(v, h);
  }
  void push_down(node *v, int h){
    if(!v || v->lazy.is_id() || h <= 0) return;
    if(v->l) propagate(v->l, h - 1, v->lazy);
    if(v->r) propagate(v->r, h - 1, v->lazy);
    v->lazy.reset();
  }
  node *merge(node *a, node *b, int h){
    if(!a || !b) return !a ? b : a;
    push_down(a, h), push_down(b, h);
    a->sz += b->sz;
    a->l = merge(a->l, b->l, h - 1);
    a->r = merge(a->r, b->r, h - 1);
    eval(a, h);
    return a;
  }
  node *split(node *v, Idx a, Idx b, Idx l, Idx r, int h){
    if(!v || b <= l || r <= a) return nullptr;
    push_down(v, h);
    if(a <= l && r <= b) return v;
    node *u = new node(*v);
    Idx m = (l >> 1) + (r >> 1) + ((l & r) & 1);
    u->l = split(u->l, a, b, l, m, h - 1);
    u->r = split(u->r, a, b, m, r, h - 1);
    if(u->l == v->l) v->l = nullptr;
    if(u->r == v->r) v->r = nullptr;
    eval(u, h);
    eval(v, h);
    return u;
  }
  node *find(node *v, Idx k, int h){
    if(!v || h == -1) return v;
    push_down(v, h);
    if((k >> h) & 1) return find(v->r, k, h - 1);
    else return find(v->l, k, h - 1);
  }
  node *insert(node *v, Idx k, int h, bool unique){
    if(h == -1){
      if(v){
        if(v->sz && unique) return v;
        v->sz++;
        return v;
      }
      return new node(1);
    }
    if(!v) v = new node(0);
    else push_down(v, h);
    if((k >> h) & 1) v->r = insert(v->r, k, h - 1, unique);
    else v->l = insert(v->l, k, h - 1, unique);
    eval(v, h);
    return v;
  }
  void erase(node *v, Idx k, int h, bool all){
    if(!v || !v->sz) return;
    if(h == -1){
      if(all) v->sz = 0;
      else v->sz = std::max(0, v->sz - 1);
      return;
    }
    push_down(v, h);
    if((k >> h) & 1) erase(v->r, k, h - 1, all);
    else erase(v->l, k, h - 1, all);
    eval(v, h);
  }
  // [a, b)の要素数
  int count_range(node *v, Idx a, Idx b, Idx l, Idx r, int h){
    if(!v || b <= l || r <= a) return 0;
    if(a <= l && r <= b) return v->sz;
    push_down(v, h);
    Idx m = (l >> 1) + (r >> 1) + ((l & r) & 1);
    Idx ret = count_range(v->l, a, b, l, m, h - 1) + count_range(v->r, a, b, m, r, h - 1);
    return ret;
  }
  void count_lcp2(node *v, Idx k, Idx l, Idx r, int h, std::vector<int> &ret){
    if(!v) return;
    if(h == -1){
      ret[bitlen] = v->sz;
      return;
    }
    push_down(v, h);
    Idx m = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < m){
      count_lcp2(v->l, k, l, m, h - 1, ret);
      if(v->r) ret[bitlen - 1 - h] += v->r->sz;
    }else{
      if(v->l) ret[bitlen - 1 - h] += v->l->sz;
      count_lcp2(v->r, k, m, r, h - 1, ret);
    }
  }
  void select(node *v, int k, int h, Idx &i){
    if(h == -1) return;
    push_down(v, h);
    int szl = v->l ? v->l->sz : 0;
    if(k < szl) select(v->l, k, h - 1, i);
    else select(v->r, k - szl, h - 1, i += ((Idx)1 << h));
  }
public:
  binary_trie_range_bit_opetarion(): root(nullptr){}
  int size(){
    return root ? root->sz : 0;
  }
  // kの数
  int find(Idx k){
    node *v = find(root, k, bitlen - 1);
    return v ? v->sz : 0;
  }
  //　すでに同じ要素があり, unique = trueの場合何もしない
  void insert(Idx k, bool unique = false){
    root = insert(root, k, bitlen - 1, unique);
  }
  //　all = trueの場合全て, そうでない場合1個だけ消す, 要素kがない場合は何もしない
  void erase(Idx k, bool all = false){
    erase(root, k, bitlen - 1, all);
  }
  // [l, r)を満たす要素の個数
  int count_range(Idx l, Idx r){
    return count_range(root, l, r, 0, inf, bitlen - 1);
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
  // ない場合はinf
  Idx min(){
    return size() ? kth_smallest(0) : inf;
  }
  // ない場合はinf
  Idx max(){
    return size() ? kth_smallest(size() - 1) : inf;
  }
  // k番目(0-indexed)に小さい要素(同じ要素も重複して数える) ない場合はinf
  Idx kth_smallest(int k){
    if(k >= size()) return inf;
    Idx ret = 0;
    select(root, k, bitlen - 1, ret);
    return ret;
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
  // 全体にor x
  void or_all(Idx x){
    propagate(root, bitlen - 1, Lazy(x, ~(Idx)0, 0));
  }
  // 全体にand x
  void and_all(Idx x){
    propagate(root, bitlen - 1, Lazy(0, x, 0));
  }
  // 全体にxor x
  void xor_all(Idx x){
    propagate(root, bitlen - 1, Lazy(0, ~(Idx)0, x));
  }
  // [l, r)を満たす要素にor x
  void or_range(Idx l, Idx r, Idx x){
    if(!root) return;
    node *u = split(root, l, r, 0, inf, bitlen - 1);
    propagate(u, bitlen - 1, Lazy(x, ~(Idx)0, 0));
    root = merge(root, u, bitlen - 1);
  }
  // [l, r)を満たす要素にand x
  void and_range(Idx l, Idx r, Idx x){
    if(!root) return;
    node *u = split(root, l, r, 0, inf, bitlen - 1);
    propagate(u, bitlen - 1, Lazy(0, x, 0));
    root = merge(root, u, bitlen - 1);
  }
  // [l, r)を満たす要素にxor x
  void xor_range(Idx l, Idx r, Idx x){
    if(!root) return;
    node *u = split(root, l, r, 0, inf, bitlen - 1);
    propagate(u, bitlen - 1, Lazy(0, ~(Idx)0, x));
    root = merge(root, u, bitlen - 1);
  }
};

// @param 0 < bitlen < (Idxのbit長)
template<typename Idx, int bitlen, typename monoid>
struct binary_trie_monoid_range_bit_opetarion{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge_val = monoid::merge;
  static constexpr Idx inf = Idx(1) << bitlen;
private:
  struct Lazy{
    Idx __or, __and, __xor;
    Lazy(): __or(0), __and(~(Idx)0), __xor(0){}
    Lazy(Idx a, Idx b, Idx c): __or(a), __and(b), __xor(c){}
    bool is_id(){return __or == 0 && __and == ~(Idx)0 && __xor == 0;}
    void reset(){__or = 0, __and = ~(Idx)0, __xor = 0;}
  };
  static constexpr Idx or_and_xor(Idx i, const Lazy &lz){
    return ((i | lz.__or) & lz.__and) ^ lz.__xor;
  }
  static constexpr void propagate_lazy(Lazy &l, const Lazy &r){
    // 例外の部分だけ, {l.or, l.and, l.xor ^ r.xor}にしたい
    Idx __exception_flag = (~r.__or) & r.__and;
    Idx a = __exception_flag & l.__or;
    Idx b = __exception_flag & l.__and;
    Idx c = __exception_flag & (l.__xor ^ r.__xor);
    l.__or  = ((l.__or & r.__or) ^ l.__xor) | r.__or;
    l.__and = ((l.__and ^ l.__xor) | r.__or) & r.__and;
    l.__xor = r.__xor;
    __exception_flag = ~__exception_flag;
    l.__or  = (l.__or  & __exception_flag) | a;
    l.__and = (l.__and & __exception_flag) | b;
    l.__xor = (l.__xor & __exception_flag) | c;
  }
  struct node{
    Val sum;
    Lazy lazy;
    node *l, *r;
    node(Val _val): sum(_val), lazy(), l(nullptr), r(nullptr){}
  };
  node *root;
  void eval(node *v, int h){
    if(h == -1) return;
    v->sum = merge_val(v->l ? v->l->sum : id(), v->r ? v->r->sum : id());
  }
  // 高さhのノードに作用
  void propagate(node *v, int h, const Lazy &lz){
    if(!v || h == -1) return;
    propagate_lazy(v->lazy, lz);
    Idx L = (or_and_xor(0, lz) >> h) & 1; // 左のノードが左右どちらになるか
    Idx R = (or_and_xor((Idx)1 << h, lz) >> h) & 1; // 右のノードが左右どちらになるか
    if(L == R){
      push_down(v, h);
      if(L){
        v->r = merge(v->l, v->r, h - 1);
        v->l = nullptr;
      }else{
        v->l = merge(v->l, v->r, h - 1);
        v->r = nullptr;
      }
    }else if(L){
      assert(!R);
      std::swap(v->l, v->r);
    }
    eval(v, h);
  }
  void push_down(node *v, int h){
    if(!v || v->lazy.is_id() || h <= 0) return;
    if(v->l) propagate(v->l, h - 1, v->lazy);
    if(v->r) propagate(v->r, h - 1, v->lazy);
    v->lazy.reset();
  }
  node *merge(node *a, node *b, int h){
    if(!a || !b) return !a ? b : a;
    push_down(a, h), push_down(b, h);
    a->sum = merge_val(a->sum, b->sum);
    a->l = merge(a->l, b->l, h - 1);
    a->r = merge(a->r, b->r, h - 1);
    eval(a, h);
    return a;
  }
  node *split(node *v, Idx a, Idx b, Idx l, Idx r, int h){
    if(!v || b <= l || r <= a) return nullptr;
    push_down(v, h);
    if(a <= l && r <= b) return v;
    node *u = new node(*v);
    Idx m = (l >> 1) + (r >> 1) + ((l & r) & 1);
    u->l = split(u->l, a, b, l, m, h - 1);
    u->r = split(u->r, a, b, m, r, h - 1);
    if(u->l == v->l) v->l = nullptr;
    if(u->r == v->r) v->r = nullptr;
    eval(u, h);
    eval(v, h);
    return u;
  }
  node *find(node *v, Idx k, int h){
    if(!v || h == -1) return v;
    push_down(v, h);
    if((k >> h) & 1) return find(v->r, k, h - 1);
    else return find(v->l, k, h - 1);
  }
  node *set(node *v, Idx k, Val x, int h){
    if(h == -1){
      if(v){
        v->sum = x;
        return v;
      }
      return new node(x);
    }
    if(!v) v = new node(0);
    else push_down(v, h);
    if((k >> h) & 1) v->r = set(v->r, k, x, h - 1);
    else v->l = set(v->l, k, x, h - 1);
    eval(v, h);
    return v;
  }
  // [a, b)のsum
  Val query(node *v, Idx a, Idx b, Idx l, Idx r, int h){
    if(!v || b <= l || r <= a) return monoid::template id<Val>();
    if(a <= l && r <= b) return v->sum;
    push_down(v, h);
    Idx m = (l >> 1) + (r >> 1) + ((l & r) & 1);
    return merge_val(query(v->l, a, b, l, m, h - 1), query(v->r, a, b, m, r, h - 1));
  }
  void query_lcp2(node *v, Idx k, Idx l, Idx r, int h, std::vector<Val> &ret){
    if(!v) return;
    if(h == -1){
      ret[bitlen] = v->sum;
      return;
    }
    push_down(v, h);
    Idx m = (l >> 1) + (r >> 1) + ((l & r) & 1);
    if(k < m){
      query_lcp2(v->l, k, l, m, h - 1, ret);
      if(v->r) ret[bitlen - 1 - h] = merge_val(ret[bitlen - 1 - h], v->r->sum);
    }else{
      if(v->l) ret[bitlen - 1 - h] = v->l->sum;
      query_lcp2(v->r, k, m, r, h - 1, ret);
    }
  }
public:
  binary_trie_monoid_range_bit_opetarion(): root(nullptr){}
  void set(Idx k, Val x){
    set(root, k, x, bitlen - 1);
  }
  Val get(Idx k){
    auto v = find(root, k, bitlen - 1);
    return v ? v->sum : id();
  }
  void update(Idx k, Val x){
    update(root, k, x, bitlen - 1);
  }
  Val query(Idx l, Idx r){
    return query(root, l, r, 0, inf, bitlen - 1);
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
    std::vector<Val> ret(bitlen + 1, id()); // [0, bitlen]
    query_lcp2(root, k, 0, inf, bitlen - 1, ret);
    return ret;
  }
  // 全体にor x
  void or_all(Idx x){
    propagate(root, bitlen - 1, Lazy(x, ~(Idx)0, 0));
  }
  // 全体にand x
  void and_all(Idx x){
    propagate(root, bitlen - 1, Lazy(0, x, 0));
  }
  // 全体にxor x
  void xor_all(Idx x){
    propagate(root, bitlen - 1, Lazy(0, ~(Idx)0, x));
  }
  // [l, r)を満たす要素にor x
  void or_range(Idx l, Idx r, Idx x){
    if(!root) return;
    node *u = split(root, l, r, 0, inf, bitlen - 1);
    propagate(u, bitlen - 1, Lazy(x, ~(Idx)0, 0));
    root = merge(root, u, bitlen - 1);
  }
  // [l, r)を満たす要素にand x
  void and_range(Idx l, Idx r, Idx x){
    if(!root) return;
    node *u = split(root, l, r, 0, inf, bitlen - 1);
    propagate(u, bitlen - 1, Lazy(0, x, 0));
    root = merge(root, u, bitlen - 1);
  }
  // [l, r)を満たす要素にxor x
  void xor_range(Idx l, Idx r, Idx x){
    if(!root) return;
    node *u = split(root, l, r, 0, inf, bitlen - 1);
    propagate(u, bitlen - 1, Lazy(0, ~(Idx)0, x));
    root = merge(root, u, bitlen - 1);
  }
};

#endif