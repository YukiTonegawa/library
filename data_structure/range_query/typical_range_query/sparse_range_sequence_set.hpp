#ifndef _SPARSE_RANGE_SEQUENCE_SET_H_
#define _SPARSE_RANGE_SEQUENCE_SET_H_
#include <vector>
#include <tuple>
#include <cassert>
#include <iostream>
#include <limits>
template<typename Sequence>
struct sparse_range_sequence_set{
private:
  using Val = typename Sequence::Val;
  using Lazy = typename Sequence::Lazy;
  static constexpr int inf = (1 << 30) - 1;
  static constexpr int minf = -inf;
  struct node{
    int h;
    node *l, *r;
    int lx, rx, minlx, maxrx;
    Val val, sum;
    Lazy lz;
    node(Lazy f, int _lx, int _rx): h(1), l(nullptr), r(nullptr), lx(_lx), rx(_rx), minlx(lx), maxrx(rx), val(Sequence::get_sum(f, lx, rx)), sum(val), lz(f){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root, *tmp_node;
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sum = v->val;
    v->minlx = v->lx;
    v->maxrx = v->rx;
    if(v->l) v->sum = Sequence::merge(v->l->sum, v->sum), v->minlx = v->l->minlx;
    if(v->r) v->sum = Sequence::merge(v->sum, v->r->sum), v->maxrx = v->r->maxrx;
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
  node *cut_leftmost(node *v){
    if(v->l){
      v->l = cut_leftmost(v->l);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  node *insert_inner(node *v, node *u){
    if(!v) return u;
    if(u->lx < v->lx){
      v->l = insert_inner(v->l, u);
    }else{
      v->r = insert_inner(v->r, u);
    }
    update(v);
    return balance(v);
  }
  node *erase_inner(node *v, int l, int r, bool &f){
    if(!v) return nullptr;
    if(v->maxrx <= l || r <= v->minlx) return v;
    int L = std::max(v->lx, l), R = std::min(v->rx, r);
    if(L == v->lx && R == v->rx){
      f = true;
      // vを消す
      if(v->r){
        v->r = cut_leftmost(v->r);
        tmp_node->l = v->l;
        tmp_node->r = v->r;
        update(tmp_node);
        return balance(tmp_node);
      }
      return v->l;
    }
    if(L < R){
      if(v->lx < l && r < v->rx){
        f = true;
        v->val = Sequence::get_sum(v->lz, v->lx, l);
        node *u = new node(v->lz, r, v->rx);
        v->rx = l;
        v->r = insert_inner(v->r, u);
        update(v);
        return balance(v);
      }else if(v->lx < l){
        v->val = Sequence::get_sum(v->lz, v->lx, l);
        v->rx = l;
        update(v);
      }else{
        v->val = Sequence::get_sum(v->lz, r, v->rx);
        v->lx = r;
        update(v);
      }
    }
    v->l = erase_inner(v->l, l, r, f);
    if(f){
      update(v);
      return balance(v);
    }
    v->r = erase_inner(v->r, l, r, f);
    update(v);
    return balance(v);
  }
  void update_inner(int l, int r, Lazy f){
    bool ok = true;
    while(ok){
      ok = false;
      root = erase_inner(root, l, r, ok);
    }
    root = insert_inner(root, new node(f, l, r));
  }
  Val query_inner(int l, int r, node *v){
    if(l <= v->minlx && v->maxrx <= r) return v->sum;
    Val Lsum = Sequence::id();
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx) return query_inner(l, r, v->l);
      Lsum = query_inner(l, v->lx, v->l);
      l = v->lx;
    }
    Val Rsum = Sequence::id();
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l) return query_inner(l, r, v->r);
      Rsum = query_inner(v->rx, r, v->r);
      r = v->rx;
    }
    Val x = Sequence::id();
    assert(r > l);
    // 更新対象の区間が現在の区間と一致するならそのままの値, 一致しないなら[l, r)の部分を再計算
    if(l == v->lx && r == v->rx) x = v->val;
    else x = Sequence::get_sum(v->lz, l, r);
    return Sequence::merge(Lsum, Sequence::merge(x, Rsum));
  }
  node *build(int l, int r, const std::vector<std::tuple<int, int, Lazy>> &v){
    int m = (l + r) / 2;
    auto [lx, rx, f] = v[m];
    node *u = new node(f, lx, rx);
    if(l < m) u->l = build(l, m, v);
    if(m + 1 < r) u->r = build(m + 1, r, v);
    update(u);
    return u;
  }
  void print(node *v){
    if(!v) return;
    if(v->l) print(v->l);
    std::cout << v->lx << " " << v->rx << " " << v->lz << '\n';
    if(v->r) print(v->r);
  }
public:
  sparse_range_sequence_set(): root(nullptr){}
  sparse_range_sequence_set(const std::vector<std::tuple<int, int, Lazy>> &v): root(build(0, v.size(), v)){}
  // [l, r)を更新
  void update(int l, int r, Lazy f){
    update_inner(l, r, f);
  }
  // [l, r)の積
  Val query(int l, int r){
    return query_inner(l, r, root);
  }
  void print(){
    print(root);
    std::cout << "------" << '\n';
  }
};
#endif