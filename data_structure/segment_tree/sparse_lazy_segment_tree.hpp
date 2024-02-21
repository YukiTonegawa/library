#ifndef _SPARSE_LAZY_SEGMENT_TREE_H_
#define _SPARSE_LAZY_SEGMENT_TREE_H_
#include <vector>
#include <cassert>
#include "../../algebraic_structure/monoid.hpp"
template<typename monoid>
struct sparse_lazy_segment_tree{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
  using I = int;
  constexpr static I inf = std::numeric_limits<I>::max() / 2;
  constexpr static I minf = std::numeric_limits<I>::min() / 2;
private:
  struct node{
    int h;
    I lx, rx, lmin, rmax;
    Val val, sum;
    bool is_reset;
    Lazy lazy, lazy_sum;
    node *l, *r;
    node(I lx, I rx, Val _val = id()): h(1), lx(lx), rx(rx),
    lmin(lx), rmax(rx), val(_val), sum(val), is_reset(false), lazy(id_lazy()), lazy_sum(lazy), l(nullptr), r(nullptr){}
    int balanace_factor(){return (l ? l->h : 0) - (r ? r->h : 0);}
  };
  node *root;
  node *make_node(I l, I r, Val x){
    return new node(l, r, x);
  }
  // 元からある区間を[a, b)とする, a <= l < r <= bを満たす区間[l, r)が新しく追加された時に
  // [a, l), [l, r), [r, b)に分割する関数が必要[1]
  // (区間が交差する場合は場合分けで処理できるため, このパターンのみ対処が必要)
  // 各ノードごとに過去に遅延伝播されたものの総和を持っておく.
  // simulate_valueは区間長を変更した時にこの情報を使って値を再計算する関数
  // split_nodeは[1]を行う関数
  Val simulate_value(node *v, I lx, I rx){
    return apply(id(), v->lazy_sum, lx, rx);
  }
  // vの区間がこれよりも大きい場合のみ使える
  std::pair<node*, node*> split_node(node *v, I l, I r){
    node *x = nullptr, *y = nullptr;
    if(v->lx < l){
      x = new node(v->lx, l, simulate_value(v, v->lx, l));
      x->lazy_sum = v->lazy_sum;
    }
    if(r < v->rx){
      y = new node(r, v->rx, simulate_value(v, r, v->rx));
      y->lazy_sum = v->lazy_sum;
    }
    v->lx = l, v->rx = r;
    v->val = simulate_value(v, l, r);
    return {x, y};
  }
  static void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->lmin = v->lx;
    v->rmax = v->rx;
    v->sum = v->val;
    if(v->l){
      v->lmin = v->l->lmin;
      v->sum = merge(v->l->sum, v->sum);
    }
    if(v->r){
      v->rmax = v->r->rmax;
      v->sum = merge(v->sum, v->r->sum);
    }
  }
  static void propagate(node *v, Lazy x){
    v->lazy = propagate_lazy(v->lazy, x);
    v->lazy_sum = propagate_lazy(v->lazy_sum, x);
    v->val = apply(v->val, x, v->lx, v->rx);
    v->sum = apply(v->sum, x, v->lmin, v->rmax);
  }
  static void reset(node *v){
    v->is_reset = true;
    v->lazy = v->lazy_sum = id_lazy();
    v->val = v->sum = id();
  }
  static void push_down(node *v){
    if(v->is_reset){
      if(v->l) reset(v->l);
      if(v->r) reset(v->r);
      v->is_reset = false;
    }
    if(v->lazy != id_lazy()){
      if(v->l) propagate(v->l, v->lazy);
      if(v->r) propagate(v->r, v->lazy);
      v->lazy = id_lazy();
    }
  }
  static node *rotate_right(node *v){
    push_down(v);
    node *l = v->l;
    push_down(l);
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  static node *rotate_left(node *v){
    push_down(v);
    node *r = v->r;
    push_down(r);
    v->r = r->l;
    r->l = v;
    update(v);
    update(r);
    return r;
  }
  node *balance(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    push_down(v);
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
  node *insert_leftmost(node *v, node *u){
    if(!v) return u;
    push_down(v);
    v->l = insert_leftmost(v->l, u);
    update(v);
    return balance(v);
  }
  node *insert_rightmost(node *v, node *u){
    if(!v) return u;
    push_down(v);
    v->r = insert_rightmost(v->r, u);
    update(v);
    return balance(v);
  }
  node *reset_inner(node *v, I l, I r){
    // 部分木の区間と一致する場合
    if(l == v->lmin && r == v->rmax){
      reset(v);
      return v;
    }
    push_down(v);
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx){
        v->l = reset_inner(v->l, l, r);
        update(v);
        return balance(v);
      }
      v->l = reset_inner(v->l, l, v->lx);
      l = v->lx;
    }
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l){
        v->r = reset_inner(v->r, l, r);
        update(v);
        return balance(v);
      }
      v->r = reset_inner(v->r, v->rx, r);
      r = v->rx;
    }
    assert(r > l);
    // 更新対象の区間が現在の区間と一致するなら値のみ更新, 一致しないなら区間を分割
    if(l != v->lx || r != v->rx){
      auto [_l, _r] = split_node(v, l, r);
      if(_l) v->l = insert_rightmost(v->l, _l);
      if(_r) v->r = insert_leftmost(v->r, _r);
    }
    v->val = id();
    v->lazy_sum = id_lazy();
    update(v);
    return balance(v);
  }
  node *update_inner(node *v, I l, I r, Lazy x){
    // 部分木の区間と一致する場合
    if(l == v->lmin && r == v->rmax){
      propagate(v, x);
      return v;
    }
    push_down(v);
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx){
        v->l = update_inner(v->l, l, r, x);
        update(v);
        return balance(v);
      }
      v->l = update_inner(v->l, l, v->lx, x);
      l = v->lx;
    }
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l){
        v->r = update_inner(v->r, l, r, x);
        update(v);
        return balance(v);
      }
      v->r = update_inner(v->r, v->rx, r, x);
      r = v->rx;
    }
    assert(r > l);
    // 更新対象の区間が現在の区間と一致するなら値のみ更新, 一致しないなら区間を分割
    if(l != v->lx || r != v->rx){
      auto [_l, _r] = split_node(v, l, r);
      if(_l) v->l = insert_rightmost(v->l, _l);
      if(_r) v->r = insert_leftmost(v->r, _r);
    }
    v->val = apply(v->val, x, l, r);
    v->lazy_sum = propagate_lazy(v->lazy_sum, x);
    update(v);
    return balance(v);
  }
  Val query_inner(node *v, I l, I r){
    if(l == v->lmin && r == v->rmax) return v->sum;
    push_down(v);
    Val left_sum = id();
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx) return query_inner(v->l, l, r);
      left_sum = query_inner(v->l, l, v->lx);
      l = v->lx;
    }
    Val right_sum = id();
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l) return query_inner(v->r, l, r);
      right_sum = query_inner(v->r, v->rx, r);
      r = v->rx;
    }
    Val res = id();
    assert(r > l);
    // 更新対象の区間が現在の区間と一致するならそのままの値, 一致しないなら[l, r)の部分を再計算
    if(l == v->lx && r == v->rx){
      res = v->val;
    }else{
      res = simulate_value(v, l, r);
    }
    return merge(left_sum, merge(res, right_sum));
  }
  void enumerate(node *v, std::vector<std::tuple<I, I, Val>> &x){
    if(v->l) enumerate(v->l, x);
    x.push_back({v->lx, v->rx, v->val});
    if(v->r) enumerate(v->r, x);
  }
  node *set_inner(node *v, I k, Lazy x){
    push_down(v);
    if(k < v->lx){
      v->l = set_inner(v->l, k, x);
    }else if(v->rx <= k){
      v->r = set_inner(v->r, k, x);
    }else{
      auto [_l, _r] = split_node(v, k, k + 1);
      if(_l) v->l = insert_rightmost(v->l, _l);
      if(_r) v->r = insert_leftmost(v->r, _r);
      v->val = apply(id(), x, k, k + 1);
      v->lazy_sum = x; // 遅延伝播の履歴をリセット
    }
    update(v);
    return balance(v);
  }
  Val get_inner(node *v, I k){
    push_down(v);
    if(v->lx <= k && k < v->rx) return simulate_value(v, k, k + 1);
    if(k < v->lx) return get_inner(v->l, k);
    return get_inner(v->r, k);
  }
  template<typename F>
  node* bisect_from_left(node *v, F &f, I l, Val &ok){
    if(!v || v->rmax <= l) return nullptr;
    Val m = merge(ok, v->sum);
    if(l <= v->lmin && !f(m)){
      ok = m;
      return nullptr;
    }
    push_down(v);
    node *u = bisect_from_left(v->l, f, l, ok);
    if(u) return u;
    if(l < v->rx){
      I lx = std::max(v->lx, l);
      Val tmp = merge(ok, simulate_value(v, lx, v->rx));
      if(f(tmp)) return v;
      else ok = tmp;
    }
    return bisect_from_left(v->r, f, l, ok);
  }
  template<typename F>
  node* bisect_from_right(node *v, F &f, I r, Val &ok){
    if(!v || r < v->lmin) return nullptr;
    Val m = merge(ok, v->sum);
    if(r >= v->rmax - 1 && !f(m)){
      ok = m;
      return nullptr;
    }
    push_down(v);
    node *u = bisect_from_right(v->r, f, r, ok);
    if(u) return u;
    if(v->lx <= r){
      I rx = std::min(v->rx, r + 1);
      Val tmp = merge(ok, simulate_value(v, v->lx, rx));
      if(f(tmp)) return v;
      else ok = tmp;
    }
    return bisect_from_right(v->l, f, r, ok);
  }
  template<typename F>
  node* bisect_from_left2(node *v, F &f, I l, Val &ok){
    if(!v || v->rmax <= l) return nullptr;
    Val m = merge(ok, v->sum);
    if(l <= v->lmin && !f(m, l, v->rmax)){
      ok = m;
      return nullptr;
    }
    push_down(v);
    node *u = bisect_from_left2(v->l, f, l, ok);
    if(u) return u;
    if(l < v->rx){
      I lx = std::max(v->lx, l);
      Val tmp = merge(ok, simulate_value(v, lx, v->rx));
      if(f(tmp, l, v->rx)) return v;
      else ok = tmp;
    }
    return bisect_from_left2(v->r, f, l, ok);
  }
  template<typename F>
  node* bisect_from_right2(node *v, F &f, I r, Val &ok){
    if(!v || r < v->lmin) return nullptr;
    Val m = merge(ok, v->sum);
    if(r >= v->rmax - 1 && !f(m, v->lmin, r + 1)){
      ok = m;
      return nullptr;
    }
    push_down(v);
    node *u = bisect_from_right2(v->r, f, r, ok);
    if(u) return u;
    if(v->lx <= r){
      I rx = std::min(v->rx, r + 1);
      Val tmp = merge(ok, simulate_value(v, v->lx, rx));
      if(f(tmp, v->lx, r + 1)) return v;
      else ok = tmp;
    }
    return bisect_from_right2(v->l, f, r, ok);
  }
public:
  sparse_lazy_segment_tree(): root(new node(minf, inf, id())){}

  // 注: 値でなく作用素
  void set(I k, Lazy x){
    assert(minf < k && k < inf);
    root = set_inner(root, k, x);
  }
  Val get(I k){
    assert(minf < k && k < inf);
    return get_inner(root, k);
  }
  void reset(I l, I r){
    assert(minf < l && r < inf);
    if(l >= r) return;
    root = reset_inner(root, l, r);
  }
  void update(I l, I r, Lazy x){
    assert(minf < l && r < inf);
    if(l >= r) return;
    root = update_inner(root, l, r, x);
  }
  Val query(I l, I r){
    assert(minf < l && r < inf);
    if(l >= r) return id();
    return query_inner(root, l, r);
  }
  Val query_all(){
    return root->sum;
  }
  std::vector<std::tuple<I, I, Val>> enumerate(){
    std::vector<std::tuple<I, I, Val>> res;
    enumerate(root, res);
    return res;
  }
  template<typename F>
  std::pair<Val, node*> bisect_from_left(I l, const F f){
    Val x = id();
    node *v = bisect_from_left(root, f, l, x);
    return {x, v};
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合はminf
  template<typename F>
  std::pair<Val, node*> bisect_from_right(I r, const F f){
    Val x = id();
    node *v = bisect_from_right(root, f, r, x);
    return {x, v};
  }
  template<typename F>
  std::pair<Val, node*> bisect_from_left2(I l, const F f){
    Val x = id();
    node *v = bisect_from_left2(root, f, l, x);
    return {x, v};
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合はminf
  template<typename F>
  std::pair<Val, node*> bisect_from_right2(I r, const F f){
    Val x = id();
    node *v = bisect_from_right2(root, f, r, x);
    return {x, v};
  }
};
#endif