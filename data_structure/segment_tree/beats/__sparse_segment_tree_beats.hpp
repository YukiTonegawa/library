#ifndef _SPARSE_SEGMENT_TREE_BEATS_H_
#define _SPARSE_SEGMENT_TREE_BEATS_H_
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>

// Val : 値1個の型, ValSum : 和の型
template<typename Val, typename ValSum>
struct sparse_segment_tree_beats{
  using I = int;
  constexpr static I inf = std::numeric_limits<I>::max();
  constexpr static I minf = std::numeric_limits<I>::min();
  constexpr static Val INF = std::numeric_limits<Val>::max();
  constexpr static Val MINF = std::numeric_limits<Val>::min();
private:
  struct Data{
    Val min, second_min, max, second_max, val;
    ValSum sum;
    I min_cnt, max_cnt;
    Data(I num, Val _val): min(_val), second_min(INF), max(_val), second_max(MINF), val(_val), sum((ValSum)_val * num), min_cnt(num), max_cnt(num){}
  };
  struct Lazy{
    Val add, lower, upper;
    Lazy(): add(0), lower(MINF), upper(INF){}
    Lazy(Val _a, Val _b, Val _c): add(_a), lower(_b), upper(_c){}
    bool is_id(){return add == 0 && lower == MINF && upper == INF;}
    void reset(){add = 0, lower = MINF, upper = INF;}
  };
  struct node{
    int h;
    I lx, rx, lmin, rmax;
    Data data;
    Lazy lazy;
    node *l, *r;
    node(I lx, I rx, Val _val): h(1), lx(lx), rx(rx), lmin(lx), rmax(rx), data(rx - lx, _val), lazy(), l(nullptr), r(nullptr){}
    int balanace_factor(){return (l ? l->h : 0) - (r ? r->h : 0);}
  };
  node *root;
  node *make_node(I l, I r, Val x){return new node(l, r, x);}
  void push_down(node *v){
    if(!v->lazy.is_id()){
      if(v->l) propagate(v->l, v->lazy);
      if(v->r) propagate(v->r, v->lazy);
      v->lazy.reset();
    }
  }
  void propagate_lazy(Lazy &a, const Lazy &b){
    if(b.add){
      a.add += b.add;
      if(a.lower != MINF) a.lower += b.add;
      if(a.upper != INF) a.upper += b.add;
    }
    if(a.upper <= b.lower) a.lower = a.upper = b.lower;
    else if(b.upper <= a.lower) a.lower = a.upper = b.upper;
    else{
      a.lower = std::max(a.lower, b.lower);
      a.upper = std::min(a.upper, b.upper);
    }
  }
  void propagate(node *v, const Lazy &x){
    propagate_lazy(v->lazy, x);
    if(x.add){
      v->data.val += x.add;
      v->data.min += x.add;
      v->data.max += x.add;
      if(v->data.second_min != INF) v->data.second_min += x.add;
      if(v->data.second_max != MINF) v->data.second_max += x.add;
      v->data.sum += (ValSum)x.add * (v->rmax - v->lmin);
    }
    if(v->data.val < x.lower) v->data.val = x.lower;
    if(v->data.val > x.upper) v->data.val = x.upper;
    if(v->data.min < x.lower){
      v->data.sum += (ValSum)(x.lower - v->data.min) * v->data.min_cnt;
      if(v->data.second_max == v->data.min) v->data.second_max = x.lower;
      else if(v->data.max == v->data.min) v->data.max = x.lower, v->data.second_max = MINF;
      v->data.min = x.lower;
    }
    if(x.upper < v->data.max){
      v->data.sum -= (ValSum)(v->data.max - x.upper) * v->data.max_cnt;
      if(v->data.second_min == v->data.max) v->data.second_min = x.upper;
      else if(v->data.min == v->data.max) v->data.min = x.upper, v->data.second_min = INF;
      v->data.max = x.upper;
    }
  }
  void update_val(Data &l, const Data &r){
    l.sum += r.sum;
    if(l.max == r.max){
      l.max_cnt += r.max_cnt;
      l.second_max = std::max(l.second_max, r.second_max);
    }else if(l.max > r.max){
      l.second_max = std::max(l.second_max, r.max);
    }else{
      l.max_cnt = r.max_cnt;
      l.second_max = std::max(l.max, r.second_max);
      l.max = r.max;
    }
    if(l.min == r.min){
      l.min_cnt += r.min_cnt;
      l.second_min = std::min(l.second_min, r.second_min);
    }else if(l.min < r.min){
      l.second_min = std::min(l.second_min, r.min);
    }else{
      l.min_cnt = r.min_cnt;
      l.second_min = std::min(l.min, r.second_min);
      l.min = r.min;
    }
  }
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->lmin = v->l ? v->l->lmin : v->lx;
    v->rmax = v->r ? v->r->rmax : v->rx;
    // 呼ばれた時点でv->data.valだけが正しいと仮定
    I len = v->rx - v->lx;
    v->data.sum = (ValSum)v->data.val * len;
    v->data.min = v->data.max = v->data.val;
    v->data.min_cnt = v->data.max_cnt = len;
    v->data.second_max = MINF;
    v->data.second_min = INF;
    if(v->l) update_val(v->data, v->l->data);
    if(v->r) update_val(v->data, v->r->data);
  }
  // vの区間がこれよりも大きい場合のみ使える
  std::pair<node*, node*> split_node(node *v, I l, I r){
    node *x = nullptr, *y = nullptr;
    if(v->lx < l) x = new node(v->lx, l, v->data.val);
    if(r < v->rx) y = new node(r, v->rx, v->data.val);
    v->lx = l, v->rx = r;
    return {x, y};
  }
  node *rotate_right(node *v){
    push_down(v);
    node *l = v->l;
    push_down(l);
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  node *rotate_left(node *v){
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
  node *update_add_inner(node *v, I l, I r, Val x){
    if(l == v->lmin && r == v->rmax){
      propagate(v, Lazy(x, MINF, INF));
      return v;
    }
    push_down(v);
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx){
        v->l = update_add_inner(v->l, l, r, x);
        update(v);
        return balance(v);
      }
      v->l = update_add_inner(v->l, l, v->lx, x);
      l = v->lx;
    }
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l){
        v->r = update_add_inner(v->r, l, r, x);
        update(v);
        return balance(v);
      }
      v->r = update_add_inner(v->r, v->rx, r, x);
      r = v->rx;
    }
    assert(r > l);
    if(l != v->lx || r != v->rx){
      auto [_l, _r] = split_node(v, l, r);
      if(_l) v->l = insert_rightmost(v->l, _l);
      if(_r) v->r = insert_leftmost(v->r, _r);
    }
    v->data.val += x;
    update(v);
    return balance(v);
  }
  node *update_chmin_inner(node *v, I l, I r, Val x){
    if(v->data.max < x) return v; // break
    // 部分木の区間と一致する場合
    if(l == v->lmin && r == v->rmax && v->data.second_max < x){ // tag
      propagate(v, Lazy(0, MINF, x));
      return v;
    }
    push_down(v);
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx){
        v->l = update_chmin_inner(v->l, l, r, x);
        update(v);
        return balance(v);
      }
      v->l = update_chmin_inner(v->l, l, v->lx, x);
      l = v->lx;
    }
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l){
        v->r = update_chmin_inner(v->r, l, r, x);
        update(v);
        return balance(v);
      }
      v->r = update_chmin_inner(v->r, v->rx, r, x);
      r = v->rx;
    }
    assert(r > l);
    if(l != v->lx || r != v->rx){
      auto [_l, _r] = split_node(v, l, r);
      if(_l) v->l = insert_rightmost(v->l, _l);
      if(_r) v->r = insert_leftmost(v->r, _r);
    }
    v->data.val = std::min(v->data.val, x);
    update(v);
    return balance(v);
  }
  node *update_chmax_inner(node *v, I l, I r, Val x){
    if(v->data.min > x) return v; // break
    // 部分木の区間と一致する場合
    if(l == v->lmin && r == v->rmax && v->data.second_min > x){ // tag
      propagate(v, Lazy(0, x, INF));
      return v;
    }
    push_down(v);
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx){
        v->l = update_chmax_inner(v->l, l, r, x);
        update(v);
        return balance(v);
      }
      v->l = update_chmax_inner(v->l, l, v->lx, x);
      l = v->lx;
    }
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l){
        v->r = update_chmax_inner(v->r, l, r, x);
        update(v);
        return balance(v);
      }
      v->r = update_chmax_inner(v->r, v->rx, r, x);
      r = v->rx;
    }
    assert(r > l);
    if(l != v->lx || r != v->rx){
      auto [_l, _r] = split_node(v, l, r);
      if(_l) v->l = insert_rightmost(v->l, _l);
      if(_r) v->r = insert_leftmost(v->r, _r);
    }
    v->data.val = std::max(v->data.val, x);
    update(v);
    return balance(v);
  }
  ValSum query_sum_inner(node *v, I l, I r){
    if(l == v->lmin && r == v->rmax) return v->data.sum;
    push_down(v);
    ValSum left_sum = 0;
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx) return query_sum_inner(v->l, l, r);
      left_sum = query_sum_inner(v->l, l, v->lx);
      l = v->lx;
    }
    ValSum right_sum = 0;
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l) return query_sum_inner(v->r, l, r);
      right_sum = query_sum_inner(v->r, v->rx, r);
      r = v->rx;
    }
    assert(r > l);
    return (ValSum)v->data.val * (r - l) + left_sum + right_sum;
  }
  Val query_min_inner(node *v, I l, I r){
    if(l == v->lmin && r == v->rmax) return v->data.min;
    push_down(v);
    Val left_sum = INF;
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx) return query_min_inner(v->l, l, r);
      left_sum = query_min_inner(v->l, l, v->lx);
      l = v->lx;
    }
    Val right_sum = INF;
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l) return query_min_inner(v->r, l, r);
      right_sum = query_min_inner(v->r, v->rx, r);
      r = v->rx;
    }
    assert(r > l);
    return std::min({v->data.val, left_sum, right_sum});
  }
  Val query_max_inner(node *v, I l, I r){
    if(l == v->lmin && r == v->rmax) return v->data.max;
    push_down(v);
    Val left_sum = MINF;
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx) return query_max_inner(v->l, l, r);
      left_sum = query_max_inner(v->l, l, v->lx);
      l = v->lx;
    }
    Val right_sum = MINF;
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l) return query_max_inner(v->r, l, r);
      right_sum = query_max_inner(v->r, v->rx, r);
      r = v->rx;
    }
    assert(r > l);
    return std::max({v->data.val, left_sum, right_sum});
  }
  node *build(const std::vector<node*> &nodes, int l, int r){
    int m = (l + r) >> 1;
    node *v = nodes[m];
    if(m > l) v->l = build(nodes, l, m);
    if(r > m + 1) v->r = build(nodes, m + 1, r);
    update(v);
    return v;
  }
  void init_points_inner(std::vector<std::pair<I, Val>> v, bool is_sorted = false){
    // 座標をソート
    if(!is_sorted) std::sort(v.begin(), v.end());
    if(v.empty()){
      root = new node(minf, inf, 0);
      return;
    }
    std::vector<node*> p;
    for(int i = 0; i < v.size(); i++){
      I prev_right = (!i ? minf : p.back()->rx);
      assert(prev_right <= v[i].first);
      if(prev_right < v[i].first) p.push_back(new node(prev_right, v[i].first, 0));
      p.push_back(new node(v[i].first, v[i].first + 1, v[i].second));
    }
    p.push_back(new node(p.back()->rx, inf, 0));
    root = build(p, 0, p.size());
  }
public:
  sparse_segment_tree_beats(): root(new node(minf, inf, 0)){}

  // 複数の重複しない点から構築, {点の位置, 値}
  // 0 <= v.size()
  void init_points(const std::vector<std::pair<I, Val>> &v, bool is_sorted = false){
    init_points_inner(v, is_sorted);
  }
  void update_add(I l, I r, Val x){
    if(l >= r) return;
    root = update_add_inner(root, l, r, x);
  }
  void update_chmin(I l, I r, Val x){
    if(l >= r) return;
    root = update_chmin_inner(root, l, r, x);
  }
  void update_chmax(I l, I r, Val x){
    if(l >= r) return;
    root = update_chmax_inner(root, l, r, x);
  }
  ValSum query_sum(I l, I r){
    if(l >= r) return 0;
    return query_sum_inner(root, l, r);
  }
  Val query_min(I l, I r){
    if(l >= r) return 0;
    return query_min_inner(root, l, r);
  }
  Val query_max(I l, I r){
    if(l >= r) return 0;
    return query_max_inner(root, l, r);
  }
  Data query_all(){
    return root->data;
  }
};
#endif