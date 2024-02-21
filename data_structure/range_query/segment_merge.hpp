#ifndef _SEGMENT_MERGE_H_
#define _SEGMENT_MERGE_H_
#include <vector>
#include <numeric>
#include <map>

template<typename Idx>
struct segment_merge{
  const Idx INF = std::numeric_limits<Idx>::max();
private:
  std::map<Idx, Idx> mp;
public:
  typename std::map<Idx, Idx>::iterator begin(){return mp.begin();}
  typename std::map<Idx, Idx>::iterator end(){return mp.end();}
  segment_merge(){}
  // [l, r)を追加して重なる区間があればマージ
  // merge_adjacent:
  // trueなら[r, x) があった時に[l, x)になる
  // falseなら[r, x)があっても[l, r), [r, x)はマージされない
  void insert(Idx l, Idx r, bool merge_adjacent=true){
    auto itr = mp.upper_bound(l);
    if(itr!=mp.begin()) itr--;
    if(itr!=mp.end()&&(merge_adjacent?itr->second<l:itr->second<=l)) itr++;
    while(itr!=mp.end()&&(merge_adjacent?itr->first<=r:itr->first<r)){
      r = std::max(r, itr->second);
      l = std::min(l, itr->first);
      itr = mp.erase(itr);
    }
    mp.emplace(l, r);
  }
  // l <= a < b <= r　となる[a, b)を全列挙
  std::vector<std::pair<Idx, Idx>> enumerate_include(Idx l, Idx r){
    auto itr = mp.lower_bound(l);
    std::vector<std::pair<Idx, Idx>> ret;
    while(itr!=mp.end()&&itr->second<=r) ret.push_back(*itr++);
    return ret;
  }
  // min(b, r) > max(a, l) となる[a, b)を全列挙
  std::vector<std::pair<Idx, Idx>> enumerate_intersect(Idx l, Idx r){
    auto itr = mp.upper_bound(l);
    std::vector<std::pair<Idx, Idx>> ret;
    while(itr!=mp.begin()&&(itr==mp.end()||itr->second>l)) itr--;
    if(itr!=mp.end()&&itr->second<=l) itr++;
    while(itr!=mp.end()&&itr->first<r) ret.push_back(*itr++);
    return ret;
  }
  // l <= a < b <= r　となる[a, b)を全削除
  void erase_include(Idx l, Idx r){
    auto itr = mp.lower_bound(l);
    while(itr!=mp.end()&&itr->second<=r) itr = mp.erase(itr);
  }
  // min(b, r) > max(a, l) となる[a, b)を全削除
  void erase_intersect(Idx l, Idx r){
    auto itr = mp.upper_bound(l);
    while(itr!=mp.begin()&&(itr==mp.end()||itr->second>l)) itr--;
    if(itr!=mp.end()&&itr->second<=l) itr++;
    while(itr!=mp.end()&&itr->first<r) itr = mp.erase(itr);
  }
  Idx isolated(Idx l, Idx r){
    Idx ret = Idx(0);
    Idx itr = l;
    for(auto &s:enumerate_intersect(l, r)){
      if(itr < s.first) ret += s.first - itr;
      itr = s.second;
    }
    return ret + r - std::min(itr, r);
  }
  // 任意の区間に含まれないa以上の　最小の非負整数
  Idx mex(Idx a = 0){
    auto itr = mp.upper_bound(a);
    if(itr == mp.begin()) return a;
    itr--;
    if(itr->second <= a) return a;
    return itr->second;
  }
  // 任意の区間に含まれないa以下の　最大の非負整数
  Idx mex_rev(Idx a = 0){
    auto itr = mp.upper_bound(a);
    if(itr == mp.begin()) return a;
    itr--;
    if(itr->second <= a) return a;
    return itr->first - 1;
  }
  // pを含む区間, ない場合は{INF, INF}
  std::pair<Idx, Idx> find(Idx p){
    auto itr = mp.upper_bound(p);
    if(itr!=mp.begin() && (--itr)->second > p) return *itr;
    else return {INF, INF};
  }
  Idx min(){
    if(mp.empty()) return INF;
    return mp.begin()->first;
  }
  Idx max(){
    if(mp.empty()) return INF;
    return (--mp.end())->second - 1;
  }
  // p, q が同じ区間に含まれるか
  bool same(Idx p, Idx q){
    std::pair<Idx, Idx> l = find(p), r = find(q);
    return l.first!=INF && l.first==r.first;
  }
  // rの要素を全部マージ(rの要素は全部消える)
  void meld(segment_merge<Idx> &r){
    if(mp.size() < r.mp.size()) mp.swap(r.mp);
    for(auto [l, r] : r.mp) insert(l, r);
    r.mp.clear();
  }
};

/*
template<typename Idx>
struct segment_merge_super{
  constexpr static Idx inf = std::numeric_limits<Idx>::max() / 2;
  constexpr static Idx minf = std::numeric_limits<Idx>::min() / 2;
  using Val = std::pair<Idx, Idx>;
  using Lazy = int;
private:
  struct node{
    int h;
    Idx lx, rx, lmin, rmax;
    Idx val;
    Val sum;
    bool is_reset;
    Lazy lazy, lazy_sum;
    node *l, *r;
    node(Idx lx, Idx rx, Idx _val = 0): h(1), lx(lx), rx(rx), lmin(lx), rmax(rx), 
    val(_val), sum{_val, rx - lx}, is_reset(false), lazy(0), lazy_sum(0), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root;
  node *make_node(Idx l, Idx r, int x){
    return new node(l, r, x);
  }
  void merge(Val &a, const Val &b){
    if(a.first < b.first) return;
    else if(a.first > b.first) a = b;
    else a.second += b.second;
  }
  // vの区間がこれよりも大きい場合のみ使える
  std::pair<node*, node*> split_node(node *v, Idx l, Idx r){
    node *x = nullptr, *y = nullptr;
    if(v->lx < l){
      x = new node(v->lx, l, v->lazy_sum);
      x->lazy_sum = v->lazy_sum;
    }
    if(r < v->rx){
      y = new node(r, v->rx, v->lazy_sum);
      y->lazy_sum = v->lazy_sum;
    }
    v->lx = l, v->rx = r;
    return {x, y};
  }
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->lmin = v->lx;
    v->rmax = v->rx;
    v->sum.first = v->val;
    v->sum.second = v->rx - v->lx;
    if(v->l){
      v->lmin = v->l->lmin;
      merge(v->sum, v->l->sum);
    }
    if(v->r){
      v->rmax = v->r->rmax;
      merge(v->sum, v->r->sum);
    }
  }
  void propagate(node *v, Lazy x){
    v->lazy += x;
    v->lazy_sum += x;
    v->val += x;
    v->sum.first += x;
  }
  void reset(node *v){
    v->is_reset = true;
    v->lazy = v->lazy_sum = 0;
    v->val = 0;
    v->sum = std::make_pair(0, v->rmax - v->lmin);
  }
  void push_down(node *v){
    if(v->is_reset){
      if(v->l) reset(v->l);
      if(v->r) reset(v->r);
      v->is_reset = false;
    }
    if(v->lazy != 0){
      if(v->l) propagate(v->l, v->lazy);
      if(v->r) propagate(v->r, v->lazy);
      v->lazy = 0;
    }
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
  node *reset_inner(node *v, Idx l, Idx r){
    if(l == v->lmin && r == v->rmax){
      reset(v);
      return v;
    }
    push_down(v);
    if(l < v->lx){
      if(r <= v->lx){
        v->l = reset_inner(v->l, l, r);
        update(v);
        return balance(v);
      }
      v->l = reset_inner(v->l, l, v->lx);
      l = v->lx;
    }
    if(v->rx < r){
      if(v->rx <= l){
        v->r = reset_inner(v->r, l, r);
        update(v);
        return balance(v);
      }
      v->r = reset_inner(v->r, v->rx, r);
      r = v->rx;
    }
    assert(r > l);
    if(l != v->lx || r != v->rx){
      auto [_l, _r] = split_node(v, l, r);
      if(_l) v->l = insert_rightmost(v->l, _l);
      if(_r) v->r = insert_leftmost(v->r, _r);
    }
    v->val = 0;
    v->lazy_sum = 0;
    update(v);
    return balance(v);
  }
  node *update_inner(node *v, Idx l, Idx r, Lazy x){
    if(l == v->lmin && r == v->rmax){
      propagate(v, x);
      return v;
    }
    push_down(v);
    if(l < v->lx){
      if(r <= v->lx){
        v->l = update_inner(v->l, l, r, x);
        update(v);
        return balance(v);
      }
      v->l = update_inner(v->l, l, v->lx, x);
      l = v->lx;
    }
    if(v->rx < r){
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
    v->val += x;
    v->lazy_sum += x;
    update(v);
    return balance(v);
  }
  Val query_inner(node *v, Idx l, Idx r){
    if(l == v->lmin && r == v->rmax) return v->sum;
    push_down(v);
    Val left_sum = std::make_pair(inf, 0);
    if(l < v->lx){
      // 左側とだけ交差する
      if(r <= v->lx) return query_inner(v->l, l, r);
      left_sum = query_inner(v->l, l, v->lx);
      l = v->lx;
    }
    Val right_sum = std::make_pair(inf, 0);
    if(v->rx < r){
      // 右側とだけ交差する
      if(v->rx <= l) return query_inner(v->r, l, r);
      right_sum = query_inner(v->r, v->rx, r);
      r = v->rx;
    }
    Val ret = std::make_pair(v->val, r - l);
    merge(ret, left_sum);
    merge(ret, right_sum);
    return ret;
  }
  Idx find_right(node *v, Idx l){
    if(!v || v->rmax <= l) return minf;
    Idx m = v->sum.first;
    if(l <= v->lmin && m > 0) return minf;
    push_down(v);
    Idx x = find_right(v->l, l);
    if(x != minf) return x;
    if(l < v->rx){
      Idx lx = std::max(v->lx, l);
      Idx tmp = v->val;
      if(tmp == 0) return lx;
    }
    return find_right(v->r, l);
  }
  Idx find_left(node *v, Idx r){
    if(!v || r < v->lmin) return minf;
    Idx m = v->sum.first;
    if(r >= v->rmax - 1 && m > 0) return minf;
    push_down(v);
    Idx x = find_left(v->r, r);
    if(x != minf) return x;
    if(r >= v->rx - 1){
      Idx rx = std::min(v->rx, r + 1);
      Idx tmp = v->val;
      if(tmp == 0) return v->rx;
    }
    return find_left(v->l, r);
  }
public:
  segment_merge_super(): root(new node(minf, inf, 0)){}
  // 区間[l, r)を追加
  void merge(Idx l, Idx r){
    if(l >= r) return;
    root = update_inner(root, l, r, 1);
  }
  // [l, r)の部分を消す. 交差する区間がある場合[l, r)の部分だけ消える
  void erase(Idx l, Idx r){
    if(l >= r) return;
    root = reset_inner(root, l, r);
  }
  // 過去に行ったmergeを取り消す
  // 操作を行った部分がすでにeraseされていると壊れる
  void undo_merge(Idx l, Idx r){
    if(l >= r) return;
    root = update_inner(root, l, r, -1);
  }
  // [l, r)中の1つ以上の区間に含まれる部分の長さ
  Idx count(Idx l, Idx r){
    if(l >= r) return 0;
    auto [min_val, min_cnt] = query_inner(root, l, r);
    assert(min_val >= 0);
    return min_val == 0 ? r - l - min_cnt : r - l; 
  }
  // [l, r)中のいずれの区間にも含まれない部分の長さ
  Idx count_not(Idx l, Idx r){
    if(l >= r) return 0;
    auto [min_val, min_cnt] = query_inner(root, l, r);
    assert(min_val >= 0);
    return min_val == 0 ? min_cnt : 0;
  }
  // kから1つ以上の区間に含まれる部分をだけ通って移動できる範囲を[l, r)とする
  // kを含む区間がない場合はinf, ある場合はlを返す
  Idx find_left(Idx k){
    Idx ret = find_left(root, k);
    return ret > k ? inf : ret;
  }
  // kから1つ以上の区間に含まれる部分をだけ通って移動できる範囲を[l, r)とする
  // kを含む区間がない場合はinf, ある場合はrを返す
  Idx find_right(Idx k){
    Idx ret = find_right(root, k);
    return ret <= k ? inf : ret;
  }
  // a以上かつkがどの区間にも含まれないような最小のk
  Idx mex(Idx a = 0){
    Idx r = find_right(a);
    return r == inf ? a : r;
  }
  bool same(Idx a, Idx b){
    if(a > b) std::swap(a, b);
    Idx r = find_right(a);
    return r != inf && b < r;
  }
};
*/
#endif
