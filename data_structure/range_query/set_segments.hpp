#ifndef _SET_SEGMENTS_H_
#define _SET_SEGMENTS_H_

#include "../basic/leftist_heap.hpp"
template<typename Idx, bool merge_adjacent = true>
struct heap_segments{
private:
  leftist_heap<std::pair<Idx, Idx>> h;
  // 左端が最小の区間をマージできる限りマージ
  void __modify(){
    assert(!h.empty());
    auto [l, r] = h.pop_min();
    while(!h.empty() && h.min().first + (!merge_adjacent) <= r){
      r = std::max(r, h.pop_min().second);
    }
    h.push({l, r});
  }
public:
  bool empty(){
    return h.empty();
  }
  // [l, r)をマージ
  // merge_adjacent: [1, 2), [2, 5)のような隣接する区間をマージするか
  void merge(Idx l, Idx r){
    h.push({l, r});
  }
  // 左端が最小の区間　
  std::pair<Idx, Idx> min(){
    __modify();
    return h.min();
  }
  // 左端が最小の区間をpopして返す
  std::pair<Idx, Idx> pop_min(){
    __modify();
    return h.pop_min();
  }
  // 任意の区間に含まれない0以上の最小要素(merge_adjacent = trueのときだけ使える)
  Idx mex(){
    static_assert(merge_adjacent);
    return empty() ? 0 : min().second;
  }
  // rの要素を全て移動(永続でないためrの要素は全て消える)
  void meld(heap_segments &r){
    h.meld(r.h);
  }
};

template<typename Idx, bool merge_adjacent = true>
struct set_segments{
  static constexpr Idx minf = std::numeric_limits<Idx>::min();
  static constexpr Idx inf = std::numeric_limits<Idx>::max();
private:
  struct node{
    int h, sz;
    Idx L, R, lensum;
    node *l, *r;
    node(Idx _L, Idx _R): h(1), sz(1), L(_L), R(_R), lensum(R - L), l(nullptr), r(nullptr){}
    int balanace_factor(){
      return (l ? l->h : 0) - (r ? r->h : 0);
    }
  };
  node *root, *tmp_node;
  int size(node *v){return v ? v->sz : 0;}
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    v->sz = (v->l ? v->l->sz : 0) + (v->r ? v->r->sz : 0) + 1;
    v->lensum = (v->R - v->L) + (v->l ? v->l->lensum : 0) + (v->r ? v->r->lensum : 0);
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
  node *leftmost(node *v){
    while(v->l) v = v->l;
    return v;
  }
  node *rightmost(node *v){
    while(v->r) v = v->r;
    return v;
  }
  std::tuple<node*, Idx, Idx> __find(node *v, Idx k){
    Idx Lmax = minf, Rmin = inf;
    while(v){
      if(v->L <= k && k < v->R){
        return {v, v->L, v->R};
      }else if(k < v->L){
        Rmin = v->L;
        v = v->l;
      }else{
        Lmax = v->R;
        v = v->r;
      }
    }
    return {nullptr, Lmax, Rmin};
  }
  Idx __kth_point(node *v, Idx k){
    while(true){
      Idx lenl = (v->l ? v->l->lensum : 0);
      Idx lenv = lenl + (v->R - v->L);
      if(lenl <= k){
        if(k < lenv) return v->L + (k - lenl);
        k -= lenv;
        v = v->r;
      }else v = v->l;
    }
    return inf;
  }
  node *__kth_segment(node *v, int k){
    while(true){
      int szl = size(v->l);
      if(szl <= k){
        if(szl == k) return v;
        k -= szl + 1;
        v = v->r;
      }else v = v->l;
    }
  }
  int __low_count(node *v, Idx x){
    int res = 0;
    while(v){
      int szl = size(v->l);
      if(x < v->R) v = v->l;
      else v = v->r, res += szl + 1;
    }
    return res;
  }
  Idx __low_count_sum(node *v, Idx x){
    Idx res = 0;
    while(v){
      if(x <= v->L){
        v = v->l;
      }else if(v->R <= x){
        res += (v->l ? v->l->lensum : 0) + (v->R - v->L);
        v = v->r;
      }else{
        return res + (v->l ? v->l->lensum : 0) + (x - v->L);
      }
    }
    return res;
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
  node *cut_rightmost(node *v){
    if(v->r){
      v->r = cut_rightmost(v->r);
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->l;
  }
  node *__insert(node *v, Idx l, Idx r){
    if(!v) return new node(l, r);
    if(l < v->L){
      v->l = __insert(v->l, l, r);
    }else{
      v->r = __insert(v->r, l, r);
    }
    update(v);
    return balance(v);
  }
  node *__erase(node *v, Idx l){
    if(!v) return nullptr;
    if(l < v->L){
      v->l = __erase(v->l, l);
    }else if(l > v->L){
      v->r = __erase(v->r, l);
    }else{
      if(v->r){
        v->r = cut_leftmost(v->r);
        tmp_node->l = v->l;
        tmp_node->r = v->r;
        free(v);
        update(tmp_node);
        return balance(tmp_node);
      }
      return v->l;
    }
    update(v);
    return balance(v);
  }
  void __merge(Idx l, Idx r){
    auto [L, R] = __erase_intersect(l - merge_adjacent, r + merge_adjacent);
    root = __insert(root, std::min(L, l), std::max(R, r));
  }
  // 消された中で{最小, 最大}
  std::pair<Idx, Idx> __erase_include(Idx l, Idx r){
    Idx emin = inf, emax = minf;
    for(auto [L, R] : enumerate_include(l, r)){
      root = __erase(root, L);
      emin = std::min(emin, L);
      emax = std::max(emax, R);
    }
    return {emin, emax};
  }
  // 消された中で{最小, 最大}
  std::pair<Idx, Idx> __erase_intersect(Idx l, Idx r){
    Idx emin = inf, emax = minf;
    for(auto [L, R] : enumerate_intersect(l, r)){
      root = __erase(root, L);
      emin = std::min(emin, L);
      emax = std::max(emax, R);
    }
    return {emin, emax};
  }
  void __enumerate_include(node *v, Idx l, Idx r, std::vector<std::pair<Idx, Idx>> &res){
    if(!v) return;
    if(v->l && l < v->L) __enumerate_include(v->l, l, r, res);
    if(l <= v->L && v->R <= r) res.push_back({v->L, v->R});
    if(v->r && v->R < r) __enumerate_include(v->r, l, r, res);
  }
  void __enumerate_intersect(node *v, Idx l, Idx r, std::vector<std::pair<Idx, Idx>> &res){
    if(!v) return;
    if(v->l && l < v->L) __enumerate_intersect(v->l, l, r, res);
    if(std::max(l, v->L) < std::min(r, v->R)) res.push_back({v->L, v->R});
    if(v->r && v->R < r) __enumerate_intersect(v->r, l, r, res);
  }
public:
  set_segments(): root(nullptr){}
  int size(){
    return size(root);
  }
  bool empty(){
    return size_sum(root) == 0;
  }
  // a, bが同じ区間に含まれるか
  bool same(Idx a, Idx b){
    auto [v, l, r] = find(a);
    return v && (l == std::get<1>(find(b)));
  }
  // kがいずれかの区間に含まれる場合 {1, L, R}
  // 含まれない場合 {0, L, R} (L, Rはkからいずれの区間にも含まれない座標だけを通って移動できる範囲, 何もない場合は{minf, inf})
  std::tuple<bool, Idx, Idx> find(Idx k){
    auto [v, L, R] = __find(root, k);
    return v ? std::make_tuple(true, L, R) : std::make_tuple(false, L, R);
  }
  std::pair<Idx, Idx> min(){
    assert(size());
    node *v = leftmost(root);
    return {v->L, v->R};
  }
  std::pair<Idx, Idx> max(){
    assert(size());
    node *v = rightmost(root);
    return {v->L, v->R};
  }
  // 任意の区間に含まれないかつa以上の最小要素
  Idx mex(Idx a = 0){
    static_assert(merge_adjacent);
    auto [v, L, R] = find(a);
    return v ? R : a;
  }
  // いずれかの区間に含まれるk番目(0-indexed)に小さい点. ない場合はinf
  Idx kth_point(Idx k){
    return __kth_point(root, k);
  }
  // k番目(0-indexed)に小さい区間. ない場合は{inf, inf}
  std::pair<Idx, Idx> kth_segment(int k){
    if(size() <= k) return {inf, inf};
    node *v = __kth_segment(root, k);
    return {v->L, v->R};
  }
  // [l, r)がr <= xであるような区間の数
  int low_count(Idx x){
    return __low_count(root, x);
  }
  // いずれかの区間に含まれ, かつx未満の座標の数
  Idx low_count_sum(Idx x){
    return __low_count_sum(root, x);
  }
  // [l, r)をマージ
  // merge_adjacent: [1, 2), [2, 5)のような隣接する区間をマージするか
  void merge(Idx l, Idx r){
    __merge(l, r);
  }
  // [l, r)に含まれる区間を削除
  void erase_include(Idx l, Idx r){
    __erase_include(l, r);
  }
  // [l, r)と少しでも交差する区間を削除
  void erase_intersect(Idx l, Idx r){
    __erase_intersect(l, r);
  }
  // [l, r)が完全に含む区間を列挙
  std::vector<std::pair<Idx, Idx>> enumerate_include(Idx l, Idx r){
    std::vector<std::pair<Idx, Idx>> res;
    __enumerate_include(root, l, r, res);
    return res;
  }
  // [l, r)と交差する区間を列挙
  std::vector<std::pair<Idx, Idx>> enumerate_intersect(Idx l, Idx r){
    std::vector<std::pair<Idx, Idx>> res;
    __enumerate_intersect(root, l, r, res);
    return res;
  }
  // 全区間を列挙
  std::vector<std::pair<Idx, Idx>> enumerate_all(){
    return enumerate_intersect(minf, inf);
  }
  void clear(){
    root = nullptr;
  }
  void swap(set_segments<Idx> &r){
    std::swap(root, r.root);
  }
  // rの要素を全て移動(永続でないためrの要素は全て消える)
  void meld(set_segments<Idx> &r){
    if(size() < r.size()) swap(r);
    for(auto [L, R] : r.enumerate_all()) merge(L, R);
    r.clear();
  }
};
#endif