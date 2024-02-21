#ifndef _DYNAMIC_FUNCTIONAL_GRAPH_H_
#define _DYNAMIC_FUNCTIONAL_GRAPH_H_
#include <vector>
#include <cassert>
#include <queue>
#include <algorithm>
#include "../../graph/dynamic_graph/link_cut_tree_plain.hpp"
#include "../../graph/dynamic_graph/link_cut_tree_path.hpp"
#include "../../minior/convert_directed_functional_graph.hpp"

struct dynamic_functional_graph{
  link_cut_tree_plain<std::pair<int, int>> lct;
  using node = link_cut_tree_plain<std::pair<int, int>>::node;
  std::vector<node*> V;
private:
  void build(const std::vector<int> &next){
    int n = next.size();
    V.resize(n);
    for(int i = 0; i < n; i++) V[i] = lct.make_node({i, i});
    std::vector<int> roots;
    for(int i = 0; i < n; i++){
      int j = next[i];
      // 連結でない場合繋ぐ, すでに連結な場合iを根としてその辺を省く
      if(!lct.is_same(V[i], V[j])){
        lct._link(V[i], V[j]);
      }else{
        V[i]->label.second = j;
        roots.push_back(i);
      }
    }
    for(int r : roots) lct.evert(V[r]);
  }
  int __move(int i, uint64_t k){
    if(!k) return i;
    int d = lct.depth(V[i]);
    if(d >= k) return lct.la(V[i], k)->label.first;
    k -= d;
    auto [rid, cid] = lct.get_root(V[i])->label;
    node *c = V[cid];
    int len = lct.depth(c) + 1;
    k %= len;
    if(!k) return rid;
    return lct.la(c, k - 1)->label.first;
  }
  int __dist_undirected(int i, int j){
    int d1 = to_loop(i), d2 = to_loop(j);
    int ret = 0;
    if(d1 > d2){
      i = move(i, d1 - d2);
      ret += d1 - d2;
      d1 = d2;
    }else{
      j = move(j, d2 - d1);
      ret += d2 - d1;
    }
    if(i == j) return ret; // 同じ枝
    ret += 2 * d1;
    int diff = lct.depth(V[i]) - lct.depth(V[j]);
    int li = len_loop(i);
    if(diff < 0) diff += li;
    if(move(i, diff) != j) return -1; // 異なる連結成分
    return ret + std::min(diff, li - diff);
  }
  int __dist_directed(int i, int j){
    int d1 = to_loop(i), d2 = to_loop(j);
    if(d2 > d1) return -1;
    i = move(i, d1 - d2);
    if(i == j) return d1 - d2;
    else if(d2) return -1; // i, jが違う枝にいる場合
    int diff = lct.depth(V[i]) - lct.depth(V[j]);
    if(diff < 0) diff += len_loop(i);
    if(move(i, diff) != j) return -1; // 異なる連結成分
    return d1 + diff;
  }
  void __change_edge_undirected(int i, int j, int a, int b){
    node *r = lct.get_root(V[i]);
    if(__move(i, 1) != j) std::swap(i, j);
    auto [rid, cid] = r->label;
    // 辺を切る
    if(rid != i){
      node *c = V[cid];
      lct.cut_from_parent(V[i]);
      if(!lct.is_same(r, c)) lct._link(r, c);
    }
    // 辺を繋ぐ
    if(lct.is_same(V[a], V[b])){
      lct.evert(V[a]);
      V[a]->label.second = b;
    }else{
      // aかbのどちらかがiと同じ連結成分
      if(!lct.is_same(V[i], V[a])) std::swap(a, b);
      node *r2 = lct.get_root(V[b]);
      lct._link(V[i], V[b]);
      lct.evert(r2);
    }
  }
  void __change_edge_directed(int i, int j){
    node *r = lct.get_root(V[i]);
    auto [rid, cid] = r->label;
    // 辺を切る
    if(rid != i){
      node *c = V[cid];
      lct.cut_from_parent(V[i]);
      if(!lct.is_same(r, c)) lct._link(r, c);
    }
    // 辺を繋ぐ
    if(lct.is_same(V[i], V[j])){
      lct.evert(V[i]);
      V[i]->label.second = j;
    }else{
      node *r2 = lct.get_root(V[j]);
      lct._link(V[i], V[j]);
      lct.evert(r2);
    }
  }
public:
  dynamic_functional_graph(){}
  dynamic_functional_graph(const std::vector<int> &next){
    build(next);
  }
  // 同じ連結成分にあるか
  bool same(int i, int j){
    return lct.is_same(V[i], V[j]);
  }
  // iからk回遷移した状態
  int move(int i, uint64_t k){
    return __move(i, k);
  }
  // iから辿り着けるループの周期
  int len_loop(int i){
    auto [rid, cid] = lct.get_root(V[i])->label;
    return lct.depth(V[cid]) + 1;
  }
  // iからループまでの距離
  int to_loop(int i){
    auto [rid, cid] = lct.get_root(V[i])->label;
    node *l = lct.lca(V[cid], V[i]);
    return lct.depth(V[i]) - lct.depth(l);
  }
  // 無向functional_graphでのi -> jの距離
  // 異なる連結成分にいる場合は-1(これ以外の場合で-1にならない)
  int dist_undirected(int i, int j){
    return __dist_undirected(i, j);
  }
  // 有向functional_graphでのi -> jの距離 移動が不可能な場合-1
  int dist_directed(int i, int j){
    return __dist_directed(i, j);
  }
  // 無向functional_graphで辺(i, j)を消して辺　(a, b)を追加する
  // この操作を行った後で, functional_graphとして壊れている場合壊れる(消す辺がない場合, 木が発生する場合など)
  // 自己辺, 多重辺があってもいい
  void change_edge_undirected(int i, int j, int a, int b){
    __change_edge_undirected(i, j, a, b);
  }
  // 有向functional_graphでiの遷移先をjに変更する
  // 自己辺があってもいい
  void change_edge_directed(int i, int j){
    __change_edge_directed(i, j);
  }
};


template<typename monoid>
struct dynamic_functional_graph_monoid{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
  static constexpr auto flip = monoid::flip;

  using node = typename link_cut_tree_path<monoid, std::pair<int, int>>::node;
  link_cut_tree_path<monoid, std::pair<int, int>> lct;
  std::vector<node*> V;
private:
  // min, maxや+のときはO(1)にできる
  Val pow_monoid(Val val, uint64_t k){
    Val ret = id();
    while(k){
      if(k & 1) ret = merge(ret, val);
      val = merge(val, val);
      k >>= 1;
    }
    return ret;
  }
  // fがfalseのままxを最大何回マージできるか
  template<typename F>
  std::pair<uint64_t, Val> bisect_monoid(Val x, F f){
    if(f(x)) return 0;
    std::vector<Val> val{x};
    while(true){
      Val y = val.back();
      y = merge(y, y);
      if(f(y)) break;
      val.push_back(y);
    }
    uint64_t cnt = 0;
    Val sum = id();
    for(int i = (int)val.size() - 1; i >= 0; i--){
      Val tmp = merge(sum, val[i]);
      if(f(tmp)) continue;
      cnt += (uint64_t)1 << i;
      sum = tmp;
    }
    return {cnt, sum};
  }
  void build(const std::vector<int> &next, const std::vector<Val> &val){
    int n = next.size();
    V.resize(n);
    for(int i = 0; i < n; i++) V[i] = lct.make_node(val[i], {i, i});
    std::vector<int> roots;
    for(int i = 0; i < n; i++){
      int j = next[i];
      // 連結でない場合繋ぐ, すでに連結な場合iを根としてその辺を省く
      if(!lct.is_same(V[i], V[j])){
        lct._link(V[i], V[j]);
      }else{
        V[i]->label.second = j;
        roots.push_back(i);
      }
    }
    for(int r : roots) lct.evert(V[r]);
  }
  int __move(int i, uint64_t k){
    if(!k) return i;
    int d = lct.depth(V[i]);
    if(d >= k) return lct.la(V[i], k)->label.first;
    k -= d;
    auto [rid, cid] = lct.get_root(V[i])->label;
    node *c = V[cid];
    int len = lct.depth(c) + 1;
    k %= len;
    if(!k) return rid;
    return lct.la(c, k - 1)->label.first;
  }
  std::pair<int, Val> __move2(int i, uint64_t k){
    if(!k) return V[i]->val;
    int d = lct.depth(V[i]);
    node *r = lct.get_root(V[i]);
    if(d >= k){
      node *v = lct.la(V[i], k);
      lct.evert(v);
      Val ret = lct.query_path(V[i]);
      lct.evert(r);
      return {v->label.first, flip(ret)};
    }
    k -= d;
    node *c = V[r->label.second];
    int len = lct.depth(c) + 1;
    Val ret = lct.query_path(V[i]);
    if(k >= len){
      ret = merge(pow_monoid(lct.query_path(c), k / len), ret);
      k %= len;
    }
    if(!k) return {r->label.first, flip(ret)};
    node *v = lct.la(c, k - 1);
    lct.evert(v);
    ret = merge(lct.query_path(c), ret);
    lct.evert(r);
    return {v->label.first, flip(ret)};
  }
  int __dist_undirected(int i, int j){
    int d1 = to_loop(i), d2 = to_loop(j);
    int ret = 0;
    if(d1 > d2){
      i = move(i, d1 - d2);
      ret += d1 - d2;
      d1 = d2;
    }else{
      j = move(j, d2 - d1);
      ret += d2 - d1;
    }
    if(i == j) return ret; // 同じ枝
    ret += 2 * d1;
    int diff = lct.depth(V[i]) - lct.depth(V[j]);
    int li = len_loop(i);
    if(diff < 0) diff += li;
    if(move(i, diff) != j) return -1; // 異なる連結成分
    return ret + std::min(diff, li - diff);
  }
  int __dist_directed(int i, int j){
    int d1 = to_loop(i), d2 = to_loop(j);
    if(d2 > d1) return -1;
    i = move(i, d1 - d2);
    if(i == j) return d1 - d2;
    else if(d2) return -1; // i, jが違う枝にいる場合
    int diff = lct.depth(V[i]) - lct.depth(V[j]);
    if(diff < 0) diff += len_loop(i);
    if(move(i, diff) != j) return -1; // 異なる連結成分
    return d1 + diff;
  }
  void __change_edge_undirected(int i, int j, int a, int b){
    node *r = lct.get_root(V[i]);
    if(__move(i, 1) != j) std::swap(i, j);
    auto [rid, cid] = r->label;
    // 辺を切る
    if(rid != i){
      node *c = V[cid];
      lct.cut_from_parent(V[i]);
      if(!lct.is_same(r, c)) lct._link(r, c);
    }
    // 辺を繋ぐ
    if(lct.is_same(V[a], V[b])){
      lct.evert(V[a]);
      V[a]->label.second = b;
    }else{
      // aかbのどちらかがiと同じ連結成分
      if(!lct.is_same(V[i], V[a])) std::swap(a, b);
      node *r2 = lct.get_root(V[b]);
      lct._link(V[i], V[b]);
      lct.evert(r2);
    }
  }
  void __change_edge_directed(int i, int j){
    node *r = lct.get_root(V[i]);
    auto [rid, cid] = r->label;
    // 辺を切る
    if(rid != i){
      node *c = V[cid];
      lct.cut_from_parent(V[i]);
      if(!lct.is_same(r, c)) lct._link(r, c);
    }
    // 辺を繋ぐ
    if(lct.is_same(V[i], V[j])){
      lct.evert(V[i]);
      V[i]->label.second = j;
    }else{
      node *r2 = lct.get_root(V[j]);
      lct._link(V[i], V[j]);
      lct.evert(r2);
    }
  }
  template<typename F>
  std::tuple<int, uint64_t, Val> __bisect(int i, F &f){
    Val s = flip(lct.query_path(V[i]));
    uint64_t d = V[i]->sz;
    node *r = lct.get_root(V[i]);
    if(f(s)){
      lct.evert(V[i]);
      node *v = lct.expose(lct.bisect_from_root(r, f));
      std::tuple<int, uint64_t, Val> ret = {v->label.first, v->sz - 1, v->sum};
      lct.evert(r);
      return ret;
    }
    node *c = V[r->label.second];
    Val loop_sum = flip(lct.query_path(c));
    uint64_t loop_len = c->sz;
    // 初めて根にたどり着くまでの{距離, 積} = {d, s}
    // ループ１周あたりの　{loop_len, loop_sum}
    auto [loop_times, sum] = bisect_monoid(loop_sum, [&](Val x){return f(merge(s, x));});
    s = merge(s, sum);
    d += loop_len * loop_times;
    lct.evert(c);
    node *v = lct.expose(lct.bisect_from_root(r, [&](Val x){return f(merge(s, x));}));
    std::tuple<int, uint64_t, Val> ret = {v->label.first, d + (v->sz - 1), merge(s, v->sum)};
    lct.evert(r);
    return ret;
  }
public:
  dynamic_functional_graph_monoid(){}
  
  dynamic_functional_graph_monoid(const std::vector<int> &next, const std::vector<Val> &val){
    build(next, val);
  }
  // 同じ連結成分にあるか
  bool same(int i, int j){
    return lct.is_same(V[i], V[j]);
  }
  // iからk回遷移した状態
  int move(int i, uint64_t k){
    return __move(i, k);
  }
  // {iからk回遷移した状態, 積} (始点を含む)
  std::pair<int, Val> move2(int i, uint64_t k){
    return __move2(i, k);
  }
  // iの値をxにする
  void set(int i, Val x){
    lct.set(V[i], x);
  }
  // iの値
  Val get(int i){
    return lct.get(V[i]);
  }
  // iから初めて値の積が初めてtrueになるような{頂点, 移動回数, そのときの値} (始点を含む)
  template<typename F>
  std::tuple<int, uint64_t, Val> bisect(int i, F f){
    return __bisect(i, f);
  }
  // iから辿り着けるループの周期
  int len_loop(int i){
    auto [rid, cid] = lct.get_root(V[i])->label;
    return lct.depth(V[cid]) + 1;
  }
  // iからループまでの距離
  int to_loop(int i){
    auto [rid, cid] = lct.get_root(V[i])->label;
    node *l = lct.lca(V[cid], V[i]);
    return lct.depth(V[i]) - lct.depth(l);
  }
  // 無向functional_graphでのi -> jの距離
  // 異なる連結成分にいる場合は-1(これ以外の場合で-1にならない)
  int dist_undirected(int i, int j){
    return __dist_undirected(i, j);
  }
  // 有向functional_graphでのi -> jの距離 移動が不可能な場合-1
  int dist_directed(int i, int j){
    return __dist_directed(i, j);
  }
  // 無向functional_graphで辺(i, j)を消して辺　(a, b)を追加する
  // この操作を行った後で, functional_graphとして壊れている場合壊れる(消す辺がない場合, 木が発生する場合など)
  // 自己辺, 多重辺があってもいい
  void change_edge_undirected(int i, int j, int a, int b){
    __change_edge_undirected(i, j, a, b);
  }
  // 有向functional_graphでiの遷移先をjに変更する
  // 自己辺があってもいい
  void change_edge_directed(int i, int j){
    __change_edge_directed(i, j);
  }
};
#endif