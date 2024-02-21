#ifndef _UNION_FIND_MONOID_H_
#define _UNION_FIND_MONOID_H_
#include <vector>
#include <cassert>
#include <queue>
#include <unordered_map>
#include "../../minior/red_black_tree_beats.hpp"
template<typename beats_struct>
struct union_find_beats{
private:
  using rbt = red_black_tree_beats<beats_struct>;
  using Val = typename beats_struct::Val;
  using Lazy = typename beats_struct::Lazy;
  using node = typename rbt::node;
  std::vector<node*> V;
  std::unordered_map<node*, int> rev;
  std::vector<std::tuple<bool, int, int>> history; // {有効か, 親, 子のサイズ}
  node *__find(int a){
    node *v = V[a];
    while(v->p) v = v->p;
    return v;
  }
public:
  union_find_beats(int n): V(n){
    for(int i = 0; i < n; i++){
      V[i] = new node(beats_struct::id_val());
      rev.emplace(V[i], i);
    }
  }
  template<typename T>
  union_find_beats(const std::vector<T> &val): V(val.size()){
    for(int i = 0; i < val.size(); i++){
      V[i] = new node(Val(val[i]));
      rev.emplace(V[i], i);
    }
  }
  void unite(int a, int b){
    node *u = __find(a);
    node *v = __find(b);
    if(u == v){
      history.push_back({false, a, 0});
      return;
    }
    history.push_back({true, a, v->sz});
    rbt::merge(u, v);
  }
  // 通常のunion-find木ではunite(a, b)後の根はaの根, bの根のうちサイズが大きい方であるが, このfindはそのルールを守らない
  // 同じ木の状態で同じ連結成分の要素から呼ばれると同じ値[0, n)を返すことのみ保証される
  int find(int a){
    node *v = V[a];
    while(v->p) v = v->p;
    while(v->l) v = v->l;
    return rev[v];
  }
  bool same(int a, int b){
    return __find(a) == __find(b);
  }
  int size(int a){
    return __find(a)->sz;
  }
  template<int id>
  void update(int a, Lazy x){
    node *v = __find(a);
    rbt::propagate(v, x);
  }
  Val query(int a){
    return __find(a)->val;
  }
  Val get(int a){
    std::vector<node*> path;
    node *v = V[a];
    while(v->p){
      v = v->p;
      path.push_back(v);
    }
    int m = path.size();
    for(int i = m - 1; i >= 0; i--) rbt::push_down(path[i]);
    return V[a]->val;
  }
  template<typename T>
  void set(int a, T x){
    std::vector<node*> path;
    node *v = V[a];
    while(v->p){
      v = v->p;
      path.push_back(v);
    }
    int m = path.size();
    for(int i = m - 1; i >= 0; i--) rbt::push_down(path[i]);
    V[a]->val = Val(x);
    for(int i = 0; i < m; i++) beats_struct::merge_val(path[i]->val, path[i]->l->val, path[i]->r->val);
  }
  void rollback(){
    assert(history.size() > 0);
    auto [q, p, s] = history.back();
    history.pop_back();
    if(!q) return;
    node *v = __find(p);
    rbt::split(v, s);
  }
  std::vector<int> enumerate(int a){
    std::vector<int> res;
    node *v = V[a];
    while(v->p) v = v->p;
    std::queue<node*> q;
    q.push(v);
    while(!q.empty()){
      v = q.front();
      q.pop();
      if(v->l){
        q.push(v->l);
        q.push(v->r);
      }else{
        res.push_back(rev[v]);
      }
    }
    return res;
  }
};
#endif
