#ifndef _ONLINE_DYNAMIC_CONNECTIVITY_MONOID_H_
#define _ONLINE_DYNAMIC_CONNECTIVITY_MONOID_H_
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <unordered_set>
#include "../../algebraic_structure/monoid.hpp"

namespace online_dynamic_connectivity_monoid_internal{
  struct surplus_edge{
    int N, maxlv;
    std::vector<std::vector<std::unordered_multiset<int>>> E;
    surplus_edge(){}
    surplus_edge(int N): N(N), maxlv(1), E(1, std::vector<std::unordered_multiset<int>>(N)){}
    void make_new_level(){
      maxlv++;
      E.push_back(std::vector<std::unordered_multiset<int>>(N));
    }
    void insert(int lv, int s, int t){
      E[lv][s].insert(t);
      E[lv][t].insert(s);
    }
    bool erase(int lv, int s, int t){
      auto itr = E[lv][s].find(t);
      if(itr == E[lv][s].end()) return false;
      E[lv][s].erase(itr);
      E[lv][t].erase(E[lv][t].find(s));
      return true;
    }
    int erase_any(int lv, int s){
      auto itr = E[lv][s].begin();
      int ret;
      if(itr == E[lv][s].end()) ret = -1;
      else{
        ret = *itr;
        E[lv][s].erase(itr);
        E[lv][ret].erase(E[lv][ret].find(s));
      }
      return ret;
    }
    bool empty(int lv, int s){
      return E[lv][s].empty();
    }
  };
  template<typename monoid, typename Val, typename Lazy>
  struct toptree{
    struct node{
      node *l, *r, *p, *light_top, *l2, *r2, *p2;
      Val val, sum, light_sum;
      Lazy lazy;
      int size, light_size;
      bool flip, is_x;
      int id;
      std::unordered_set<int> t;
      node *heavy_x, *heavy_y, *light_x, *light_y;
      node(Val val, int id): l(nullptr), r(nullptr), p(nullptr), light_top(nullptr), l2(nullptr), r2(nullptr), p2(nullptr),
      val(val), sum(val), light_sum(val), lazy(monoid::template id2<Lazy>()), size(1), light_size(1), flip(false), is_x(false), 
      id(id), heavy_x(nullptr), heavy_y(nullptr), light_x(nullptr), light_y(nullptr){}
      bool is_root_heavy(){return !p || (p->l != this && p->r != this);}
      bool is_root_light(){return !p2 || (p2->l2 != this && p2->r2 != this);}
    };
    node *make_node(Val val, int id){return new node(val, id);}
    toptree(){}
    
    void update_heavy(node *v){
      v->size = 1;
      v->sum = v->val;
      if(v->l){
        v->size += v->l->size;
        v->sum = monoid::template merge<Val>(v->sum, v->l->sum);
      }
      if(v->r){
        v->size += v->r->size;
        v->sum = monoid::template merge<Val>(v->sum, v->r->sum);
      }
      if(v->light_top){
        v->size += v->light_top->light_size;
        v->sum = monoid::template merge<Val>(v->sum, v->light_top->light_sum);
      }
      v->heavy_x = v->is_x ? v : nullptr;
      v->heavy_y = !v->t.empty() ? v : nullptr;
      if(!v->heavy_x){
        if(v->l && v->l->heavy_x) v->heavy_x = v->l->heavy_x;
        else if(v->r && v->r->heavy_x) v->heavy_x = v->r->heavy_x;
        else if(v->light_top && v->light_top->light_x) v->heavy_x = v->light_top->light_x;
      }
      if(!v->heavy_y){
        if(v->l && v->l->heavy_y) v->heavy_y = v->l->heavy_y;
        else if(v->r && v->r->heavy_y) v->heavy_y = v->r->heavy_y;
        else if(v->light_top && v->light_top->light_y) v->heavy_y = v->light_top->light_y;
      }
    }
    void update_light(node *v){
      v->light_size = v->size;
      v->light_sum = v->sum;
      v->light_x = v->heavy_x;
      v->light_y = v->heavy_y;
      if(v->l2){
        v->light_size += v->l2->light_size;
        v->light_sum = monoid::template merge<Val>(v->light_sum, v->l2->light_sum);
        if(v->l2->light_x) v->light_x = v->l2->light_x;
        if(v->l2->light_y) v->light_y = v->l2->light_y;
      }
      if(v->r2){
        v->light_size += v->r2->light_size;
        v->light_sum = monoid::template merge<Val>(v->light_sum, v->r2->light_sum);
        if(v->r2->light_x) v->light_x = v->r2->light_x;
        if(v->r2->light_y) v->light_y = v->r2->light_y;
      }
    }
    void push_down(node *v){
      if(v->lazy != monoid::template id2<Lazy>()){
        if(v->l) propagate(v->l, v->lazy);
        if(v->r) propagate(v->r, v->lazy);
        if(v->light_top) propagate(v->light_top, v->lazy);
        if(v->l2) propagate(v->l2, v->lazy);
        if(v->r2) propagate(v->r2, v->lazy);
        v->lazy = monoid::template id2<Lazy>();
      }
      if(v->flip){
        if(v->l) flip(v->l);
        if(v->r) flip(v->r);
        v->flip = false;
      }
    }
    void flip(node *v){
      std::swap(v->l, v->r);
      v->flip ^= 1;
    }
    static void propagate(node *v, Lazy x){
      v->lazy = monoid::template propagate<Lazy>(v->lazy, x);
      v->val = monoid::template apply<Val, Lazy>(v->val, x, 0, 1);
      v->sum = monoid::template apply<Val, Lazy>(v->sum, x, 0, v->size);
      v->light_sum = monoid::template apply<Val, Lazy>(v->light_sum, x, 0, v->light_size);
    }
    void rotate_right_heavy(node *v){
      node *p = v->p, *pp = p->p;
      if((p->l = v->r)) v->r->p = p;
      v->r = p, p->p = v;
      update_heavy(p), update_heavy(v);
      if((v->p = pp)){
        if(pp->l == p) pp->l = v;
        if(pp->r == p) pp->r = v;
        update_heavy(pp);
      }
    }
    void rotate_left_heavy(node *v){
      node *p = v->p, *pp = p->p;
      if((p->r = v->l)) v->l->p = p;
      v->l = p, p->p = v;
      update_heavy(p), update_heavy(v);
      if((v->p = pp)){
        if(pp->l == p) pp->l = v;
        if(pp->r == p) pp->r = v;
        update_heavy(pp);
      }
    }
    void splay_heavy(node *v){
      node *u = nullptr;
      push_down(v);
      while(!v->is_root_heavy()){
        node *p = v->p;
        if(p->is_root_heavy()){
          u = p;
          splay_light(u);
          push_down(p), push_down(v);
          if(p->l == v) rotate_right_heavy(v);
          else rotate_left_heavy(v);
        }else{
          node *pp = p->p;
          if(pp->is_root_heavy()) u = pp, splay_light(u);
          push_down(pp), push_down(p), push_down(v);
          if(pp->l == p){
            if(p->l == v) rotate_right_heavy(p);
            else rotate_left_heavy(v);
            rotate_right_heavy(v);
          }else{
            if(p->r == v) rotate_left_heavy(p);
            else rotate_right_heavy(v);
            rotate_left_heavy(v);
          }
        }
      }
      if(u){
        node *l2 = u->l2, *r2 = u->r2, *p2 = u->p2;
        u->l2 = u->r2 = u->p2 = nullptr;
        if((v->l2 = l2)) l2->p2 = v;
        if((v->r2 = r2)) r2->p2 = v;
        update_light(v);
        if((v->p2 = p2)){
          p2->light_top = v;
          update_heavy(p2);
        }
      }
    }
    void rotate_right_light(node *v){
      node *p = v->p2, *pp = p->p2;
      if((p->l2 = v->r2)) v->r2->p2 = p;
      v->r2 = p, p->p2 = v;
      update_light(p), update_light(v);
      if((v->p2 = pp)){
        if(pp->l2 == p) pp->l2 = v;
        if(pp->r2 == p) pp->r2 = v;
        update_light(pp);
      }
    }
    void rotate_left_light(node *v){
      node *p = v->p2, *pp = p->p2;
      if((p->r2 = v->l2)) v->l2->p2 = p;
      v->l2 = p, p->p2 = v;
      update_light(p), update_light(v);
      if((v->p2 = pp)){
        if(pp->l2 == p) pp->l2 = v;
        if(pp->r2 == p) pp->r2 = v;
        update_light(pp);
      }
    }
    void splay_light(node *v){
      push_down(v);
      while(!v->is_root_light()){
        node *p = v->p2;
        if(p->is_root_light()){
          push_down(p), push_down(v);
          if(p->l2 == v) rotate_right_light(v);
          else rotate_left_light(v);
        }else{
          node *pp = p->p2;
          push_down(pp), push_down(p), push_down(v);
          if(pp->l2 == p){
            if(p->l2 == v) rotate_right_light(p);
            else rotate_left_light(v);
            rotate_right_light(v);
          }else{
            if(p->r2 == v) rotate_left_light(p);
            else rotate_right_light(v);
            rotate_left_light(v);
          }
        }
      }
      if(v->p2){
        v->p2->light_top = v;
        update_heavy(v->p2);
      }
    }
    void insert_light(node *p, node *c){
      push_down(p);
      push_down(c);
      if((c->l2 = p->light_top)) p->light_top->p2 = c;
      update_light(c);
      p->light_top = c;
      c->p2 = p;
      update_heavy(p);
    }
    void erase_light(node *p, node *c){
      splay_light(c);
      node *l = c->l2, *r = c->r2;
      c->l2 = c->r2 = c->p2 = nullptr;
      if(l && r){
        l->p2 = r->p2 = nullptr;
        while(l->r2) l = l->r2;
        splay_light(l);
        l->r2 = r;
        r->p2 = l;
        update_light(l);
      }else if(r) l = r;
      if(l) l->p2 = p;
      p->light_top = l;
      update_heavy(p);
    }
    void swap_light(node *p, node *a, node *b){
      push_down(a);
      splay_light(b);
      node *l = b->l2, *r = b->r2;
      b->l2 = b->r2 = b->p2 = nullptr;
      if((a->l2 = l)) l->p2 = a;
      if((a->r2 = r)) r->p2 = a;
      if((a->p2 = p)) p->light_top = a;
      update_light(a);
    }
    node *expose(node *v){
      node *c = nullptr;
      for(node *u = v; u; u = u->p){
        splay_heavy(u);
        if(c && u->r) swap_light(u, u->r, c);
        else if(c) erase_light(u, c);
        else if(u->r) insert_light(u, u->r);
        u->r = c;
        update_heavy(u);
        c = u;
      }
      splay_heavy(v);
      return c;
    }
    // vを根にする
    node *evert(node *v){
      expose(v);
      flip(v);
      push_down(v);
      return v;
    }
    // 非連結であることが保証されている場合
    void _link(node *p, node *c, bool is_k){
      evert(c);
      expose(p);
      c->p = p;
      p->r = c;
      if(is_k){
        if(p->id < c->id) p->t.insert(c->id);
        else{
          c->t.insert(p->id);
          update_heavy(c);
        }
      }
      update_heavy(p);
    }
    // cの親との辺を切る, false: 辺を切れなかった
    bool _cut(node *c){
      expose(c);
      node *p = c->l;
      if(p == nullptr) return false;
      c->l = p->p = nullptr;
      update_heavy(c);
      return true;
    }
    node *right_most(node *v){
      push_down(v);
      while(v->r){
        v = v->r;
        push_down(v);
      }
      return v;
    }
    bool cut(node *u, node *v){
      if(!u || !v) return false;
      evert(u);
      expose(v);
      if(!v->l) return false;
      node *rm = right_most(v->l);
      if(rm != u){
        splay_heavy(rm);
        return false;
      }
      if(v->id < u->id) v->t.erase(u->id);
      else u->t.erase(v->id);
      v->l->p = nullptr;
      v->l = nullptr;
      update_heavy(v);
      splay_heavy(u);
      return true;
    }
    int size(node *v){
      expose(v);
      return v->size;
    }
    node *get_root(node *v){
      expose(v);
      while(v->l){
        push_down(v);
        v = v->l;
      }
      splay_heavy(v);
      return v;
    }
    bool is_same(node *u, node *v){
      if(!u || !v) return false;
      return get_root(u) == get_root(v);
    }
    std::vector<std::pair<int, int>> enumerate_level_k_edge(node *v){
      std::vector<std::pair<int, int>> ret;
      expose(v);
      while(v->heavy_y){
        v = v->heavy_y;
        expose(v);
        for(int t: v->t) ret.push_back({v->id, t});
        v->t.clear();
        update_heavy(v);
      }
      return ret;
    }
    void set(node *v, Val x){
      expose(v);
      v->val = x;
      update_heavy(v);
    }
    Val get(node *v){
      expose(v);
      return v->val;
    }
    // vを含む連結成分にxを作用させる
    void update_cc(node *v, Lazy x){
      expose(v);
      propagate(v, x);
    }
    // vを含む連結成分のモノイド積　
    Val query_cc(node *v){
      expose(v);
      return v->sum;
    }
  };
};

template<typename monoid, typename Val, typename Lazy>
struct online_dynamic_connectivity_monoid{
private:
  using tt = online_dynamic_connectivity_monoid_internal::toptree<monoid, Val, Lazy>;
  using node = typename tt::node;
  int N, maxlv;
  online_dynamic_connectivity_monoid_internal::surplus_edge E;
  tt F;
  std::vector<std::vector<node*>> ptr;
  void make_new_level(){
    ptr.push_back(std::vector<node*>(N));
    for(int i = 0; i < N; i++ ) ptr[maxlv][i] = F.make_node(monoid::template id1<Val>(), i);
    E.make_new_level();
    maxlv++;
  }
  void on(node *v){
    F.expose(v);
    v->is_x = true;
    F.update_heavy(v);
  }
  void off(node *v){
    F.expose(v);
    v->is_x = false;
    F.update_heavy(v);
  }
  // a, bを連結にできるか判定
  bool replace(int a, int b, int k){
    if(k == -1) return false;
    if(F.size(ptr[k][a]) > F.size(ptr[k][b])) std::swap(a, b);
    node *A = ptr[k][a];
    auto e = F.enumerate_level_k_edge(A);
    for(auto [x, y] : e){
      if(k == maxlv - 1) make_new_level();
      F._link(ptr[k + 1][x], ptr[k + 1][y], true);
    }
    F.expose(A);
    while(A->heavy_x){
      A = A->heavy_x;
      if(k == maxlv - 1) make_new_level();
      node *A2 = ptr[k + 1][A->id];
      while(!E.empty(k, A->id)){
        int t = E.erase_any(k, A->id);
        if(E.empty(k, t)) off(ptr[k][t]);
        if(!F.is_same(A, ptr[k][t])){
          for(int i = 0; i <= k; i++){
            F._link(ptr[i][A->id], ptr[i][t], i == k);
          }
          return true;
        }
        if(E.empty(k + 1, A->id)) on(A2);
        if(E.empty(k + 1, t)) on(ptr[k + 1][t]);
        E.insert(k + 1, A->id, t);
      }
      off(A);
    }
    return replace(a, b, k - 1);
  }
public:
  int num_cc;

  online_dynamic_connectivity_monoid(int N): N(N), maxlv(1), E(N), ptr(1, std::vector<node*>(N)), num_cc(N){
    for(int i = 0; i < N; i++) ptr[0][i] = F.make_node(monoid::template id1<Val>(), i);
  }
  online_dynamic_connectivity_monoid(const std::vector<Val> &v): N(v.size()), maxlv(1), E(N), ptr(1, std::vector<node*>(N)), num_cc(N){
    for(int i = 0; i < N; i++) ptr[0][i] = F.make_node(v[i], i);
  }
  // 辺(a, b)を追加する　, すでに辺(a, b)があっても追加する
  // 0: すでに連結だった, 1: 新たに連結になった
  bool link(int a, int b){
    node *A = ptr[0][a], *B = ptr[0][b];
    if(F.is_same(A, B)){
      if(E.empty(0, a)) on(A);
      if(E.empty(0, b)) on(B);
      E.insert(0, a, b);
      return false;
    }else{
      F._link(A, B, true);
      num_cc--;
      return true;
    }
  }
  // 辺(a, b)を切る 辺が無かった: 0, 辺が橋だった: 1, まだ連結: 2
  int cut(int a, int b){
    node *A = ptr[0][a], *B = ptr[0][b];
    if(!F.is_same(A, B)) return 0;
    if(!F.cut(A, B)){
      for(int i = 0; i < maxlv; i++){
        if(E.erase(i, a, b)){
          if(E.empty(i, a)) off(ptr[i][a]);
          if(E.empty(i, b)) off(ptr[i][b]);
          return 2;
        }
      }
      return 0;
    }else{
      int k = 0;
      for(int i = 1; i < maxlv; i++){
        if(F.cut(ptr[i][a], ptr[i][b])) k = i;
        else break;
      }
      int ret = replace(a, b, k) ? 2 : 1;
      if(ret == 1) num_cc++;
      return ret;
    }
  }
  // aを含む連結成分のサイズ
  int size(int a){
    return F.size(ptr[0][a]);
  }
  // a, bが連結か
  bool same(int a, int b){
    return F.is_same(ptr[0][a], ptr[0][b]);
  }
  // 連結成分の数
  int connected_component(){
    return num_cc;
  }
  // aの値をxにする
  void set(int a, Val x){
    F.set(ptr[0][a], x);
  }
  // aの値を取得
  Val get(int a){
    return F.get(ptr[0][a]);
  } 
  // aを含む連結成分にxを作用させる
  void update_cc(int a, Lazy x){
    F.update_cc(ptr[0][a], x);
  }
  // aを含む連結成分のモノイド積　
  Val query_cc(int a){
    return F.query_cc(ptr[0][a]);
  }
};
#endif