#ifndef _DYNAMIC_KD_TREE_H_
#define _DYNAMIC_KD_TREE_H_
#include "../bit_sequence/bit_operation.hpp"
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

// insert(x, y)
// erase(x, y)
// find_rect(lx, rx, ly, ry)
// range_freq(lx, rx, ly, ry)

// 0 < Dim <= 31
template<typename Key, int Dim>
struct kd_tree{
  using point = std::array<Key, Dim>;
private:
  static constexpr float alpha = 0.6;
  struct node{
    point p, lp, rp;
    uint32_t dim_flag;
    int current_dim;
    bool is_alive;
    int size;
    node *l, *m, *r;
    node(){}
    node(const point &p): p(p), lp(p), rp(p), dim_flag(((uint32_t)1 << Dim) - 1), is_alive(true), size(1), l(nullptr), m(nullptr), r(nullptr){}
  };
  node *make_node(const point &p){return new node(p);}
  int next_dim(uint32_t dim_flag, int current_dim){
    if(!dim_flag) return -1;
    int res = find_next_32bit(dim_flag, current_dim + 1);
    if(res == -1) res = __builtin_ctz(dim_flag);
    return res;
  }
  int max_node_count;
  node *root;
  inline int size(node *v){return !v ? 0 : v->size;}
  inline bool check_alpha_weight_balanced(node *v){
    int threshold = alpha * v->size;
    return size(v->l) <= threshold && size(v->r) <= threshold;
  }
  void update(node *v){
    v->size = v->is_alive;
    for(int i = 0; i < Dim; i++) v->lp[i] = v->rp[i] = v->p[i];
    if(v->l){
      v->size += v->l->size;
      for(int i = 0; i < Dim; i++){
        v->lp[i] = std::min(v->lp[i], v->l->lp[i]);
        v->rp[i] = std::max(v->rp[i], v->l->rp[i]);
      }
    }
    if(v->m){
      v->size += v->m->size;
      for(int i = 0; i < Dim; i++){
        v->lp[i] = std::min(v->lp[i], v->m->lp[i]);
        v->rp[i] = std::max(v->rp[i], v->m->rp[i]);
      }
    }
    if(v->r){
      v->size += v->r->size;
      for(int i = 0; i < Dim; i++){
        v->lp[i] = std::min(v->lp[i], v->r->lp[i]);
        v->rp[i] = std::max(v->rp[i], v->r->rp[i]);
      }
    }
  }
  void subtree_enumerate(node *v, std::vector<node*> &nodes, int &cnt){
    if(v->l && v->l->size) subtree_enumerate(v->l, nodes, cnt);
    if(v->m && v->m->size) subtree_enumerate(v->m, nodes, cnt);
    if(v->is_alive) nodes[cnt++] = v;
    if(v->r && v->r->size) subtree_enumerate(v->r, nodes, cnt);
  }
  node* rebalance(node *v){
    if(!v || !v->size) return nullptr;
    std::vector<node*> ret(v->size);
    int cnt = 0;
    subtree_enumerate(v, ret, cnt);
    return build(0, v->size, ret, v->dim_flag, v->current_dim);
  }
  node *insert(node *v, point &p){
    if(!v) return make_node(p);
    if(v->p == p) return v;
    else if(v->p[v->current_dim] > p[v->current_dim]) v->l = insert(v->l, p);
    else if(v->p[v->current_dim] == p[v->current_dim]) v->m = insert(v->m, p);
    else v->r = insert(v->r, p);
    if(!check_alpha_weight_balanced(v)) return rebalance(v);
    update(v);
    return v;
  }
  void erase(node *v, point &p){
    if(!v) return;
    if(v->p == p) v->is_alive = false;
    else if(v->p[v->current_dim] > p[v->current_dim]) v->l = erase(v->l, p);
    else if(v->p[v->current_dim] == p[v->current_dim]) v->m = erase(v->m, p);
    else v->r = erase(v->r, p);
    update(v);
  }
  node *build(int l, int r, std::vector<node*> &nodes, uint32_t dim_flag, int current_dim){
    int m = (l + r) / 2;
    // current_dimで分割
    std::nth_element(nodes.begin() + l, nodes.begin() + m,  nodes.begin() + r, [current_dim](node *a, node *b){
      return a->p[current_dim] < b->p[current_dim];
    });
    Key border = nodes[m]->p[current_dim];
    // [中央値より小さい], [中央値], [中央値より大きい]
    //     [l, m0),    [m0, m1),    [m1, r)
    int m0 = m, m1 = m;
    int l0 = l, r0 = r - 1;
    while(m0 > l0){
      while(m0 > l0 && nodes[l0]->p[current_dim] < border) l0++;
      while(m0 > l0 && nodes[m0]->p[current_dim] == border) m0--;
      if(m0 != l0) std::swap(nodes[l0], nodes[m0]);
    }
    while(r0 > m1){
      while(r0 > m1 && nodes[r0]->p[current_dim] > border) r0--;
      while(r0 > m1 && nodes[m1]->p[current_dim] == border) m1++;
      if(m1 != r0) std::swap(nodes[r0], nodes[m1]);
    }
    m1++;
    node *v = nodes[m0];
    v->current_dim = current_dim;
    v->dim_flag = dim_flag;
    int d1 = next_dim(dim_flag, current_dim);
    if(l < m0) v->l = build(l, m0, nodes, dim_flag, d1);
    else v->l = nullptr;
    if(m0 + 1 < m1){
      int dim_flag_mid = dim_flag ^ (1 << current_dim);
      v->m = build(m0 + 1, m1, nodes, dim_flag_mid, next_dim(dim_flag_mid, current_dim));
    }else v->m = nullptr;
    if(m1 < r) v->r = build(m1, r, nodes, dim_flag, d1);
    else v->r = nullptr;
    update(v);
    return v;
  }
  // l, rが点vpを含むか
  bool contain_point(const point &vp, const point &l, const point &r){
    for(int i = 0; i < Dim; i++) if(vp[i] < l[i] || r[i] <= vp[i]) return false;
    return true;
  }
  // l, rが領域vlp, vrpを含むか
  // 0 : 完全に含む
  // 1 : 全く重ならない
  // 2 : 部分的に重なる
  int contain_rect(const point &vlp, const point &vrp, const point &l, const point &r){
    bool f = true, g = true;
    for(int i = 0; i < Dim; i++){
      if(f && (vlp[i] < l[i] || r[i] <= vrp[i])) f = false;
      if(g && (vrp[i] < l[i] || r[i] <= vlp[i])) g = false;
    }
    if(f) return 0;
    return !g ? 1 : 2;
  }
  int range_freq(const node *v, const point &l, const point &r){
    if(!v || !v->size) return 0;
    int f = contain_rect(v->lp, v->rp, l, r);
    if(f == 1) return 0;
    if(f == 0) return v->size;
    int res = (v->is_alive && contain_point(v->p, l, r));
    return res + range_freq(v->l, l, r) + range_freq(v->m, l, r) + range_freq(v->r, l, r);
  }
  // debug
  int max_depth_finder(node *v, int d){
    if(!v || !v->is_alive) return d;
    return std::max({max_depth_finder(v->l, d + 1), max_depth_finder(v->m, d + 1), max_depth_finder(v->r, d + 1)});
  }
  int max_dep(){
    return max_depth_finder(root, 0);
  }
public:
  kd_tree(): max_node_count(0), root(nullptr){}
  kd_tree(const std::vector<point> &v){
    int n = v.size();
    if(!n){
      root = nullptr;
      return;
    }
    std::vector<node*> nodes(n);
    for(int i = 0; i < n; i++) nodes[i] = make_node(v[i]);
    root = build(0, n, nodes, ((uint32_t)1 << Dim) - 1, 0);
  }
  int size(){
    return size(root);
  }
  void insert(point p){
    root = insert(root, p);
    max_node_count = std::max(max_node_count, size(root));
  }
  void erase(point p){
    erase(root, p);
    if(size(root) <= alpha * max_node_count) root = rebalance(root);
  }
  int range_freq(point l, point r){
    return range_freq(root, l, r);
  }
};
#endif