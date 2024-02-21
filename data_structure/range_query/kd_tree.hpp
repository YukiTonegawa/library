#ifndef _KD_TREE_H_
#define _KD_TREE_H_
#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>
#include "../../algebraic_structure/monoid.hpp"

template<typename Idx, typename Val>
struct kd_tree{
  struct point{
    Idx x, y;
    Val z;
    point(){}
    point(Idx x, Idx y, Val z): x(x), y(y), z(z){}
    Idx access(bool dim)const{
      return !dim ? x : y;
    }
  };
  struct node{
    point p;
    int sz;
    Idx lx, rx, ly, ry;
    bool alive;
    node *l, *r, *u;
    node(Idx x, Idx y, Val z): p(x, y, z), sz(1), lx(x), rx(x), ly(y), ry(y), alive(true), l(nullptr), r(nullptr), u(nullptr){}
    node(point p): p(p), sz(1), lx(p.x), rx(p.x), ly(p.y), ry(p.y), alive(true), l(nullptr), r(nullptr), u(nullptr){}
  };
  static constexpr Idx inf = std::numeric_limits<Idx>::max() / 2;
  static constexpr Idx minf = std::numeric_limits<Idx>::min() / 2;
private:
  node *root;
  void update(node *v){
    v->sz = v->alive;
    if(v->alive){
      v->lx = v->rx = v->p.x;
      v->ly = v->ry = v->p.y;
    }else{
      v->lx = v->ly = inf;
      v->rx = v->ry = minf;
    }
    if(v->l){
      v->sz += v->l->sz;
      v->lx = std::min(v->lx, v->l->lx);
      v->rx = std::max(v->rx, v->l->rx);
      v->ly = std::min(v->ly, v->l->ly);
      v->ry = std::max(v->ry, v->l->ry);
    }
    if(v->r){
      v->sz += v->r->sz;
      v->lx = std::min(v->lx, v->r->lx);
      v->rx = std::max(v->rx, v->r->rx);
      v->ly = std::min(v->ly, v->r->ly);
      v->ry = std::max(v->ry, v->r->ry);
    }
  }
  node *build(std::vector<point> &v, int l, int r, bool dim){
    int m = (l + r) / 2;
    // current_dimで分割
    std::nth_element(v.begin() + l, v.begin() + m, v.begin() + r, [dim](const point &a, const point &b){
      return a.access(dim) < b.access(dim);
    });
    /*
    Idx border = v[m].access(dim);
    int m0 = m, l0 = l;
    while(m0 > l0){
      while(m0 > l0 && v[l0].access(dim) < border) l0++;
      while(m0 > l0 && v[m0].access(dim) == border) m0--;
      if(m0 != l0) std::swap(v[l0], v[m0]);
    }
    m = m0;
    */
    node *u = new node(v[m]);
    if(l < m){
      if((u->l = build(v, l, m, !dim))){
        u->l->u = u;
      }
    }
    if(m + 1 < r){
      if((u->r = build(v, m + 1, r, !dim))){
        u->r->u = u;
      }
    }
    update(u);
    return u;
  }
  int __range_freq(node *v, Idx lx, Idx rx, Idx ly, Idx ry){
    if(!v || !v->sz || v->rx < lx || v->lx >= rx || v->ry < ly || v->ly >= ry) return 0;
    if(lx <= v->lx && v->rx < rx && ly <= v->ly && v->ry < ry) return v->sz;
    int ret = 0;
    if(v->alive && lx <= v->p.x && v->p.x < rx && ly <= v->p.y && v->p.y < ry){
      ret = 1;
    }
    return ret + __range_freq(v->l, lx, rx, ly, ry) + __range_freq(v->r, lx, rx, ly, ry);
  }
  void __range_find(node *v, Idx lx, Idx rx, Idx ly, Idx ry, std::vector<node*> &ret){
    if(!v || !v->sz || v->rx < lx || v->lx >= rx || v->ry < ly || v->ly >= ry) return;
    if(v->alive && lx <= v->p.x && v->p.x < rx && ly <= v->p.y && v->p.y < ry){
      ret.push_back(v);
    }
    __range_find(v->l, lx, rx, ly, ry, ret);
    __range_find(v->r, lx, rx, ly, ry, ret);
  }
  void __range_find_manhattan(node *v, Idx x, Idx y, Idx dist, std::vector<node*> &ret){
    if(!v || !v->sz) return;
    Idx min_diff_x = (x < v->lx ? v->lx - x : x > v->rx ? x - v->rx : 0); // 部分木内に含まれる任意の点のx座標は少なくともこの距離離れている
    Idx min_diff_y = (y < v->ly ? v->ly - y : y > v->ry ? y - v->ry : 0); // 部分木内に含まれる任意の点のy座標は少なくともこの距離離れている
    if(min_diff_x + min_diff_y > dist) return;
    if(v->alive && abs(x - v->p.x) + abs(y - v->p.y) <= dist){
      ret.push_back(v);
    }
    __range_find_manhattan(v->l, x, y, dist, ret);
    __range_find_manhattan(v->r, x, y, dist, ret);
  }
  template<typename Dist>
  void __range_find_euclid(node *v, Idx x, Idx y, Dist dist2, std::vector<node*> &ret){
    if(!v || !v->sz) return;
    Dist min_diff_x = (x < v->lx ? v->lx - x : x > v->rx ? x - v->rx : 0); // 部分木内に含まれる任意の点のx座標は少なくともこの距離離れている
    Dist min_diff_y = (y < v->ly ? v->ly - y : y > v->ry ? y - v->ry : 0); // 部分木内に含まれる任意の点のy座標は少なくともこの距離離れている
    if(min_diff_x * min_diff_x + min_diff_y * min_diff_y > dist2) return;
    if(v->alive && (Dist)(x - v->p.x) * (x - v->p.x) + (Dist)(y - v->p.y) * (y - v->p.y) <= dist2){
      ret.push_back(v);
    }
    __range_find_euclid(v->l, x, y, dist2, ret);
    __range_find_euclid(v->r, x, y, dist2, ret);
  }
  void __erase(node *v){
    if(!v || !v->alive) return;
    v->alive = false;
    while(v){
      update(v);
      v = v->u;
    }
  }
public:
  kd_tree(): root(nullptr){}
  kd_tree(std::vector<point> v): root(nullptr){
    if(!v.empty()){
      root = build(v, 0, v.size(), 0);
    }
  }
  kd_tree(const std::vector<std::tuple<Idx, Idx, Val>> &v): root(nullptr){
    if(!v.empty()){
      std::vector<point> p;
      for(auto [x, y, z] : v) p.push_back(point(x, y, z));
      root = build(p, 0, p.size(), 0);
    }
  }
  int range_freq(Idx lx, Idx rx, Idx ly, Idx ry){
    return __range_freq(root, lx, rx, ly, ry);
  }
  // k := 返す点の数 として O(√N + k)
  std::vector<node*> range_find(Idx lx, Idx rx, Idx ly, Idx ry){
    std::vector<node*> ret;
    __range_find(root, lx, rx, ly, ry, ret);
    return ret;
  }
  // 点(x, y)からマンハッタン距離がdist以内の点
  std::vector<node*> range_find_manhattan(Idx x, Idx y, Idx dist){
    std::vector<node*> ret;
    __range_find_manhattan(root, x, y, dist, ret);
    return ret;
  }
  // 点(x, y)からユークリッド距離がdist以内の点
  template<typename Dist>
  std::vector<node*> range_find_euclid(Idx x, Idx y, Dist dist){
    std::vector<node*> ret;
    __range_find_euclid(root, x, y, dist * dist, ret);
    return ret;
  }
  // O(logN)
  void erase(node *v){
    __erase(v);
  }
};

template<typename Idx, typename monoid>
struct kd_tree_monoid{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;

  struct point{
    Idx x, y;
    Val z;
    point(){}
    point(Idx x, Idx y, Val z): x(x), y(y), z(z){}
    Idx access(bool dim)const{
      return !dim ? x : y;
    }
  };
  struct node{
    point p;
    int sz;
    Idx lx, rx, ly, ry;
    Val sum;
    bool alive;
    node *l, *r, *u;
    node(Idx x, Idx y, Val z): p(x, y, z), sz(1), lx(x), rx(x), ly(y), ry(y), sum(z), alive(true), l(nullptr), r(nullptr), u(nullptr){}
    node(point p): p(p), sz(1), lx(p.x), rx(p.x), ly(p.y), ry(p.y), sum(p.z), alive(true), l(nullptr), r(nullptr), u(nullptr){}
  };
  static constexpr Idx inf = std::numeric_limits<Idx>::max() / 2;
  static constexpr Idx minf = std::numeric_limits<Idx>::min() / 2;
private:
  node *root;
  void update(node *v){
    v->sz = v->alive;
    if(v->alive){
      v->lx = v->rx = v->p.x;
      v->ly = v->ry = v->p.y;
      v->sum = v->p.z;
    }else{
      v->lx = v->ly = inf;
      v->rx = v->ry = minf;
      v->sum = id();
    }
    if(v->l){
      v->sz += v->l->sz;
      v->lx = std::min(v->lx, v->l->lx);
      v->rx = std::max(v->rx, v->l->rx);
      v->ly = std::min(v->ly, v->l->ly);
      v->ry = std::max(v->ry, v->l->ry);
      v->sum = merge(v->l->sum, v->sum);
    }
    if(v->r){
      v->sz += v->r->sz;
      v->lx = std::min(v->lx, v->r->lx);
      v->rx = std::max(v->rx, v->r->rx);
      v->ly = std::min(v->ly, v->r->ly);
      v->ry = std::max(v->ry, v->r->ry);
      v->sum = merge(v->sum, v->r->sum);
    }
  }
  node *build(std::vector<point> &v, int l, int r, bool dim){
    int m = (l + r) / 2;
    // current_dimで分割
    std::nth_element(v.begin() + l, v.begin() + m, v.begin() + r, [dim](const point &a, const point &b){
      return a.access(dim) < b.access(dim);
    });
    /*
    Idx border = v[m].access(dim);
    int m0 = m, l0 = l;
    while(m0 > l0){
      while(m0 > l0 && v[l0].access(dim) < border) l0++;
      while(m0 > l0 && v[m0].access(dim) == border) m0--;
      if(m0 != l0) std::swap(v[l0], v[m0]);
    }
    m = m0;
    */
    node *u = new node(v[m]);
    if(l < m){
      if((u->l = build(v, l, m, !dim))){
        u->l->u = u;
      }
    }
    if(m + 1 < r){
      if((u->r = build(v, m + 1, r, !dim))){
        u->r->u = u;
      }
    }
    update(u);
    return u;
  }
  int __range_freq(node *v, Idx lx, Idx rx, Idx ly, Idx ry){
    if(!v || !v->sz || v->rx < lx || v->lx >= rx || v->ry < ly || v->ly >= ry) return 0;
    if(lx <= v->lx && v->rx < rx && ly <= v->ly && v->ry < ry) return v->sz;
    int ret = 0;
    if(v->alive && lx <= v->p.x && v->p.x < rx && ly <= v->p.y && v->p.y < ry){
      ret = 1;
    }
    return ret + __range_freq(v->l, lx, rx, ly, ry) + __range_freq(v->r, lx, rx, ly, ry);
  }
  Val __query(node *v, Idx lx, Idx rx, Idx ly, Idx ry){
    if(!v || !v->sz || v->rx < lx || v->lx >= rx || v->ry < ly || v->ly >= ry) return id();
    if(lx <= v->lx && v->rx < rx && ly <= v->ly && v->ry < ry) return v->sum;
    Val ret = id();
    if(v->alive && lx <= v->p.x && v->p.x < rx && ly <= v->p.y && v->p.y < ry){
      ret = v->p.z;
    }
    return merge(ret, merge(__query(v->l, lx, rx, ly, ry) + __query(v->r, lx, rx, ly, ry)));
  }
  void __range_find(node *v, Idx lx, Idx rx, Idx ly, Idx ry, std::vector<node*> &ret){
    if(!v || !v->sz || v->rx < lx || v->lx >= rx || v->ry < ly || v->ly >= ry) return;
    if(v->alive && lx <= v->p.x && v->p.x < rx && ly <= v->p.y && v->p.y < ry){
      ret.push_back(v);
    }
    __range_find(v->l, lx, rx, ly, ry, ret);
    __range_find(v->r, lx, rx, ly, ry, ret);
  }
  void __range_find_manhattan(node *v, Idx x, Idx y, Idx dist, std::vector<node*> &ret){
    if(!v || !v->sz) return;
    Idx min_diff_x = (x < v->lx ? v->lx - x : x > v->rx ? x - v->rx : 0); // 部分木内に含まれる任意の点のx座標は少なくともこの距離離れている
    Idx min_diff_y = (y < v->ly ? v->ly - y : y > v->ry ? y - v->ry : 0); // 部分木内に含まれる任意の点のy座標は少なくともこの距離離れている
    if(min_diff_x + min_diff_y > dist) return;
    if(v->alive && abs(x - v->p.x) + abs(y - v->p.y) <= dist){
      ret.push_back(v);
    }
    __range_find_manhattan(v->l, x, y, dist, ret);
    __range_find_manhattan(v->r, x, y, dist, ret);
  }
  Val __query_manhattan(node *v, Idx x, Idx y, Idx dist){
    if(!v || !v->sz) return id();
    Idx min_diff_x = (x < v->lx ? v->lx - x : x > v->rx ? x - v->rx : 0); // 部分木内に含まれる任意の点のx座標は少なくともこの距離離れている
    Idx min_diff_y = (y < v->ly ? v->ly - y : y > v->ry ? y - v->ry : 0); // 部分木内に含まれる任意の点のy座標は少なくともこの距離離れている
    if(min_diff_x + min_diff_y > dist) return id();
    Idx max_diff_x = std::max(abs(x - v->lx), abs(x - v->rx));
    Idx max_diff_y = std::max(abs(y - v->ly), abs(y - v->ry));
    if(max_diff_x + max_diff_y <= dist) return v->sum;
    Val ret = id();
    if(v->alive && abs(x - v->p.x) + abs(y - v->p.y) <= dist){
      ret = v->p.z;
    }
    return merge(ret, merge(__query_manhattan(v->l, x, y, dist), __query_manhattan(v->r, x, y, dist)));
  }
  template<typename Dist>
  void __range_find_euclid(node *v, Idx x, Idx y, Dist dist2, std::vector<node*> &ret){
    if(!v || !v->sz) return;
    Dist min_diff_x = (x < v->lx ? v->lx - x : x > v->rx ? x - v->rx : 0); // 部分木内に含まれる任意の点のx座標は少なくともこの距離離れている
    Dist min_diff_y = (y < v->ly ? v->ly - y : y > v->ry ? y - v->ry : 0); // 部分木内に含まれる任意の点のy座標は少なくともこの距離離れている
    if(min_diff_x * min_diff_x + min_diff_y * min_diff_y > dist2) return;
    if(v->alive && (Dist)(x - v->p.x) * (x - v->p.x) + (Dist)(y - v->p.y) * (y - v->p.y) <= dist2){
      ret.push_back(v);
    }
    __range_find_euclid(v->l, x, y, dist2, ret);
    __range_find_euclid(v->r, x, y, dist2, ret);
  }
  template<typename Dist>
  Val __query_euclid(node *v, Idx x, Idx y, Dist dist2){
    if(!v || !v->sz) return id();
    Dist min_diff_x = (x < v->lx ? v->lx - x : x > v->rx ? x - v->rx : 0); // 部分木内に含まれる任意の点のx座標は少なくともこの距離離れている
    Dist min_diff_y = (y < v->ly ? v->ly - y : y > v->ry ? y - v->ry : 0); // 部分木内に含まれる任意の点のy座標は少なくともこの距離離れている
    if(min_diff_x * min_diff_x + min_diff_y * min_diff_y > dist2) return id();
    Dist max_diff_x = std::max(abs(x - v->lx), abs(x - v->rx));
    Dist max_diff_y = std::max(abs(y - v->ly), abs(y - v->ry));
    if(max_diff_x * max_diff_x + max_diff_y + max_diff_y <= dist2) return v->sum;
    Val ret = id();
    if(v->alive && (Dist)(x - v->p.x) * (x - v->p.x) + (Dist)(y - v->p.y) * (y - v->p.y) <= dist2){
      ret = v->p.z;
    }
    return merge(ret, merge(__query_euclid<Dist>(v->l, x, y, dist2), __query_euclid<Dist>(v->r, x, y, dist2)));
  }
  void __erase(node *v){
    if(!v || !v->alive) return;
    v->alive = false;
    while(v){
      update(v);
      v = v->u;
    }
  }
  void __set_val(node *v, Val z){
    if(!v || !v->alive) return;
    v->p->z = z;
    while(v){
      update(v);
      v = v->u;
    }
  }
public:
  kd_tree_monoid(): root(nullptr){}
  kd_tree_monoid(std::vector<point> v): root(nullptr){
    if(!v.empty()){
      root = build(v, 0, v.size(), 0);
    }
  }
  kd_tree_monoid(const std::vector<std::tuple<Idx, Idx, Val>> &v): root(nullptr){
    if(!v.empty()){
      std::vector<point> p;
      for(auto [x, y, z] : v) p.push_back(point(x, y, z));
      root = build(p, 0, p.size(), 0);
    }
  }
  int range_freq(Idx lx, Idx rx, Idx ly, Idx ry){
    return __range_freq(root, lx, rx, ly, ry);
  }
  // O(√N)
  Val query(Idx lx, Idx rx, Idx ly, Idx ry){
    return __query(root, lx, rx, ly, ry);
  }
  // k := 返す点の数 として O(√N + k)
  std::vector<node*> range_find(Idx lx, Idx rx, Idx ly, Idx ry){
    std::vector<node*> ret;
    __range_find(root, lx, rx, ly, ry, ret);
    return ret;
  }
  // 点(x, y)からマンハッタン距離がdist以内の点
  std::vector<node*> range_find_manhattan(Idx x, Idx y, Idx dist){
    std::vector<node*> ret;
    __range_find_manhattan(root, x, y, dist, ret);
    return ret;
  }
  Val query_manhattan(Idx x, Idx y, Idx dist){
    return __query_manhattan(root, x, y, dist);
  }
  // 点(x, y)からユークリッド距離がdist以内の点
  template<typename Dist>
  std::vector<node*> range_find_euclid(Idx x, Idx y, Dist dist){
    std::vector<node*> ret;
    __range_find_euclid(root, x, y, dist * dist, ret);
    return ret;
  }
  template<typename Dist>
  Val query_euclid(Idx x, Idx y, Dist dist){
    return __query_euclid<Dist>(x, y, dist);
  }
  // O(logN)
  void erase(node *v){
    __erase(v);
  }
  // O(logN)
  void set_val(node *v, Val z){
    __set_val(v, z);
  }
};

template<typename Idx, typename Val>
struct kd_tree_3d{
  struct point{
    Idx x, y, z;
    Val w;
    point(){}
    point(Idx x, Idx y, Idx z, Val w): x(x), y(y), z(z), w(w){}
    Idx access(int dim)const{
      return dim == 0 ? x : dim == 1 ? y : z;
    }
  };
  struct node{
    point p;
    int sz;
    Idx lx, rx, ly, ry, lz, rz;
    bool alive;
    node *l, *r, *u;
    node(Idx x, Idx y, Idx z, Val w): p(x, y, z, w), sz(1), lx(x), rx(x), ly(y), ry(y), lz(z), rz(z), alive(true), l(nullptr), r(nullptr), u(nullptr){}
    node(point p): p(p), sz(1), lx(p.x), rx(p.x), ly(p.y), ry(p.y), lz(p.z), rz(p.z), alive(true), l(nullptr), r(nullptr), u(nullptr){}
  };
  static constexpr Idx inf = std::numeric_limits<Idx>::max() / 2;
  static constexpr Idx minf = std::numeric_limits<Idx>::min() / 2;
private:
  node *root;
  void update(node *v){
    v->sz = v->alive;
    if(v->alive){
      v->lx = v->rx = v->p.x;
      v->ly = v->ry = v->p.y;
      v->lz = v->rz = v->p.z;
    }else{
      v->lx = v->ly = v->lz = inf;
      v->rx = v->ry = v->rz = minf;
    }
    if(v->l){
      v->sz += v->l->sz;
      v->lx = std::min(v->lx, v->l->lx);
      v->rx = std::max(v->rx, v->l->rx);
      v->ly = std::min(v->ly, v->l->ly);
      v->ry = std::max(v->ry, v->l->ry);
      v->lz = std::min(v->lz, v->l->lz);
      v->rz = std::max(v->rz, v->l->rz);
    }
    if(v->r){
      v->sz += v->r->sz;
      v->lx = std::min(v->lx, v->r->lx);
      v->rx = std::max(v->rx, v->r->rx);
      v->ly = std::min(v->ly, v->r->ly);
      v->ry = std::max(v->ry, v->r->ry);
      v->lz = std::min(v->lz, v->r->lz);
      v->rz = std::max(v->rz, v->r->rz);
    }
  }
  node *build(std::vector<point> &v, int l, int r, int dim){
    int m = (l + r) / 2;
    // current_dimで分割
    std::nth_element(v.begin() + l, v.begin() + m, v.begin() + r, [dim](const point &a, const point &b){
      return a.access(dim) < b.access(dim);
    });
    /*
    Idx border = v[m].access(dim);
    int m0 = m, l0 = l;
    while(m0 > l0){
      while(m0 > l0 && v[l0].access(dim) < border) l0++;
      while(m0 > l0 && v[m0].access(dim) == border) m0--;
      if(m0 != l0) std::swap(v[l0], v[m0]);
    }
    m = m0;
    */
    node *u = new node(v[m]);
    if(l < m){
      if((u->l = build(v, l, m, (dim + 1) % 3))){
        u->l->u = u;
      }
    }
    if(m + 1 < r){
      if((u->r = build(v, m + 1, r, (dim + 1) % 3))){
        u->r->u = u;
      }
    }
    update(u);
    return u;
  }
  int __range_freq(node *v, Idx lx, Idx rx, Idx ly, Idx ry, Idx lz, Idx rz){
    if(!v || !v->sz || v->rx < lx || v->lx >= rx || v->ry < ly || v->ly >= ry || v->rz < lz || v->lz >= rz) return 0;
    if(lx <= v->lx && v->rx < rx && ly <= v->ly && v->ry < ry && lz <= v->lz && v->rz < rz) return v->sz;
    int ret = 0;
    if(v->alive && lx <= v->p.x && v->p.x < rx && ly <= v->p.y && v->p.y < ry && lz <= v->p.z && v->p.z < rz){
      ret = 1;
    }
    return ret + __range_freq(v->l, lx, rx, ly, ry, lz, rz) + __range_freq(v->r, lx, rx, ly, ry, lz, rz);
  }
  void __range_find(node *v, Idx lx, Idx rx, Idx ly, Idx ry, Idx lz, Idx rz, std::vector<node*> &ret){
    if(!v || !v->sz || v->rx < lx || v->lx >= rx || v->ry < ly || v->ly >= ry || v->rz < lz || v->lz >= rz) return;
    if(v->alive && lx <= v->p.x && v->p.x < rx && ly <= v->p.y && v->p.y < ry && lz <= v->p.z && v->p.z < rz){
      ret.push_back(v);
    }
    __range_find(v->l, lx, rx, ly, ry, lz, rz, ret);
    __range_find(v->r, lx, rx, ly, ry, lz, rz, ret);
  }
  void __range_find_manhattan(node *v, Idx x, Idx y, Idx z, Idx dist, std::vector<node*> &ret){
    if(!v || !v->sz) return;
    Idx min_diff_x = (x < v->lx ? v->lx - x : x > v->rx ? x - v->rx : 0); // 部分木内に含まれる任意の点のx座標は少なくともこの距離離れている
    Idx min_diff_y = (y < v->ly ? v->ly - y : y > v->ry ? y - v->ry : 0); // 部分木内に含まれる任意の点のy座標は少なくともこの距離離れている
    Idx min_diff_z = (z < v->lz ? v->lz - z : z > v->rz ? z - v->rz : 0); // 部分木内に含まれる任意の点のz座標は少なくともこの距離離れている
    if(min_diff_x + min_diff_y + min_diff_z > dist) return;
    if(v->alive && abs(x - v->p.x) + abs(y - v->p.y) + abs(z - v->p.z) <= dist){
      ret.push_back(v);
    }
    __range_find_manhattan(v->l, x, y, z, dist, ret);
    __range_find_manhattan(v->r, x, y, z, dist, ret);
  }
  template<typename Dist>
  void __range_find_euclid(node *v, Idx x, Idx y, Idx z, Dist dist2, std::vector<node*> &ret){
    if(!v || !v->sz) return;
    Dist min_diff_x = (x < v->lx ? v->lx - x : x > v->rx ? x - v->rx : 0); // 部分木内に含まれる任意の点のx座標は少なくともこの距離離れている
    Dist min_diff_y = (y < v->ly ? v->ly - y : y > v->ry ? y - v->ry : 0); // 部分木内に含まれる任意の点のy座標は少なくともこの距離離れている
    Dist min_diff_z = (z < v->lz ? v->lz - z : z > v->rz ? z - v->rz : 0); // 部分木内に含まれる任意の点のz座標は少なくともこの距離離れている
    if(min_diff_x * min_diff_x + min_diff_y * min_diff_y + min_diff_z * min_diff_z > dist2) return;
    if(v->alive && (Dist)(x - v->p.x) * (x - v->p.x) + (Dist)(y - v->p.y) * (y - v->p.y) + (Dist)(z - v->p.z) * (z - v->p.z) <= dist2){
      ret.push_back(v);
    }
    __range_find_euclid(v->l, x, y, z, dist2, ret);
    __range_find_euclid(v->r, x, y, z, dist2, ret);
  }
  void __erase(node *v){
    if(!v || !v->alive) return;
    v->alive = false;
    while(v){
      update(v);
      v = v->u;
    }
  }
public:
  kd_tree_3d(): root(nullptr){}
  kd_tree_3d(std::vector<point> v): root(nullptr){
    if(!v.empty()){
      root = build(v, 0, v.size(), 0);
    }
  }
  kd_tree_3d(const std::vector<std::tuple<Idx, Idx, Idx, Val>> &v): root(nullptr){
    if(!v.empty()){
      std::vector<point> p;
      for(auto [x, y, z, w] : v) p.push_back(point(x, y, z, w));
      root = build(p, 0, p.size(), 0);
    }
  }
  int range_freq(Idx lx, Idx rx, Idx ly, Idx ry, Idx lz, Idx rz){
    return __range_freq(root, lx, rx, ly, ry, lz, rz);
  }
  std::vector<node*> range_find(Idx lx, Idx rx, Idx ly, Idx ry, Idx lz, Idx rz){
    std::vector<node*> ret;
    __range_find(root, lx, rx, ly, ry, lz, rz, ret);
    return ret;
  }
  // 点(x, y)からマンハッタン距離がdist以内の点
  std::vector<node*> range_find_manhattan(Idx x, Idx y, Idx z, Idx dist){
    std::vector<node*> ret;
    __range_find_manhattan(root, x, y, z, dist, ret);
    return ret;
  }
  // 点(x, y)からユークリッド距離がdist以内の点
  template<typename Dist>
  std::vector<node*> range_find_euclid(Idx x, Idx y, Idx z, Dist dist){
    std::vector<node*> ret;
    __range_find_euclid(root, x, y, z, dist * dist, ret);
    return ret;
  }
  // O(logN)
  void erase(node *v){
    __erase(v);
  }
};
template<typename Idx, typename Val>
struct kd_tree_rectangle{
  struct point{
    Idx lx, rx, ly, ry;
    Val z;
    point(){}
    point(Idx __lx, Idx __rx, Idx __ly, Idx __ry, Val __z): lx(__lx), rx(__rx - 1), ly(__ly), ry(__ry - 1), z(__z){assert(lx <= rx && ly <= ry);}
    Idx access(int dim)const{
      return dim < 2 ? (dim == 0 ? lx : rx) : (dim == 2 ? ly : ry);
    }
  };
  struct node{
    point p;
    int sz;
    Idx llx, rrx, lly, rry; // lxの最小値, rxの最大値, lyの最小値, ryの最大値
    bool alive;
    node *l, *r, *u;
    node(Idx lx, Idx rx, Idx ly, Idx ry, Val z): p(lx, rx, ly, ry, z), sz(1), llx(lx), rrx(rx), lly(ly), rry(ry), alive(true), l(nullptr), r(nullptr), u(nullptr){}
    node(point p): p(p), sz(1), llx(p.lx), rrx(p.rx), lly(p.ly), rry(p.ry), alive(true), l(nullptr), r(nullptr), u(nullptr){}
  };
  static constexpr Idx inf = std::numeric_limits<Idx>::max() / 2;
  static constexpr Idx minf = std::numeric_limits<Idx>::min() / 2;
private:
  node *root;
  void update(node *v){
    v->sz = v->alive;
    if(v->alive){
      v->llx = v->p.lx;
      v->rrx = v->p.rx;
      v->lly = v->p.ly;
      v->rry = v->p.ry;
    }else{
      v->llx = v->lly = inf;
      v->rrx = v->rry = minf;
    }
    if(v->l){
      v->sz += v->l->sz;
      v->llx = std::min(v->llx, v->l->llx);
      v->rrx = std::max(v->rrx, v->l->rrx);
      v->lly = std::min(v->lly, v->l->lly);
      v->rry = std::max(v->rry, v->l->rry);
    }
    if(v->r){
      v->sz += v->r->sz;
      v->llx = std::min(v->llx, v->r->llx);
      v->rrx = std::max(v->rrx, v->r->rrx);
      v->lly = std::min(v->lly, v->r->lly);
      v->rry = std::max(v->rry, v->r->rry);
    }
  }
  node *build(std::vector<point> &v, int l, int r, int dim){
    int m = (l + r) / 2;
    // current_dimで分割
    std::nth_element(v.begin() + l, v.begin() + m, v.begin() + r, [dim](const point &a, const point &b){
      return a.access(dim) < b.access(dim);
    });
    /*
    Idx border = v[m].access(dim);
    int m0 = m, l0 = l;
    while(m0 > l0){
      while(m0 > l0 && v[l0].access(dim) < border) l0++;
      while(m0 > l0 && v[m0].access(dim) == border) m0--;
      if(m0 != l0) std::swap(v[l0], v[m0]);
    }
    m = m0;
    */
    node *u = new node(v[m]);
    if(l < m && (u->l = build(v, l, m, (dim + 1) & 3))){
      u->l->u = u;
    }
    if(m + 1 < r && (u->r = build(v, m + 1, r, (dim + 1) & 3))){
      u->r->u = u;
    }
    update(u);
    return u;
  }
  void __point(node *v, Idx x, Idx y, std::vector<node*> &res){
    if(!v || !v->sz || v->rrx < x || v->llx > x || v->rry < y || v->lly > y) return;
    if(v->alive && v->p.lx <= x && x <= v->p.rx && v->p.ly <= y && y <= v->p.ry){
      res.push_back(v);
    }
    __point(v->l, x, y, res);
    __point(v->r, x, y, res);
  }
  /*
  void __rectangle(node *v, Idx lx, Idx rx, Idx ly, Idx ry, std::vector<node*> &ret){
    if(!v || !v->sz || v->rx < lx || v->lx >= rx || v->ry < ly || v->ly >= ry || v->rz < lz || v->lz >= rz) return;
    if(v->alive && lx <= v->p.x && v->p.x < rx && ly <= v->p.y && v->p.y < ry && lz <= v->p.z && v->p.z < rz){
      ret.push_back(v);
    }
    __rectangle(v->l, lx, rx, ly, ry, lz, rz, ret);
    __rectangle(v->r, lx, rx, ly, ry, lz, rz, ret);
  }
  */
  void __erase(node *v){
    if(!v || !v->alive) return;
    v->alive = false;
    while(v){
      update(v);
      v = v->u;
    }
  }
public:
  kd_tree_rectangle(): root(nullptr){}
  kd_tree_rectangle(const std::vector<std::tuple<Idx, Idx, Idx, Idx, Val>> &v): root(nullptr){
    if(!v.empty()){
      std::vector<point> p;
      for(auto [lx, rx, ly, ry, z] : v) p.push_back(point(lx, rx, ly, ry, z));
      root = build(p, 0, p.size(), 0);
    }
  }
  std::vector<node*> fp(Idx x, Idx y){
    std::vector<node*> res;
    __point(root, x, y, res);
    return res;
  }
  /*
  std::vector<node*> range_find(Idx lx, Idx rx, Idx ly, Idx ry, Idx lz, Idx rz){
    std::vector<node*> ret;
    __range_find(root, lx, rx, ly, ry, lz, rz, ret);
    return ret;
  }
  */
  void erase(node *v){
    __erase(v);
  }
};
#endif