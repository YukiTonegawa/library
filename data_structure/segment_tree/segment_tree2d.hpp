#ifndef _SEGMENT_TREE2D_H_
#define _SEGMENT_TREE2D_H_
#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <array>
#include <iostream>
#include <cassert>
#include <limits>
#include "../../algebraic_structure/monoid.hpp"
#include "segment_tree.hpp"

template<typename monoid>
struct segment_tree2d{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
  int N, N2, M, M2;
  std::vector<std::vector<Val>> sum;
  int ceil_pow2(int y){
    int x = 0;
    while ((1U << x) < (unsigned int)(y)) x++;
    return x;
  };
  segment_tree2d(): N(0), M(0){}
  segment_tree2d(int n, int m): N(n), N2(1 << ceil_pow2(N)), M(m), M2(1 << ceil_pow2(M)), sum(2 * N2 - 1, std::vector<Val>(2 * M2 - 1, id())){}
  segment_tree2d(const std::vector<std::vector<Val>> &v): N(v.size()), N2(1 << ceil_pow2(N)), sum(2 * N2 - 1){
    if(!N) return;
    M = v[0].size(), M2 = (1 << ceil_pow2(M));
    for(int i = 0; i < 2 * N2 - 1; i++) sum[i].resize(2 * M2 - 1);
    for(int i = 0; i < N; i++) for(int j = 0; j < M; j++) sum[N2 - 1 + i][M2 - 1 - j] = v[i][j];
    for(int i = N2 - 2; i >= 0; i--){
      for(int j = M2 - 1; j < 2 * M2 - 1; j++){
        sum[i][j] = merge(sum[i * 2 + 1][j], sum[i * 2 + 2][j]);
      }
    }
    for(int i = 0; i < 2 * N2 - 1; i++){
      for(int j = M2 - 2; j >= 0; j--){
        sum[i][j] = merge(sum[i][j * 2 + 1], sum[i][j * 2 + 2]);
      }
    }
  }
  int depth(){
    return N;
  }
  int width(){
    return M;
  }
  void set(int i, int j, Val x){
    assert(0 <= i && i < N);
    assert(0 <= j && j < M);
    i += N2 - 1, j += M2 - 1;
    while(true){
      sum[i][j] = (i >= N2 - 1 ? x : merge(sum[i * 2 + 1][j], sum[i * 2 + 2][j]));
      int tmp_j = j;
      while(tmp_j){
        tmp_j = (tmp_j - 1) / 2;
        sum[i][tmp_j] = merge(sum[i][tmp_j * 2 + 1], sum[i][tmp_j * 2 + 2]);
      }
      if(!i) return;
      i = (i - 1) / 2;
    }
  }
  Val get(int i, int j){
    assert(0 <= i && i < N);
    assert(0 <= j && j < M);
    return sum[N2 - 1 + i][M2 - 1 + j];
  }
  Val query(int lx, int rx, int ly, int ry){
    lx = std::max(lx, 0), rx = std::min(rx, N);
    ly = std::max(ly, 0), ry = std::min(ry, M);
    assert(lx <= rx && ly <= ry);
    lx += N2, rx += N2;
    ly += M2, ry += M2;
    Val ret = id();
    while(lx < rx){
      if(lx & 1){
        int tmp_ly = ly, tmp_ry = ry;
        while(tmp_ly < tmp_ry){
          if(tmp_ly & 1) ret = merge(ret, sum[lx - 1][(tmp_ly++) - 1]);
          if(tmp_ry & 1) ret = merge(ret, sum[lx - 1][(--tmp_ry) - 1]);
          tmp_ly >>= 1;
          tmp_ry >>= 1;
        }
        lx++;
      }
      if(rx & 1){
        rx--;
        int tmp_ly = ly, tmp_ry = ry;
        while(tmp_ly < tmp_ry){
          if(tmp_ly & 1) ret = merge(ret, sum[rx - 1][(tmp_ly++) - 1]);
          if(tmp_ry & 1) ret = merge(ret, sum[rx - 1][(--tmp_ry) - 1]);
          tmp_ly >>= 1;
          tmp_ry >>= 1;
        }
      }
      lx >>= 1;
      rx >>= 1;
    }
    return ret;
  }
  Val query_all(){
    return sum[0][0];
  }
};

template<typename Idx, typename monoid>
struct compressed_segment_tree2d{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge = monoid::merge;
private:
  using point = std::tuple<Idx, Idx, Val>;
  using segment_tree_1d = segment_tree<monoid>;
  std::vector<Idx> x, y;
  struct node{
    int lx, rx;
    node *l, *r;
    std::vector<std::array<int, 4>> next_idx;
    segment_tree_1d seg;
    node(int lx, int rx): lx(lx), rx(rx - 1), l(nullptr), r(nullptr){}
  };
  node *root;
  void build(node *v, int l, int r, const std::vector<int> &list,const std::vector<point> &p){
    std::vector<Val> value_list(list.size());
    for(int i = 0; i < list.size(); i++) value_list[i] = std::get<2>(p[list[i]]);
    v->seg = segment_tree_1d(value_list);
    if(r - l < 2) return;
    int mid = (l + r) / 2;
    Idx split_x = x[mid];
    v->l = new node(l, mid);
    v->r = new node(mid, r);
    int lsz = 0, rsz = 0;
    std::vector<int> left_list, right_list;
    for(int p_idx : list){
      if(std::get<0>(p[p_idx]) < split_x){
        v->next_idx.push_back(std::array<int, 4>{lsz + 1, lsz, rsz, rsz});
        left_list.push_back(p_idx);
        lsz++;
      }else{
        v->next_idx.push_back(std::array<int, 4>{lsz, lsz, rsz + 1, rsz});
        right_list.push_back(p_idx);
        rsz++;
      }
    }
    build(v->l, l, mid, left_list, p);
    build(v->r, mid, r, right_list, p);
  }
  void update(node *v, Idx ux, int uy, Val uz){
    Idx a = x[v->lx], b = x[v->rx];
    if(ux < a || b < ux) return;
    if(a <= ux && ux <= b) v->seg.update(uy, uz);
    if(v->l) update(v->l, ux, v->next_idx[uy][1], uz);
    if(v->r) update(v->r, ux, v->next_idx[uy][3], uz);
  }
  Val query(node *v, Idx lx, Idx rx, int ly, int ry){
    if(!v) return id();
    Idx a = x[v->lx], b = x[v->rx];
    if(rx <= a || b < lx || lx >= rx || ly >= ry) return id();
    if(lx <= a && b < rx) return v->seg.query(ly, ry);
    return merge(query(v->l, lx, rx, v->next_idx[ly][1], v->next_idx[ry - 1][0]),
             query(v->r, lx, rx, v->next_idx[ly][3], v->next_idx[ry - 1][2]));
  }
public:
  compressed_segment_tree2d(){}
  compressed_segment_tree2d(std::vector<point> p){
    std::sort(p.begin(), p.end(), [](point &A, point &B){return std::get<1>(A) < std::get<1>(B);});
    for(int i = 0; i < p.size(); i++){
      x.push_back(std::get<0>(p[i]));
      y.push_back(std::get<1>(p[i]));
    }
    std::sort(x.begin(), x.end());
    x.erase(std::unique(x.begin(), x.end()), x.end());
    std::vector<int> list(p.size());
    std::iota(list.begin(), list.end(), 0);
    root = new node(0, x.size());
    build(root, 0, x.size(), list, p);
  }
  void update(Idx ux, Idx uy, Val uz){
    int y_idx = std::lower_bound(y.begin(), y.end(), uy) - y.begin();
    update(root, ux, y_idx, uz);
  }
  Val query(Idx lx, Idx rx, Idx ly, Idx ry){
    int ly_idx = std::lower_bound(y.begin(), y.end(), ly) - y.begin();
    int ry_idx = std::lower_bound(y.begin(), y.end(), ry) - y.begin();
    return query(root, lx, rx, ly_idx, ry_idx);
  }
};
#endif