#ifndef _ARITHMETIC_ADD_RANGE_MIN_H_
#define _ARITHMETIC_ADD_RANGE_MIN_H_
#include <vector>
#include <limits>
#include <algorithm>
#include <cassert>

template<typename Val, typename Val2>
struct linear_add_range_min{
  static constexpr Val inf = std::numeric_limits<Val>::max() / 2;
  int N, M;
  int ceil_pow2(int y){
    int x = 0;
    while ((1U << x) < (unsigned int)(y)) x++;
    return x;
  };
  struct node{
    int L, R, bl, br;
    Val bly, bry;
    Val a, b;
    // 葉
    node(int i, Val y): bl(i), br(i), bly(y), bry(y), a(0), b(0){}
    // 葉以外(pullで初期化)
    node(): a(0), b(0){}
  };
  std::vector<node> nodes;
  void propagate(int k, Val a, Val b){
    assert(k < nodes.size());
    nodes[k].a += a;
    nodes[k].b += b;
    nodes[k].bly += a * nodes[k].bl + b;
    nodes[k].bry += a * nodes[k].br + b;
  }
  void push_down(int k){
    if(nodes[k].a == 0 && nodes[k].b == 0) return;
    propagate(k * 2 + 1, nodes[k].a, nodes[k].b);
    propagate(k * 2 + 2, nodes[k].a, nodes[k].b);
    nodes[k].a = nodes[k].b = 0;
  }
  void pull(int k){
    assert(k < M - 1);
    int l = k * 2 + 1, r = k * 2 + 2;
    int Mid = (nodes[k].L + nodes[k].R) / 2;
    while((l < M - 1) || (r < M - 1)){
      if(l < M - 1) push_down(l);
      if(r < M - 1) push_down(r);
      int a = nodes[l].bl, b = nodes[l].br, c = nodes[r].bl, d = nodes[r].br;
      Val ay = nodes[l].bly, by = nodes[l].bry, cy = nodes[r].bly, dy = nodes[r].bry;
      //Val2 s1 = cross(a, ay, b, by, c, cy);
      Val2 s1 = (Val2)(by - ay) * (c - a) - (Val2)(cy - ay) * (b - a);
      if(a != b && s1 > 0){
        l = l * 2 + 1;
      }else if(c != d && (Val2)(cy - by) * (d - b) - (Val2)(dy - by) * (c - b) > 0){
        // cross(b, by, c, cy, d, dy) > 0
        r = r * 2 + 2;
      }else if(a == b){
        r = r * 2 + 1;
      }else if(c == d){
        l = l * 2 + 2;
      }else{
        // Val2 s2 = cross(b, by, a, ay, d, dy);
        Val2 s2 = (Val2)(ay - by) * (d - b) - (Val2)(dy - by) * (a - b);
        //assert(s1 + s2 >= 0);
        if(s1 + s2 == 0 || s1 * (d - Mid) < s2 * (Mid - c)){
          l = l * 2 + 2;
        }else{
          r = r * 2 + 1;
        }
      }
    }
    nodes[k].bl = l - (M - 1);
    nodes[k].bly = nodes[l].bly;
    nodes[k].br = r - (M - 1);
    nodes[k].bry = nodes[r].bly;
  }
  void set(int k, int i, Val x){
    if(M - 1 <= k){
      nodes[k].bly = nodes[k].bry = x;
      return;
    }
    push_down(k);
    int mid = (nodes[k].L + nodes[k].R) / 2;
    if(i < mid) set(k * 2 + 1, i, x);
    else set(k * 2 + 2, i, x);
    pull(k);
  }
  void update(int k, int l, int r, Val a, Val b){
    if(nodes[k].R <= l || r <= nodes[k].L) return;
    if(l <= nodes[k].L && nodes[k].R <= r){
      propagate(k, a, b);
      return;
    }
    push_down(k);
    update(k * 2 + 1, l, r, a, b);
    update(k * 2 + 2, l, r, a, b);
    pull(k);
  }
  Val query(int k, int l, int r){
    if(nodes[k].R <= l || r <= nodes[k].L) return inf;
    if(l <= nodes[k].L && nodes[k].R <= r){
      while(k < M - 1){
        push_down(k);
        if(nodes[k].bly <= nodes[k].bry) k = k * 2 + 1;
        else k = k * 2 + 2;
      }
      return nodes[k].bly;
    }
    push_down(k);
    return std::min(query(k * 2 + 1, l, r), query(k * 2 + 2, l, r));
  }
  std::pair<Val, int> query_idx(int k, int l, int r){
    if(nodes[k].R <= l || r <= nodes[k].L) return {inf, -1};
    if(l <= nodes[k].L && nodes[k].R <= r){
      while(k < M - 1){
        push_down(k);
        if(nodes[k].bly <= nodes[k].bry) k = k * 2 + 1;
        else k = k * 2 + 2;
      }
      return {nodes[k].bly, nodes[k].bl};
    }
    push_down(k);
    return std::min(query_idx(k * 2 + 1, l, r), query_idx(k * 2 + 2, l, r));
  }
  Val query_assume_add(int k, int l, int r, Val a, Val b){
    if(nodes[k].R <= l || r <= nodes[k].L) return inf;
    if(l <= nodes[k].L && nodes[k].R <= r){
      while(k < M - 1){
        push_down(k);
        if(nodes[k].bly <= nodes[k].bry + a * (nodes[k].br - nodes[k].bl)) k = k * 2 + 1;
        else k = k * 2 + 2;
      }
      return nodes[k].bly + a * nodes[k].bl + b;
    }
    push_down(k);
    return std::min(query_assume_add(k * 2 + 1, l, r, a, b), query_assume_add(k * 2 + 2, l, r, a, b));
  }
  std::pair<int, Val> binary_search_left(int k, int l, Val a, Val b, Val x, Val &mn){
    if(nodes[k].R <= a) return {-1, mn};
    if(k < M - 1) push_down(k);
    if(l <= nodes[k].L){
      int v = k;
      while(v < M - 1){
        push_down(v);
        if(nodes[v].bly <= nodes[v].bry + a * (nodes[v].br - nodes[v].bl)) v = v * 2 + 1;
        else v = v * 2 + 2;
      }
      Val tmp = std::min(mn, nodes[v].bly + a * nodes[v].bl + b);
      if(tmp > x) return {-1, tmp};
      if(k >= M - 1) return {nodes[k].L, tmp};
    }
    auto left = bisect_from_left(k * 2 + 1, l, a, b, x, mn);
    if(left.first != -1) return left;
    mn = std::min(left.second, mn);
    return binary_search_left(k * 2 + 2, l, a, b, x, mn);
  }
  std::pair<int, Val> binary_search_right(int k, int r, Val a, Val b, Val x, Val &mn){
    if(r <= nodes[k].L) return {-1, mn};
    if(k < M - 1) push_down(k);
    if(nodes[k].R <= r){
      int v = k;
      while(v < M - 1){
        push_down(v);
        if(nodes[v].bly < nodes[v].bry + a * (nodes[v].br - nodes[v].bl)) v = v * 2 + 1;
        else v = v * 2 + 2;
      }
      Val tmp = std::min(mn, nodes[v].bly + a * nodes[v].bl + b);
      if(tmp > x) return {-1, tmp};
      if(k >= M - 1) return {nodes[k].L, tmp};
    }
    auto right = binary_search_right(k * 2 + 2, r, a, b, x, mn);
    if(right.first != -1) return right;
    mn = std::min(right.second, mn);
    return binary_search_right(k * 2 + 1, r, a, b, x, mn);
  }
public:
  linear_add_range_min(const std::vector<Val>& val): N(val.size()), M(1 << ceil_pow2(N)), nodes(2 * M - 1){
    for(int i = 0; i < M; i++){
      nodes[i + M - 1] = node(i, (i < N ? val[i] : inf));
      nodes[i + M - 1].L = i;
      nodes[i + M - 1].R = i + 1;
    }
    for(int i = M - 2; i >= 0; i--){
      nodes[i].L = nodes[i * 2 + 1].L;
      nodes[i].R = nodes[i * 2 + 2].R;
      pull(i);
    }
  }
  void set(int i, Val x){
    set(0, i, x);
  }
  Val get(int i){
    return query(i, i + 1);
  }
  // [l, r)にax+bを足す
  void update(int l, int r, Val a, Val b){
    update(0, l, r, a, b);
  }
  // [l, r)のmin
  Val query(int l, int r){
    return query(0, l, r);
  }
  // [l, r)の{min, そのインデックス}
  std::pair<Val, int> query_idx(int l, int r){
    return query_idx(0, l, r);
  }
  // [l, r)にax+bを足したと仮定してのmin
  Val query_assume_add(int l, int r, Val a, Val b){
    return query_assume_add(0, l, r, a, b);
  }
  /*
  // [i, r) (i + 1 < r)の下側凸包でiの次の点
  std::pair<int, Val> next(int i, int r){

  }
  // [l, i] (l < i)の下側凸包でiの前の点
  std::pair<int, Val> prev(int i, int l){

  }
  */
  // 全体にax+bを足したと仮定して値がx以下になるようなl以上の最小インデックス, ない場合は-1
  std::pair<int, Val> binary_search_left(int l, Val a, Val b, Val x){
    Val mn = inf;
    return binary_search_left(0, l, a, b, x, mn);
  }
  // 全体にax+bを足したと仮定して値がx以下になるようなr以下の最大インデックス, ない場合は-1
  int binary_search_right(int r, Val a, Val b, Val x){
    Val mn = inf;
    return binary_search_right(0, r, a, b, x, mn);
  }
};
#endif