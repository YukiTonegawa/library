#ifndef _BEATS_CHMIN_ADD_MAXSUM_H_
#define _BEATS_CHMIN_ADD_MAXSUM_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include <limits>
// 以下の操作が行える
//  xi += k in [l, r)
//  yi += k in [l, r)
//  chmin(xi, k) in [l, r)
//  chmin(yi, k) in [l, r)
//  get max(xi + yi) in [l, r)
// ならしO(log^2N)
template<typename Val>
struct beats_chmin_add_maxsum{
  static constexpr Val inf = std::numeric_limits<Val>::max();
  static constexpr Val minf = std::numeric_limits<Val>::min();
private:
  int N, M;
  int ceil_pow2(int y){
    int x = 0;
    while ((1U << x) < (unsigned int)(y)) x++;
    return x;
  };
  struct T{
    Val maxx, maxy;
    Val smaxx, smaxy;
    std::array<Val, 4> maxz;
    Val addx, addy, upperx, uppery;
    T(): maxx(minf), maxy(minf), smaxx(minf), smaxy(minf), addx(0), addy(0), upperx(inf), uppery(inf){maxz.fill(minf);}
  };
  std::vector<T> val;
  void propagate(int k, int len, const T &lz){
    Val add = 0;
    if(lz.addx){
      add += lz.addx;
      val[k].addx += lz.addx;
      if(val[k].upperx != inf) val[k].upperx += lz.addx;
      val[k].maxx += lz.addx;
      if(val[k].smaxx != minf) val[k].smaxx += lz.addx;
    }
    if(lz.addy){
      add += lz.addy;
      val[k].addy += lz.addy;
      if(val[k].uppery != inf) val[k].uppery += lz.addy;
      val[k].maxy += lz.addy;
      if(val[k].smaxy != minf) val[k].smaxy += lz.addy;
    }
    if(val[k].maxz[0] != minf) val[k].maxz[0] += add;
    if(val[k].maxz[1] != minf) val[k].maxz[1] += add;
    if(val[k].maxz[2] != minf) val[k].maxz[2] += add;
    if(val[k].maxz[3] != minf) val[k].maxz[3] += add;
    if(lz.upperx < val[k].maxx){
      val[k].upperx = lz.upperx;
      if(val[k].maxz[1] != minf) val[k].maxz[1] -= (val[k].maxx - lz.upperx);
      if(val[k].maxz[3] != minf) val[k].maxz[3] -= (val[k].maxx - lz.upperx);
      val[k].maxx = lz.upperx;
    }
    if(lz.uppery < val[k].maxy){
      val[k].uppery = lz.uppery;
      if(val[k].maxz[2] != minf) val[k].maxz[2] -= (val[k].maxy - lz.uppery);
      if(val[k].maxz[3] != minf) val[k].maxz[3] -= (val[k].maxy - lz.uppery);
      val[k].maxy = lz.uppery;
    }
  }
  void merge_val(T &v, const T &l, const T &r){
    v.maxz.fill(minf);
    v.maxz[0] = std::max(l.maxz[0], r.maxz[0]);
    if(l.maxx == r.maxx){
      v.maxx = l.maxx;
      v.smaxx = std::max(l.smaxx, r.smaxx);
      v.maxz[1] = std::max(l.maxz[1], r.maxz[1]);
    }else if(l.maxx > r.maxx){
      v.maxx = l.maxx;
      v.smaxx = std::max(l.smaxx, r.maxx);
      v.maxz[1] = l.maxz[1];
      v.maxz[0] = std::max(v.maxz[0], r.maxz[1]);
    }else{
      v.maxx = r.maxx;
      v.smaxx = std::max(l.maxx, r.smaxx);
      v.maxz[0] = std::max(v.maxz[0], l.maxz[1]);
      v.maxz[1] = r.maxz[1];
    }
    if(l.maxy == r.maxy){
      v.maxy = l.maxy;
      v.smaxy = std::max(l.smaxy, r.smaxy);
      v.maxz[2] = std::max(l.maxz[2], r.maxz[2]);
    }else if(l.maxy > r.maxy){
      v.maxy = l.maxy;
      v.smaxy = std::max(l.smaxy, r.maxy);
      v.maxz[2] = l.maxz[2];
      v.maxz[0] = std::max(v.maxz[0], r.maxz[2]);
    }else{
      v.maxy = r.maxy;
      v.smaxy = std::max(l.maxy, r.smaxy);
      v.maxz[0] = std::max(v.maxz[0], l.maxz[2]);
      v.maxz[2] = r.maxz[2];
    }
    int flag = (l.maxx == v.maxx) + ((l.maxy == v.maxy) << 1);
    v.maxz[flag] = std::max(v.maxz[flag], l.maxz[3]);
    flag = (r.maxx == v.maxx) + ((r.maxy == v.maxy) << 1);
    v.maxz[flag] = std::max(v.maxz[flag], r.maxz[3]);
  }
  void push_down(int k, int len){
    if(M - 1 <= k || (val[k].addx == 0 && val[k].addy == 0 && val[k].upperx == inf && val[k].uppery == inf)) return;
    len /= 2;
    propagate(k * 2 + 1, len, val[k]);
    propagate(k * 2 + 2, len, val[k]);
    val[k].addx = val[k].addy = 0, val[k].upperx = inf, val[k].uppery = inf;
  }
  void set(int a, const T &x, int k, int l, int r){
    if(r - l == 1){
      val[k].maxx = x.maxx;
      val[k].maxy = x.maxy;
      val[k].maxz[3] = x.maxx + x.maxy;
      return;
    }
    push_down(k, r - l);
    int mid = (l + r) >> 1;
    if(a < mid) set(a, x, k * 2 + 1, l, mid);
    else set(a, x, k * 2 + 2, mid, r);
    merge_val(val[k], val[k * 2 + 1], val[k * 2 + 2]);
  }
  std::pair<Val, Val> get(int a, int k, int l, int r){
    if(r - l == 1) return {val[k].maxx, val[k].maxy};
    push_down(k, r - l);
    int mid = (l + r) >> 1;
    if(a < mid) return get(a, k * 2 + 1, l, mid);
    else return get(a, k * 2 + 2, mid, r);
  }
  void update_chmin_x(int a, int b, const T &x, int k, int l, int r){
    if(r <= a || b <= l || val[k].maxx <= x.upperx) return;
    if(a <= l && r <= b && val[k].smaxx < x.upperx){
      propagate(k, r - l, x);
      return;
    }
    push_down(k, r - l);
    update_chmin_x(a, b, x, k * 2 + 1, l, (l + r) / 2);
    update_chmin_x(a, b, x, k * 2 + 2, (l + r) / 2, r);
    merge_val(val[k], val[k * 2 + 1], val[k * 2 + 2]);
  }
  void update_chmin_y(int a, int b, const T &x, int k, int l, int r){
    if(r <= a || b <= l || val[k].maxy <= x.uppery) return;
    if(a <= l && r <= b && val[k].smaxy < x.uppery){
      propagate(k, r - l, x);
      return;
    }
    push_down(k, r - l);
    update_chmin_y(a, b, x, k * 2 + 1, l, (l + r) / 2);
    update_chmin_y(a, b, x, k * 2 + 2, (l + r) / 2, r);
    merge_val(val[k], val[k * 2 + 1], val[k * 2 + 2]);
  }
  void update_add(int a, int b, const T &x, int k, int l, int r){
    if(r <= a || b <= l) return;
    if(a <= l && r <= b){
      propagate(k, r - l, x);
      return;
    }
    push_down(k, r - l);
    update_add(a, b, x, k * 2 + 1, l, (l + r) / 2);
    update_add(a, b, x, k * 2 + 2, (l + r) / 2, r);
    merge_val(val[k], val[k * 2 + 1], val[k * 2 + 2]);
  }
  Val query_maxsum(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return minf;
    if(a <= l && r <= b) return std::max({val[k].maxz[0], val[k].maxz[1], val[k].maxz[2], val[k].maxz[3]});
    push_down(k, r - l);
    return std::max(query_maxsum(a, b, k * 2 + 1, l, (l + r) / 2), query_maxsum(a, b, k * 2 + 2, (l + r) / 2, r));
  }
  Val query_max_x(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return minf;
    if(a <= l && r <= b) return val[k].maxx;
    push_down(k, r - l);
    return std::max(query_max_x(a, b, k * 2 + 1, l, (l + r) / 2), query_max_x(a, b, k * 2 + 2, (l + r) / 2, r));
  }
  Val query_max_y(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return minf;
    if(a <= l && r <= b) return val[k].maxy;
    push_down(k, r - l);
    return std::max(query_max_y(a, b, k * 2 + 1, l, (l + r) / 2), query_max_y(a, b, k * 2 + 2, (l + r) / 2, r));
  }
  std::pair<Val, Val> query_max_x2(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return {minf, minf};
    if(a <= l && r <= b) return {val[k].maxx, std::max(val[k].maxz[1], val[k].maxz[3]) - val[k].maxx};
    push_down(k, r - l);
    return std::max(query_max_x2(a, b, k * 2 + 1, l, (l + r) / 2), query_max_x2(a, b, k * 2 + 2, (l + r) / 2, r));
  }
  std::pair<Val, Val> query_max_y2(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return {minf, minf};
    if(a <= l && r <= b) return {val[k].maxy, std::max(val[k].maxz[2], val[k].maxz[3]) - val[k].maxy};
    push_down(k, r - l);
    return std::max(query_max_y(a, b, k * 2 + 1, l, (l + r) / 2), query_max_y(a, b, k * 2 + 2, (l + r) / 2, r));
  }
public:
  beats_chmin_add_maxsum(): M(0){}
  beats_chmin_add_maxsum(int n, Val x0, Val y0): N(n), M(1 << ceil_pow2(n)), val(2 * M - 1){
    for(int i = 0; i < n; i++){
      val[M - 1 + i].maxx = x0;
      val[M - 1 + i].maxy = y0;
      val[M - 1 + i].maxz[3] = x0 + y0;
    }
    for(int i = M - 2; i >= 0; i--) merge_val(val[i], val[i * 2 + 1], val[i * 2 + 2]);
  }
  beats_chmin_add_maxsum(const std::vector<std::pair<Val, Val>> &v): N(v.size()), M(1 << ceil_pow2(N)), val(2 * M - 1){
    for(int i = 0; i < v.size(); i++){
      val[M - 1 + i].maxx = v[i].first;
      val[M - 1 + i].maxy = v[i].second;
      val[M - 1 + i].maxz[3] = v[i].first + v[i].second;
    }
    for(int i = M - 2; i >= 0; i--) merge_val(val[i], val[i * 2 + 1], val[i * 2 + 2]);
  }
  // {xi, yi} <- {x, y}
  void set(int k, Val x, Val y){
    T lz;
    lz.maxx = x;
    lz.maxy = y;
    set(k, x, 0, 0, M);
  }
  // {xi, yi}
  std::pair<Val, Val> get(int k){
    return get(k, 0, 0, M);
  }
  // xi in [a, b)をchmin(xi, k)
  void update_chmin_x(int a, int b, Val k){
    T lz;
    lz.upperx = k;
    update_chmin_x(a, b, lz, 0, 0, M);
  }
  // yi in [a, b)をchmin(yi, k)
  void update_chmin_y(int a, int b, Val k){
    T lz;
    lz.uppery = k;
    update_chmin_y(a, b, lz, 0, 0, M);
  }
  // xi in [a, b) += x
  // yi in [a, b) += y
  void update_add(int a, int b, Val x, Val y){
    T lz;
    lz.addx = x;
    lz.addy = y;
    update_add(a, b, lz, 0, 0, M);
  }
  // xi in [0, N)をchmin(xi, k)
  void update_chmin_x_all(Val k){
    T lz;
    lz.upperx = k;
    update_chmin_x(0, N, lz, 0, 0, M);
  }
  // yi in [0, N)をchmin(yi, k)
  void update_chmin_y_all(Val k){
    T lz;
    lz.uppery = k;
    update_chmin_y(0, N, lz, 0, 0, M);
  }
  // xi in [0, N) += x
  // yi in [0, N) += y
  void update_add_all(Val x, Val y){
    T lz;
    lz.addx = x;
    lz.addy = y;
    update_add(0, N, lz, 0, 0, M);
  }
  // max(xi) in [a, b)
  Val query_max_x(int a, int b){
    return query_max_x(a, b, 0, 0, M);
  }
  // max(yi) in [a, b)
  Val query_max_y(int a, int b){
    return query_max_y(a, b, 0, 0, M);
  }
  // max({xi, yi}) in [a, b) (xが大きいものが優先)
  std::pair<Val, Val> query_max_x2(int a, int b){
    return query_max_x2(a, b, 0, 0, M);
  }
  // max({xi, yi}) in [a, b) (yが大きいものが優先)
  std::pair<Val, Val> query_max_y2(int a, int b){
    auto res = query_max_y2(a, b, 0, 0, M);
    return {res.second, res.first};
  }
  // max(xi + yi) in [a, b)
  Val query_maxsum(int a, int b){
    return query_maxsum(a, b, 0, 0, M);
  }
  // max(xi) in [0, N)
  Val query_max_x_all(){
    return query_max_x(0, M, 0, 0, M);
  }
  // max(yi) in [0, N)
  Val query_max_y_all(){
    return query_max_y(0, M, 0, 0, M);
  }
  // max(xi + yi) in [0, N)
  Val query_maxsum_all(){
    return query_maxsum(0, M, 0, 0, M);
  }
};
#endif