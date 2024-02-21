#ifndef _BEATS_CHMIN_ADD_H_
#define _BEATS_CHMIN_ADD_H_
#include <vector>
#include <limits>
#include <algorithm>
#include <cassert>
// Val := 1要素の値を扱う型
// ValSum := 区間和を扱う型
// ならしO(log^2N)
template<typename Val, typename ValSum>
struct beats_chmin_add{
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
    Val max, second_max;
    ValSum sum;
    int max_cnt;
    Val add, upper;
    T(): max(minf), second_max(minf), sum(0), max_cnt(0), add(0), upper(inf){}
  };
  std::vector<T> val;
  void propagate(int k, int len, const T &lz){
    if(lz.add){
      val[k].add += lz.add;
      if(val[k].upper != inf) val[k].upper += lz.add;
      val[k].max += lz.add;
      if(val[k].second_max != minf) val[k].second_max += lz.add;
      val[k].sum += (ValSum)lz.add * len;
    }
    if(lz.upper < val[k].max){
      val[k].upper = lz.upper;
      val[k].sum -= (ValSum)(val[k].max - lz.upper) * val[k].max_cnt;
      val[k].max = lz.upper;
    }
  }
  void merge_val(T &v, const T &l, const T &r){
    v.sum = l.sum + r.sum;
    if(l.max == r.max){
      v.max = l.max;
      v.max_cnt = l.max_cnt + r.max_cnt;
      v.second_max = std::max(l.second_max, r.second_max);
    }else if(l.max > r.max){
      v.max = l.max;
      v.max_cnt = l.max_cnt;
      v.second_max = std::max(l.second_max, r.max);
    }else{
      v.max = r.max;
      v.max_cnt = r.max_cnt;
      v.second_max = std::max(l.max, r.second_max);
    }
  }
  void push_down(int k, int len){
    if(M - 1 <= k || (val[k].add == 0 && val[k].upper == inf)) return;
    len /= 2;
    propagate(k * 2 + 1, len, val[k]);
    propagate(k * 2 + 2, len, val[k]);
    val[k].add = 0, val[k].upper = inf;
  }
  void set(int a, Val x, int k, int l, int r){
    if(r - l == 1){
      val[k].max = val[k].sum = x;
      return;
    }
    push_down(k, r - l);
    int mid = (l + r) >> 1;
    if(a < mid) set(a, x, k * 2 + 1, l, mid);
    else set(a, x, k * 2 + 2, mid, r);
    merge_val(val[k], val[k * 2 + 1], val[k * 2 + 2]);
  }
  Val get(int a, int k, int l, int r){
    if(r - l == 1) return val[k].max;
    push_down(k, r - l);
    int mid = (l + r) >> 1;
    if(a < mid) return get(a, k * 2 + 1, l, mid);
    else return get(a, k * 2 + 2, mid, r);
  }
  void update_chmin(int a, int b, const T &x, int k, int l, int r){
    if(r <= a || b <= l || val[k].max <= x.upper) return;
    if(a <= l && r <= b && val[k].second_max < x.upper){
      propagate(k, r - l, x);
      return;
    }
    push_down(k, r - l);
    update_chmin(a, b, x, k * 2 + 1, l, (l + r) / 2);
    update_chmin(a, b, x, k * 2 + 2, (l + r) / 2, r);
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
  Val query_max(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return minf;
    if(a <= l && r <= b) return val[k].max;
    push_down(k, r - l);
    return std::max(query_max(a, b, k * 2 + 1, l, (l + r) / 2), query_max(a, b, k * 2 + 2, (l + r) / 2, r));
  }
  ValSum query_sum(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return 0;
    if(a <= l && r <= b) return val[k].sum;
    push_down(k, r - l);
    return query_sum(a, b, k * 2 + 1, l, (l + r) / 2) + query_sum(a, b, k * 2 + 2, (l + r) / 2, r);
  }
public:
  beats_chmin_add(): M(0){}
  beats_chmin_add(int n, Val x): N(n), M(1 << ceil_pow2(n)), val(2 * M - 1){
    for(int i = 0; i < n; i++){
      val[M - 1 + i].max = val[M - 1 + i].sum = x;
      val[M - 1 + i].max_cnt = 1;
    }
    for(int i = M - 2; i >= 0; i--) merge_val(val[i], val[i * 2 + 1], val[i * 2 + 2]);
  }
  template<typename H>
  beats_chmin_add(const std::vector<H> &v): N(v.size()), M(1 << ceil_pow2(N)), val(2 * M - 1){
    for(int i = 0; i < v.size(); i++){
      val[M - 1 + i].max = val[M - 1 + i].sum = v[i];
      val[M - 1 + i].max_cnt = 1;
    }
    for(int i = M - 2; i >= 0; i--) merge_val(val[i], val[i * 2 + 1], val[i * 2 + 2]);
  }
  void set(int k, Val x){
    set(k, x, 0, 0, M);
  }
  Val get(int k){
    return get(k, 0, 0, M);
  }
  void update_chmin(int a, int b, Val x){
    T lz;
    lz.upper = x;
    update_chmin(a, b, lz, 0, 0, M);
  }
  void update_add(int a, int b, Val x){
    T lz;
    lz.add = x;
    update_add(a, b, lz, 0, 0, M);
  }
  void update_chmin_all(Val x){
    T lz;
    lz.upper = x;
    update_chmin(0, N, lz, 0, 0, M);
  }
  void update_add_all(Val x){
    T lz;
    lz.add = x;
    update_add(0, N, lz, 0, 0, M);
  }
  Val query_max(int a, int b){
    return query_max(a, b, 0, 0, M);
  }
  ValSum query_sum(int a, int b){
    return query_sum(a, b, 0, 0, M);
  }
  Val query_max_all(){
    return query_max(0, M, 0, 0, M);
  }
  ValSum query_sum_all(){
    return query_sum(0, M, 0, 0, M);
  }
};
#endif