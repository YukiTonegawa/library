#ifndef _ABSTRACT_SEGMENT_TREE_BEATS_H_
#define _ABSTRACT_SEGMENT_TREE_BEATS_H_
#include <vector>
#include <algorithm>
#include <cassert>

template<typename beats_struct>
struct segment_tree_beats{
private:
  using Val = typename beats_struct::Val;
  using Lazy = typename beats_struct::Lazy;
  int N, M;
  std::vector<Val> val;
  std::vector<Lazy> lazy;
  std::vector<bool> is_id;

  void push_down(int k, int l, int r){
    if(M - 1 <= k || is_id[k]) return;
    int mid = (l + r) >> 1;
    propagate(k * 2 + 1, l, mid, lazy[k]);
    propagate(k * 2 + 2, mid, r, lazy[k]);
    is_id[k] = true;
  }
  void propagate(int k, int l, int r, const Lazy &x){
    if(k < M - 1){
      if(is_id[k]){
        lazy[k] = x;
        is_id[k] = false;
      }else{
        beats_struct::propagate_lazy(lazy[k], x);
      }
    }
    beats_struct::apply(val[k], x, l, r);
  }
  void set(int a, Val x, int k, int l, int r){
    if(r - l == 1){
      val[k] = x;
      return;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid) set(a, x, k * 2 + 1, l, mid);
    else set(a, x, k * 2 + 2, mid, r);
    beats_struct::merge_val(val[k], val[k * 2 + 1], val[k * 2 + 2]);
  }
  void get(int a, int k, int l, int r, Val &ans){
    if(r - l == 1){
      ans = val[k];
      return;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid) get(a, k * 2 + 1, l, mid, ans);
    else get(a, k * 2 + 2, mid, r, ans);
  }
  template<int id>
  void update(int a, int b, const Lazy &x, int k, int l, int r){
    if(r <= a || b <= l) return;
    if(r <= N && beats_struct::template break_check<id>(val[k], x)) return;
    if(a <= l && r <= b && beats_struct::template tag_check<id>(val[k], x)){
      propagate(k, l, r, x);
      return;
    }
    push_down(k, l, r);
    update<id>(a, b, x, k * 2 + 1, l, (l + r) / 2);
    update<id>(a, b, x, k * 2 + 2, (l + r) / 2, r);
    beats_struct::merge_val(val[k], val[k * 2 + 1], val[k * 2 + 2]);
  }
  // lazyをpropagateする直前に変形するタイプ
  template<int id>
  void update2(int a, int b, const Lazy &x, int k, int l, int r){
    if(r <= a || b <= l) return;
    if(r <= N && beats_struct::template break_check<id>(val[k], x)) return;
    if(a <= l && r <= b && beats_struct::template tag_check<id>(val[k], x)){
      Lazy y = beats_struct::template transform_lazy<id>(val[k], x);
      propagate(k, l, r, y);
      return;
    }
    push_down(k, l, r);
    update2<id>(a, b, x, k * 2 + 1, l, (l + r) / 2);
    update2<id>(a, b, x, k * 2 + 2, (l + r) / 2, r);
    beats_struct::merge_val(val[k], val[k * 2 + 1], val[k * 2 + 2]);
  }
  void query(int a, int b, int k, int l, int r, Val &ans){
    if(r <= a || b <= l) return;
    if(a <= l && r <= b){
      Val tmp = beats_struct::id_val();
      beats_struct::merge_val(tmp, ans, val[k]);
      ans = tmp;
      return;
    }
    push_down(k, l, r);
    query(a, b, k * 2 + 1, l, (l + r) / 2, ans);
    query(a, b, k * 2 + 2, (l + r) / 2, r, ans);
  }
  int c2(int x){
    int y = 1;
    while(y < x) y <<= 1;
    return y;
  }
public:
  segment_tree_beats(): N(0), M(0){}
  template<typename T>
  segment_tree_beats(const std::vector<T> &v): N(v.size()), M(c2(N)), val(2 * M - 1, Val()), lazy(M - 1, Lazy()), is_id(M - 1, true){
    for(int i = 0; i < N; i++) val[M - 1 + i] = Val(v[i]);
    for(int i = N; i < M; i++) val[M - 1 + i] = beats_struct::id_val();
    for(int i = M - 2; i >= 0; i--) beats_struct::merge_val(val[i], val[i * 2 + 1], val[i * 2 + 2]);
  }
  void set(int k, Val x){
    set(k, x, 0, 0, M);
  }
  Val get(int k){
    Val ans = beats_struct::id_val();
    get(k, 0, 0, M, ans);
    return ans;
  }
  template<int id>
  void update(int l, int r, Lazy x){
    update<id>(l, r, x, 0, 0, M);
  }
  template<int id>
  void update2(int l, int r, Lazy x){
    update2<id>(l, r, x, 0, 0, M);
  }
  Val query(int l, int r){
    Val ans = beats_struct::id_val();
    query(l, r, 0, 0, M, ans);
    return ans;
  }
};
/*
template<typename T>
struct add_chmin_chmax{
  static constexpr T inf = std::numeric_limits<T>::max();
  static constexpr T minf = std::numeric_limits<T>::min();
  struct Val{
    T min, second_min, max, second_max, sum;
    int min_cnt, max_cnt;
    Val(): min(inf), second_min(inf), max(minf), second_max(minf), sum(0), min_cnt(0), max_cnt(0){}
    Val(T x): min(x), second_min(inf), max(x), second_max(minf), sum(x), min_cnt(1), max_cnt(1){}
  };
  struct Lazy{
    T add, lower, upper;
    Lazy(): add(0), lower(minf), upper(inf){}
    Lazy(T a, T b, T c): add(a), lower(b), upper(c){}
  };
  static Val id_val(){
    return Val();
  }
  // l, rをマージしてvに代入
  static void merge_val(Val &v, const Val &l, const Val &r){
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
    if(l.min == r.min){
      v.min = l.min;
      v.min_cnt = l.min_cnt + r.min_cnt;
      v.second_min = std::min(l.second_min, r.second_min);
    }else if(l.min < r.min){
      v.min = l.min;
      v.min_cnt = l.min_cnt;
      v.second_min = std::min(l.second_min, r.min);
    }else{
      v.min = r.min;
      v.min_cnt = r.min_cnt;
      v.second_min = std::min(l.min, r.second_min);
    }
  }
  static void apply(Val &a, const Lazy &b, int l, int r){
    if(b.add){
      a.min += b.add;
      if(a.second_min != inf) a.second_min += b.add;
      a.max += b.add;
      if(a.second_max != minf) a.second_max += b.add;
      a.sum += b.add * (r - l);
    }
    if(a.min < b.lower){
      a.sum += (b.lower - a.min) * a.min_cnt;
      if(a.second_max == a.min) a.second_max = b.lower;
      else if(a.max == a.min) a.max = b.lower, a.second_max = minf;
      a.min = b.lower;
    }
    if(b.upper < a.max){
      a.sum -= (a.max - b.upper) * a.max_cnt;
      if(a.second_min == a.max) a.second_min = b.upper;
      else if(a.min == a.max) a.min = b.upper, a.second_min = inf;
      a.max = b.upper;
    }
  }
  static void propagate_lazy(Lazy &a, const Lazy &b){
    if(b.add){
      a.add += b.add;
      if(a.lower != minf) a.lower += b.add;
      if(a.upper != inf) a.upper += b.add;
    }
    if(a.upper <= b.lower) a.lower = a.upper = b.lower;
    else if(b.upper <= a.lower) a.lower = a.upper = b.upper;
    else{
      a.lower = std::max(a.lower, b.lower);
      a.upper = std::min(a.upper, b.upper);
    }
  }
  // chmin
  template<int id, std::enable_if_t<id == 0>* = nullptr>
  static bool break_check(const Val &v, const Lazy &x){
    return v.max < x.upper;
  }
  template<int id, std::enable_if_t<id == 0>* = nullptr>
  static bool tag_check(const Val &v, const Lazy &x){
    return v.second_max < x.upper;
  }
  // chmax
  template<int id, std::enable_if_t<id == 1>* = nullptr>
  static bool break_check(const Val &v, const Lazy &x){
    return v.min > x.lower;
  }
  template<int id, std::enable_if_t<id == 1>* = nullptr>
  static bool tag_check(const Val &v, const Lazy &x){
    return v.second_min > x.lower;
  }
  // add
  template<int id, std::enable_if_t<id == 2>* = nullptr>
  static bool break_check(const Val &v, const Lazy &x){
    return false;
  }
  template<int id, std::enable_if_t<id == 2>* = nullptr>
  static bool tag_check(const Val &v, const Lazy &x){
    return true;
  }
};


template<typename T>
struct abc256ex_set_div_sum{
  static constexpr T minf = std::numeric_limits<T>::min();
  struct Val{
    T sum, unique; // unique != minfなら区間がその値でユニーク
    Val(): sum(0), unique(minf){}
    Val(T x): sum(x), unique(x){}
  };
  struct Lazy{
    T setx; 
    Lazy(): setx(minf){}
    Lazy(T x): setx(x){}
  };
  static Val id_val(){
    return Val();
  }
  // l, rをマージしてvに代入
  static void merge_val(Val &v, const Val &l, const Val &r){
    v.sum = l.sum + r.sum;
    if(l.unique == r.unique) v.unique = l.unique;
    else v.unique = minf;
  }
  static void apply(Val &a, const Lazy &b, int l, int r){
    a.sum = b.setx * (r - l);
    a.unique = b.setx;
  }
  static void propagate_lazy(Lazy &a, const Lazy &b){
    a = b;
  }
  // set
  template<int id, std::enable_if_t<id == 0>* = nullptr>
  static bool break_check(const Val &v, const Lazy &x){
    return false;
  }
  template<int id, std::enable_if_t<id == 0>* = nullptr>
  static bool tag_check(const Val &v, const Lazy &x){
    return true;
  }
  // div
  template<int id, std::enable_if_t<id == 1>* = nullptr>
  static bool break_check(const Val &v, const Lazy &x){
    return v.unique == 0;
  }
  template<int id, std::enable_if_t<id == 1>* = nullptr>
  static bool tag_check(const Val &v, const Lazy &x){
    return v.unique != minf;
  }
  template<int id, std::enable_if_t<id == 1>* = nullptr>
  static Lazy transform_lazy(const Val &v, const Lazy &x){
    return Lazy(v.unique / x.setx);
  }
};
template<typename T>
struct set_mod_min_max_sum{
  static constexpr T inf = std::numeric_limits<T>::max();
  static constexpr T minf = std::numeric_limits<T>::min();
  struct Val{
    T sum, min, max;
    Val(): sum(0), min(inf), max(minf){}
    Val(T x): sum(x), min(x), max(x){}
  };
  struct Lazy{
    T setx; 
    Lazy(): setx(minf){}
    Lazy(T x): setx(x){}
  };
  static Val id_val(){
    return Val();
  }
  // l, rをマージしてvに代入
  static void merge_val(Val &v, const Val &l, const Val &r){
    v.sum = l.sum + r.sum;
    v.min = std::min(l.min, r.min);
    v.max = std::max(l.max, r.max);
  }
  static void apply(Val &a, const Lazy &b, int l, int r){
    a.sum = b.setx * (r - l);
    a.min = a.max = b.setx;
  }
  static void propagate_lazy(Lazy &a, const Lazy &b){
    a = b;
  }
  // set
  template<int id, std::enable_if_t<id == 0>* = nullptr>
  static bool break_check(const Val &v, const Lazy &x){
    return false;
  }
  template<int id, std::enable_if_t<id == 0>* = nullptr>
  static bool tag_check(const Val &v, const Lazy &x){
    return true;
  }
  // mod
  template<int id, std::enable_if_t<id == 1>* = nullptr>
  static bool break_check(const Val &v, const Lazy &x){
    return v.max < x.setx;
  }
  template<int id, std::enable_if_t<id == 1>* = nullptr>
  static bool tag_check(const Val &v, const Lazy &x){
    return v.min == v.max;
  }
  template<int id, std::enable_if_t<id == 1>* = nullptr>
  static Lazy transform_lazy(const Val &v, const Lazy &x){
    return Lazy(v.min % x.setx);
  }
};
*/

#endif