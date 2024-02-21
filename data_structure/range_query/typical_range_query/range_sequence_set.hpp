#ifndef _RANGE_SEQUENCE_SET_H_
#define _RANGE_SEQUENCE_SET_H_
#include <vector>
#include <cassert>
#include <limits>

template<typename Sequence>
struct range_sequence_set{
private:
  using Val = typename Sequence::Val;
  using Lazy = typename Sequence::Lazy;
  struct node{
    Val sum;
    bool is_id;
    Lazy lazy;
  };
  int N, M;
  std::vector<node> v;
  static constexpr int ceil_pow2(int n){
    int m = 1;
    while(m < n) m <<= 1;
    return m;
  }
  void push_down(int k, int l, int r){
    if(v[k].is_id) return;
    int mid = (l + r) >> 1;
    auto [y, z] = Sequence::split_lazy(v[k].lazy, l, mid, r);
    propagate(k * 2 + 1, l, mid, y);
    propagate(k * 2 + 2, mid, r, z);
    v[k].is_id = true;
  }
  void propagate(int k, int l, int r, Lazy x){
    v[k].sum = Sequence::get_sum(x, l, r);
    v[k].is_id = false;
    v[k].lazy = x;
  }
  void set(int a, Val x, int k, int l, int r){
    if(r - l == 1){
      v[k].sum = x;
      return;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid) set(a, x, k * 2 + 1, l, mid);
    else set(a, x, k * 2 + 2, mid, r);
    v[k].sum = Sequence::merge(v[k * 2 + 1].sum, v[k * 2 + 2].sum);
  }
  Val get(int a, int k, int l, int r){
    if(r - l == 1) return v[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid) return get(a, k * 2 + 1, l, mid);
    else return get(a, k * 2 + 2, mid, r);
  }
  void update(int a, int b, Lazy x, int k, int l, int r){
    if(r <= a || b <= l) return;
    if(a <= l && r <= b){
      propagate(k, l, r, x);
      return;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    auto [y, z] = Sequence::split_lazy(x, l, mid, r);
    update(a, b, y, k * 2 + 1, l, mid);
    update(a, b, z, k * 2 + 2, mid, r);
    v[k].sum = Sequence::merge(v[k * 2 + 1].sum, v[k * 2 + 2].sum);
  }
  Val query(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return Sequence::id();
    if(a <= l && r <= b) return v[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return Sequence::merge(query(a, b, k * 2 + 1, l, mid), query(a, b, k * 2 + 2, mid, r));
  }
public:
  range_sequence_set(){}
  range_sequence_set(int n, Val val): N(n), M(ceil_pow2(n)), v(M * 2 - 1){
    for(int i = 0; i < N; i++) v[M - 1 + i].sum = val;
    for(int i = N; i < M; i++) v[M - 1 + i].sum = Sequence::id();
    for(int i = M - 2; i >= 0; i--) v[i].sum = Sequence::merge(v[i * 2 + 1].sum, v[i * 2 + 2].sum);
    for(int i = 0; i < 2 * M - 1; i++) v[i].is_id = true;
  }
  template<typename T>
  range_sequence_set(const std::vector<T> &val): N(val.size()),  M(ceil_pow2(N)), v(M * 2 - 1){
    for(int i = 0; i < N; i++) v[M - 1 + i].sum = Val(val[i]);
    for(int i = N; i < M; i++) v[M - 1 + i].sum = Sequence::id();
    for(int i = M - 2; i >= 0; i--) v[i].sum = Sequence::merge(v[i * 2 + 1].sum, v[i * 2 + 2].sum);
    for(int i = 0; i < 2 * M - 1; i++) v[i].is_id = true;
  }
  // val[k] <- x
  void set(int k, Val x){
    set(k, x, 0, 0, M);
  }
  // val[k]
  Val get(int k){
    return get(k, 0, 0, M);
  }
  // [l, r)を更新
  void update(int l, int r, Lazy x){
    update(l, r, x, 0, 0, M);
  }
  // sum[a, b)
  Val query(int a, int b){
    return query(a, b, 0, 0, M);
  }
  Val query_all(){
    return v[0].sum;
  }
};
template<typename T>
struct arithmetic_set_min{
  using Val = T; // 値の型
  using Lazy = std::pair<T, T>; // 数列を表す型
  static Val id(){
    return std::numeric_limits<T>::max();
  }
  static Val merge(Val a, Val b){
    return std::min(a, b);
  }
  // xを区間[l, r)に作用させた時の区間のsum
  static Val get_sum(Lazy x, int l, int r){
    return (x.first < 0 ? (r - 1) * x.first : l * x.first) + x.second;
  }
  // [l, r)にxを作用させることと[l, m)にy, [m, r)にzを作用させることが等価になるような作用素{y, z}
  static std::pair<Lazy, Lazy> split_lazy(Lazy x, int l, int m, int r){
    return {x, x};
  }
};
template<typename T>
struct arithmetic_set_max{
  using Val = T; // 値の型
  using Lazy = std::pair<T, T>; // 数列を表す型
  static Val id(){
    return std::numeric_limits<T>::min();
  }
  static Val merge(Val a, Val b){
    return std::max(a, b);
  }
  // xを区間[l, r)に作用させた時の区間のsum
  static Val get_sum(Lazy x, int l, int r){
    return (x.first >= 0 ? (r - 1) * x.first : l * x.first) + x.second;
  }
  // [l, r)にxを作用させることと[l, m)にy, [m, r)にzを作用させることが等価になるような作用素{y, z}
  static std::pair<Lazy, Lazy> split_lazy(Lazy x, int l, int m, int r){
    return {x, x};
  }
};

template<typename T>
struct arithmetic_set_sum{
  using Val = T; // 値の型
  using Lazy = std::pair<T, T>; // 数列を表す型
  static Val id(){
    return 0;
  }
  static Val merge(Val a, Val b){
    return a + b;
  }
  // xを区間[l, r)に作用させた時の区間のsum
  static Val get_sum(Lazy x, int l, int r){
    return x.first * (((Val)r * (r - 1) - ((Val)l * (l - 1))) / 2) + x.second * (Val)(r - l);
  }
  // [l, r)にxを作用させることと[l, m)にy, [m, r)にzを作用させることが等価になるような作用素{y, z}
  static std::pair<Lazy, Lazy> split_lazy(Lazy x, int l, int m, int r){
    return {x, x};
  }
};

template<typename mint>
std::pair<mint, mint> composite_pow(int k, mint a, mint b){
  mint A = 1, B = 0;
  while(k){
    if(k & 1) B += A * b, A *= a;
    k >>= 1;
    b += a * b;
    a *= a;
  }
  return {A, B};
}
template<typename T>
struct range_set_composite{
  using Val = T; // 値の型
  using Lazy = T; // 数列を表す型
  static Val id(){
    return {1, 0};
  }
  static Val merge(Val a, Val b){
    return {a.first * b.first, b.first * a.second + b.second};
  }
  // xを区間[l, r)に作用させた時の区間のsum
  static Val get_sum(Lazy x, int l, int r){
    return composite_pow(r - l, x.first, x.second);
  }
  // [l, r)にxを作用させることと[l, m)にy, [m, r)にzを作用させることが等価になるような作用素{y, z}
  static std::pair<Lazy, Lazy> split_lazy(Lazy x, int l, int m, int r){
    return {x, x};
  }
};
#endif