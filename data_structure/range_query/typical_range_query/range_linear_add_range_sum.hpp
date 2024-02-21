#ifndef _RANGE_ARITHMETIC_ADD_RANGE_SUM_H_
#define _RANGE_ARITHMETIC_ADD_RANGE_SUM_H_
#include <vector>
template<typename Val>
struct range_linear_add_range_sum{
private:
  using Lazy = std::pair<Val, Val>;
  Lazy z = {0, 0};
  struct node{
    Val sum;
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
    if(v[k].lazy == z) return;
    int mid = (l + r) >> 1;
    propagate(k * 2 + 1, l, mid, v[k].lazy);
    propagate(k * 2 + 2, mid, r, {v[k].lazy.first, v[k].lazy.first * (mid - l) + v[k].lazy.second});
    v[k].lazy = z;
  }
  void propagate(int k, int l, int r, Lazy x){
    r -= l;
    v[k].sum += x.first * ((Val)r * (r - 1) / 2) + x.second * r;
    v[k].lazy.first += x.first, v[k].lazy.second += x.second;
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
    v[k].sum = v[k * 2 + 1].sum + v[k * 2 + 2].sum;
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
    update(a, b, x, k * 2 + 1, l, mid);
    update(a, b, {x.first, x.second + x.first * (mid - l)}, k * 2 + 2, mid, r);
    v[k].sum = v[k * 2 + 1].sum + v[k * 2 + 2].sum;
  }
  Val query(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return 0;
    if(a <= l && r <= b) return v[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return query(a, b, k * 2 + 1, l, mid) + query(a, b, k * 2 + 2, mid, r);
  }
public:
  range_linear_add_range_sum(){}
  range_linear_add_range_sum(int n, Val val): N(n), M(ceil_pow2(n)), v(M * 2 - 1){
    for(int i = 0; i < N; i++) v[M - 1 + i].sum = val;
    for(int i = N; i < M; i++) v[M - 1 + i].sum = 0;
    for(int i = M - 2; i >= 0; i--) v[i].sum = v[i * 2 + 1].sum + v[i * 2 + 2].sum;
    for(int i = 0; i < 2 * M - 1; i++) v[i].lazy = z;
  }
  template<typename T>
  range_linear_add_range_sum(const std::vector<T> &val): N(val.size()),  M(ceil_pow2(N)), v(M * 2 - 1){
    for(int i = 0; i < N; i++) v[M - 1 + i].sum = Val(val[i]);
    for(int i = N; i < M; i++) v[M - 1 + i].sum = 0;
    for(int i = M - 2; i >= 0; i--) v[i].sum = v[i * 2 + 1].sum + v[i * 2 + 2].sum;
    for(int i = 0; i < 2 * M - 1; i++) v[i].lazy = z;
  }
  // val[k] <- x
  void set(int k, Val x){
    set(k, x, 0, 0, M);
  }
  // val[k]
  Val get(int k){
    return get(k, 0, 0, M);
  }
  // [l, r)にa * k + bを足す
  // kは0から始める
  // {b, a + b, 2a + b, ..., (r - l - 1)a + b}
  void update(int l, int r, Val a, Val b){
    update(l, r, {a, b - a * l}, 0, 0, M);
  }
  // sum[a, b)
  Val query(int a, int b){
    return query(a, b, 0, 0, M);
  }
  Val query_all(){
    return v[0].sum;
  }
};
#endif