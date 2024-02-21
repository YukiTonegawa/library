#ifndef _RANGE_BIT_OPERATION_H_
#define _RANGE_BIT_OPERATION_H_

#include <vector>
#include <array>
#include <cassert>

template<typename Val>
struct range_bit_operation_plain{
private:
  static constexpr Val inf = ~(Val)0;
  struct Lazy{
    Val __or, __and, __xor;
    Lazy(): __or(0), __and(inf), __xor(0){}
    Lazy(Val a, Val b, Val c): __or(a), __and(b), __xor(c){}
    bool is_id(){return __or == 0 && __and == inf && __xor == 0;}
    void reset(){__or = 0, __and = inf, __xor = 0;}
  };
  struct node{
    Val val;
    Lazy lazy;
    node(): lazy(){}
  };
  int N, M;
  static constexpr int ceil_pow2(int n){
    int m = 1;
    while(m < n) m <<= 1;
    return m;
  }
  std::vector<node> V;
  Val or_and_xor(Val i, const Lazy &lz){
    return ((i | lz.__or) & lz.__and) ^ lz.__xor;
  }
  void propagate_lazy(Lazy &l, const Lazy &r){
    // 例外の部分だけ, {l.or, l.and, l.xor ^ r.xor}にしたい
    Val __exception_flag = (~r.__or) & r.__and;
    Val a = __exception_flag & l.__or;
    Val b = __exception_flag & l.__and;
    Val c = __exception_flag & (l.__xor ^ r.__xor);
    l.__or  = ((l.__or & r.__or) ^ l.__xor) | r.__or;
    l.__and = ((l.__and ^ l.__xor) | r.__or) & r.__and;
    l.__xor = r.__xor;
    __exception_flag = ~__exception_flag;
    l.__or  = (l.__or  & __exception_flag) | a;
    l.__and = (l.__and & __exception_flag) | b;
    l.__xor = (l.__xor & __exception_flag) | c;
  }
  void push_down(int k, int l, int r){
    if(V[k].lazy.is_id()) return;
    int mid = (l + r) >> 1;
    propagate(k * 2 + 1, l, mid, V[k].lazy);
    propagate(k * 2 + 2, mid, r, V[k].lazy);
    V[k].lazy.reset();
  }
  void propagate(int k, int l, int r, Lazy x){
    if(M - 1 <= k && k < N + M - 1) V[k].val = or_and_xor(V[k].val, x);
    else propagate_lazy(V[k].lazy, x);
  }
  void set(int a, Val x, int k, int l, int r){
    if(r - l == 1){
      V[k].val = x;
      return;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid){
      set(a, x, k * 2 + 1, l, mid);
    }else{
      set(a, x, k * 2 + 2, mid, r);
    }
  }
  Val get(int a, int k, int l, int r){
    if(r - l == 1){
      return V[k].val;
    }
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
    update(a, b, x, k * 2 + 2, mid, r);
  }
public:
  range_bit_operation_plain(): N(0), M(0){}
  range_bit_operation_plain(int n): N(n), M(ceil_pow2(N)), V(2 * M - 1, node()){}
  template<typename T>
  range_bit_operation_plain(const std::vector<T> &v): N(v.size()), M(ceil_pow2(N)), V(2 * M - 1){
    for(int i = 0; i < N; i++) V[i + M - 1].val = v[i];
  }
  int size(){
    return N;
  }
  // val[k] <- x, O(logN)
  void set(int k, Val x){
    set(k, x, 0, 0, M);
  }
  // val[k], O(logN)
  Val get(int k){
    return get(k, 0, 0, M);
  }
  // O(logN)
  void update_or(int l, int r, Val x){
    assert(0 <= l && r <= size());
    update(l, r, Lazy(x, inf, 0), 0, 0, M);
  }
  // O(logN)
  void update_and(int l, int r, Val x){
    assert(0 <= l && r <= size());
    update(l, r, Lazy(0, x, 0), 0, 0, M);
  }
  // O(logN)
  void update_xor(int l, int r, Val x){
    assert(0 <= l && r <= size());
    update(l, r, Lazy(0, inf, x), 0, 0, M);
  }
  // O(logN)
  void update_set(int l, int r, Val x){
    assert(0 <= l && r <= size());
    update(l, r, Lazy(x, x, 0), 0, 0, M);
  }
  // O(1)
  void update_or_all(Val x){
    update(0, M, Lazy(x, inf, 0), 0, 0, M);
  }
  // O(1)
  void update_and_all(Val x){
    update(0, M, Lazy(0, x, 0), 0, 0, M);
  }
  // O(1)
  void update_xor_all(Val x){
    update(0, M, Lazy(0, inf, x), 0, 0, M);
  }
};

template<typename Val = unsigned int, typename ValSum = unsigned long long, int bitlen = 30>
struct range_bit_operation{
private:
  static constexpr Val inf = ~(Val)0;
  struct Lazy{
    Val __or, __and, __xor;
    Lazy(): __or(0), __and(inf), __xor(0){}
    Lazy(Val a, Val b, Val c): __or(a), __and(b), __xor(c){}
    bool is_id(){return __or == 0 && __and == inf && __xor == 0;}
    void reset(){__or = 0, __and = inf, __xor = 0;}
  };
  struct node{
    std::array<int, bitlen> sum;
    Lazy lazy;
    node(): lazy(){sum.fill(0);}
  };
  int N, M;
  static constexpr int ceil_pow2(int n){
    int m = 1;
    while(m < n) m <<= 1;
    return m;
  }
  std::vector<node> V;
  Val or_and_xor(Val i, const Lazy &lz){
    return ((i | lz.__or) & lz.__and) ^ lz.__xor;
  }
  void propagate_lazy(Lazy &l, const Lazy &r){
    // 例外の部分だけ, {l.or, l.and, l.xor ^ r.xor}にしたい
    Val __exception_flag = (~r.__or) & r.__and;
    Val a = __exception_flag & l.__or;
    Val b = __exception_flag & l.__and;
    Val c = __exception_flag & (l.__xor ^ r.__xor);
    l.__or  = ((l.__or & r.__or) ^ l.__xor) | r.__or;
    l.__and = ((l.__and ^ l.__xor) | r.__or) & r.__and;
    l.__xor = r.__xor;
    __exception_flag = ~__exception_flag;
    l.__or  = (l.__or  & __exception_flag) | a;
    l.__and = (l.__and & __exception_flag) | b;
    l.__xor = (l.__xor & __exception_flag) | c;
  }
  void push_down(int k, int l, int r){
    if(V[k].lazy.is_id()) return;
    int mid = (l + r) >> 1;
    propagate(k * 2 + 1, l, mid, V[k].lazy);
    propagate(k * 2 + 2, mid, r, V[k].lazy);
    V[k].lazy.reset();
  }
  void propagate(int k, int l, int r, Lazy x){
    Val zero = or_and_xor(0, x);
    Val one = or_and_xor(inf, x);
    r = std::min(r, N); // 全体更新をした時などにN以上のインデックスが更新されないようrを丸める
    int len = r - l;
    for(int i = 0; i < bitlen; i++, zero >>= 1, one >>= 1){
      V[k].sum[i] = (zero & 1) * (len - V[k].sum[i]) + (one & 1) * V[k].sum[i];
    }
    propagate_lazy(V[k].lazy, x);
  }
  void set(int a, Val x, int k, int l, int r){
    if(r - l == 1){
      for(int i = 0; i < bitlen; i++, x >>= 1) V[k].sum[i] = (x & 1);
      return;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid){
      set(a, x, k * 2 + 1, l, mid);
    }else{
      set(a, x, k * 2 + 2, mid, r);
    }
    eval(k);
  }
  Val get(int a, int k, int l, int r){
    if(r - l == 1){
      Val ret = 0;
      for(int i = 0; i < bitlen; i++) ret += V[k].sum[i] << i;
      return ret;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid) return get(a, k * 2 + 1, l, mid);
    else return get(a, k * 2 + 2, mid, r);
  }
  void eval(int k){
    assert(k < M - 1);
    for(int i = 0; i < bitlen; i++) V[k].sum[i] = V[k * 2 + 1].sum[i] + V[k * 2 + 2].sum[i];
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
    update(a, b, x, k * 2 + 2, mid, r);
    eval(k);
  }
  Val query_or(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return 0;
    if(a <= l && r <= b){
      Val ret = 0;
      for(int i = 0; i < bitlen; i++){
        ret += (Val)(V[k].sum[i] != 0) << i;
      }
      return ret;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return query_or(a, b, k * 2 + 1, l, mid) | query_or(a, b, k * 2 + 2, mid, r);
  }
  Val query_and(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return inf;
    if(a <= l && r <= b){
      Val ret = 0;
      for(int i = 0; i < bitlen; i++){
        ret += (Val)(V[k].sum[i] == r - l) << i;
      }
      return ret;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return query_and(a, b, k * 2 + 1, l, mid) & query_and(a, b, k * 2 + 2, mid, r);
  }
  Val query_xor(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return 0;
    if(a <= l && r <= b){
      Val ret = 0;
      for(int i = 0; i < bitlen; i++){
        ret += (Val)(V[k].sum[i] & 1) << i;
      }
      return ret;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return query_xor(a, b, k * 2 + 1, l, mid) ^ query_xor(a, b, k * 2 + 2, mid, r);
  }
  ValSum query_sum(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return 0;
    if(a <= l && r <= b){
      ValSum ret = 0;
      for(int i = 0; i < bitlen; i++){
        ret += (ValSum)(V[k].sum[i]) * ((Val)1 << i);
      }
      return ret;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return query_sum(a, b, k * 2 + 1, l, mid) + query_sum(a, b, k * 2 + 2, mid, r);
  }
  int query_popcount(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return 0;
    if(a <= l && r <= b){
      int ret = 0;
      for(int i = 0; i < bitlen; i++){
        ret += V[k].sum[i];
      }
      return ret;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return query_popcount(a, b, k * 2 + 1, l, mid) + query_popcount(a, b, k * 2 + 2, mid, r);
  }
  template<typename F>
  int bisect_from_left(int a, const F &f, std::array<int, bitlen> &ok, int k, int l, int r){
    if(r <= a) return -1;
    if(k < M - 1) push_down(k, l, r);
    if(a <= l){
      for(int i = 0; i < bitlen; i++) ok[i] += V[k].sum[i];
      if(!f(ok)) return -1;
      if(k >= M - 1) return l;
      for(int i = 0; i < bitlen; i++) ok[i] -= V[k].sum[i];
    }
    int mid = (l + r) >> 1;
    auto left = bisect_from_left(a, f, ok, 2 * k + 1, l, mid);
    if(left != -1) return left;
    return bisect_from_left(a, f, ok, 2 * k + 2, mid, r);
  }
  template<typename F>
  int bisect_from_right(int a, const F &f, std::array<int, bitlen> &ok, int k, int l, int r){
    if(a <= l) return -1;
    if(k < M - 1) push_down(k, l, r);
    if(r <= a){
      for(int i = 0; i < bitlen; i++) ok[i] += V[k].sum[i];
      if(!f(ok)) return -1;
      if(k >= M - 1) return l;
      for(int i = 0; i < bitlen; i++) ok[i] -= V[k].sum[i];
    }
    int mid = (l + r) >> 1;
    auto right = bisect_from_right(a, f, ok, 2 * k + 2, mid, r);
    if(right != -1) return right;
    return bisect_from_right(a, f, ok, 2 * k + 1, l, mid);
  }
  template<typename F>
  int bisect_from_left2(int a, const F &f, std::array<int, bitlen> &ok, int k, int l, int r){
    if(r <= a) return -1;
    if(k < M - 1) push_down(k, l, r);
    if(a <= l){
      for(int i = 0; i < bitlen; i++) ok[i] += V[k].sum[i];
      if(!f(ok, a, r)) return -1;
      if(k >= M - 1) return l;
      for(int i = 0; i < bitlen; i++) ok[i] -= V[k].sum[i];
    }
    int mid = (l + r) >> 1;
    auto left = bisect_from_left2(a, f, ok, 2 * k + 1, l, mid);
    if(left != -1) return left;
    return bisect_from_left2(a, f, ok, 2 * k + 2, mid, r);
  }
  template<typename F>
  int bisect_from_right2(int a, const F &f, std::array<int, bitlen> &ok, int k, int l, int r){
    if(a <= l) return -1;
    if(k < M - 1) push_down(k, l, r);
    if(r <= a){
      for(int i = 0; i < bitlen; i++) ok[i] += V[k].sum[i];
      if(!f(ok, l, a)) return -1;
      if(k >= M - 1) return l;
      for(int i = 0; i < bitlen; i++) ok[i] -= V[k].sum[i];
    }
    int mid = (l + r) >> 1;
    auto right = bisect_from_right2(a, f, ok, 2 * k + 2, mid, r);
    if(right != -1) return right;
    return bisect_from_right2(a, f, ok, 2 * k + 1, l, mid);
  }
public:
  range_bit_operation(): N(0), M(0){
    static_assert(bitlen > 0);
    static_assert(sizeof(Val) * 8 >= bitlen);
  }
  range_bit_operation(int n): N(n), M(ceil_pow2(N)), V(2 * M - 1, node()){
    static_assert(bitlen > 0);
    static_assert(sizeof(Val) * 8 >= bitlen);
  }
  template<typename T>
  range_bit_operation(const std::vector<T> &v): N(v.size()), M(ceil_pow2(N)), V(2 * M - 1){
    static_assert(bitlen > 0);
    static_assert(sizeof(Val) * 8 >= bitlen);
    for(int i = 0; i < N; i++){
      for(int j = 0; j < bitlen; j++){
        V[i + M - 1].sum[j] = (v[i] >> j) & 1;
      }
    }
    for(int i = M - 2; i >= 0; i--) eval(i);
  }
  int size(){
    return N;
  }
  // val[k] <- x, O(logN + bitlen)
  void set(int k, Val x){
    set(k, x, 0, 0, M);
  }
  // val[k], O(logN + bitlen)
  Val get(int k){
    return get(k, 0, 0, M);
  }
  // O(logN * bitlen)
  void update_or(int l, int r, Val x){
    assert(0 <= l && r <= size());
    update(l, r, Lazy(x, inf, 0), 0, 0, M);
  }
  // O(logN * bitlen)
  void update_and(int l, int r, Val x){
    assert(0 <= l && r <= size());
    update(l, r, Lazy(0, x, 0), 0, 0, M);
  }
  // O(logN * bitlen)
  void update_xor(int l, int r, Val x){
    assert(0 <= l && r <= size());
    update(l, r, Lazy(0, inf, x), 0, 0, M);
  }
  // O(logN * bitlen)
  void update_set(int l, int r, Val x){
    assert(0 <= l && r <= size());
    update(l, r, Lazy(x, x, 0), 0, 0, M);
  }
  // O(logN * bitlen)
  Val query_or(int l, int r){
    return query_or(l, r, 0, 0, M);
  }
  // O(logN * bitlen)
  Val query_and(int l, int r){
    return query_and(l, r, 0, 0, M);
  }
  // O(logN * bitlen)
  Val query_xor(int l, int r){
    return query_xor(l, r, 0, 0, M);
  }
  // O(logN * bitlen)
  ValSum query_sum(int l, int r){
    return query_sum(l, r, 0, 0, M);
  }
  // O(logN * bitlen)
  int query_popcount(int l, int r){
    return query_popcount(l, r, 0, 0, M);
  }
  // O(bitlen)
  void update_or_all(Val x){
    update(0, M, Lazy(x, inf, 0), 0, 0, M);
  }
  // O(bitlen)
  void update_and_all(Val x){
    update(0, M, Lazy(0, x, 0), 0, 0, M);
  }
  // O(bitlen)
  void update_xor_all(Val x){
    update(0, M, Lazy(0, inf, x), 0, 0, M);
  }
  // O(bitlen)
  Val query_or_all(){
    return query_or(0, M, 0, 0, M);
  }
  // O(bitlen)
  Val query_and_all(){
    return query_and(0, M, 0, 0, M);
  }
   // O(bitlen)
  Val query_xor_all(){
    return query_xor(0, M, 0, 0, M);
  }
  // O(bitlen)
  ValSum query_sum_all(){
    return query_sum(0, M, 0, 0, M);
  }
  // O(bitlen)
  int query_popcount_all(){
    return query_popcount(0, M, 0, 0, M);
  }
  // f(sum[l, r])が初めてtrueになる最大のl, 無い場合は-1
  template<typename F>
  int bisect_from_left(int l, const F &f){
    std::array<int, bitlen> ok;
    ok.fill(0);
    return bisect_from_left(l, f, ok, 0, 0, M);
  }
  // 与える関数が f(x, l, r)
  template<typename F>
  int bisect_from_left2(int l, const F &f){
    std::array<int, bitlen> ok;
    ok.fill(0);
    return bisect_from_left2(l, f, ok, 0, 0, M);
  }
  // f(sum[l, r])が初めてtrueになる最小のr, 無い場合は-1
  template<typename F>
  int bisect_from_right(int r, const F &f){
    std::array<int, bitlen> ok;
    ok.fill(0);
    return bisect_from_right(r, f, ok, 0, 0, M);
  }
  template<typename F>
  int bisect_from_right2(int r, const F &f){
    std::array<int, bitlen> ok;
    ok.fill(0);
    return bisect_from_right2(r, f, ok, 0, 0, M);
  }
};
#endif