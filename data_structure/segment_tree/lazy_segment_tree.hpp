#ifndef _LAZY_SEGMENT_TREE_H_
#define _LAZY_SEGMENT_TREE_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include <array>
#include <queue>
#include "../../algebraic_structure/monoid.hpp"

template<typename monoid>
struct lazy_segment_tree{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
private:
  int N, M;
  struct node{
    Val sum;
    Lazy lazy;
    node(Val val = id()): sum(val), lazy(id_lazy()){}
  };
  static constexpr int ceil_pow2(int n){
    int m = 1;
    while(m < n) m <<= 1;
    return m;
  }
  std::vector<node> V;
  inline void push_down(int k, int l, int r){
    if(V[k].lazy == id_lazy()) return;
    int mid = (l + r) >> 1;
    propagate(k * 2 + 1, l, mid, V[k].lazy);
    propagate(k * 2 + 2, mid, r, V[k].lazy);
    V[k].lazy = id_lazy();
  }
  inline void propagate(int k, int l, int r, Lazy x){
    V[k].sum = apply(V[k].sum, x, l, r);
    V[k].lazy = propagate_lazy(V[k].lazy, x);
  }
  Val set(int a, Val x, int k, int l, int r){
    if(r - l == 1) return V[k].sum = x;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid){
      return V[k].sum = merge(set(a, x, k * 2 + 1, l, mid), V[k * 2 + 2].sum);
    }else{
      return V[k].sum = merge(V[k * 2 + 1].sum, set(a, x, k * 2 + 2, mid, r));
    }
  }
  Val get(int a, int k, int l, int r){
    if(r - l == 1) return V[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid) return get(a, k * 2 + 1, l, mid);
    else return get(a, k * 2 + 2, mid, r);
  }
  Val update(int a, int b, Lazy x, int k, int l, int r){
    if(r <= a || b <= l) return V[k].sum;
    if(a <= l && r <= b){
      propagate(k, l, r, x);
      return V[k].sum;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return V[k].sum = merge(update(a, b, x, k * 2 + 1, l, mid), update(a, b, x, k * 2 + 2, mid, r));
  }
  Val query(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return id();
    if(a <= l && r <= b) return V[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return merge(query(a, b, k * 2 + 1, l, mid), query(a, b, k * 2 + 2, mid, r));
  }
  template<typename F>
  inline std::pair<int, Val> bisect_from_left(int a, const F &f, Val val, int k, int l, int r){
    if(r <= a) return {-1, val};
    if(k < M - 1) push_down(k, l, r);
    if(a <= l){
      Val merged = merge(val, V[k].sum);
      if(!f(merged)) return {-1, merged};
      if(k >= M - 1) return {l, merged};
    }
    int mid = (l + r) >> 1;
    auto left = bisect_from_left(a, f, val, 2 * k + 1, l, mid);
    if(left.first != -1) return left;
    return bisect_from_left(a, f, left.second, 2 * k + 2, mid, r);
  }
  template<typename F>
  inline std::pair<int, Val> bisect_from_left2(int a, const F &f, Val val, int k, int l, int r){
    if(r <= a) return {-1, val};
    if(k < M - 1) push_down(k, l, r);
    if(a <= l){
      Val merged = merge(val, V[k].sum);
      if(!f(merged, l, r)) return {-1, merged};
      if(k >= M - 1) return {l, merged};
    }
    int mid = (l + r) >> 1;
    auto left = bisect_from_left2(a, f, val, 2 * k + 1, l, mid);
    if(left.first != -1) return left;
    return bisect_from_left2(a, f, left.second, 2 * k + 2, mid, r);
  }
  template<typename F>
  inline std::pair<int, Val> bisect_from_right(int a, const F &f, Val val, int k, int l, int r){
    if(a <= l) return {-1, val};
    if(k < M - 1) push_down(k, l, r);
    if(r <= a){
      Val merged = merge(val, V[k].sum);
      if(!f(merged)) return {-1, merged};
      if(k >= M - 1) return {l, merged};
    }
    int mid = (l + r) >> 1;
    auto right = bisect_from_right(a, f, val, 2 * k + 2, mid, r);
    if(right.first != -1) return right;
    return bisect_from_right(a, f, right.second, 2 * k + 1, l, mid);
  }
  template<typename F>
  inline std::pair<int, Val> bisect_from_right2(int a, const F &f, Val val, int k, int l, int r){
    if(a <= l) return {-1, val};
    if(k < M - 1) push_down(k, l, r);
    if(r <= a){
      Val merged = merge(val, V[k].sum);
      if(!f(merged, l, r)) return {-1, merged};
      if(k >= M - 1) return {l, merged};
    }
    int mid = (l + r) >> 1;
    auto right = bisect_from_right2(a, f, val, 2 * k + 2, mid, r);
    if(right.first != -1) return right;
    return bisect_from_right2(a, f, right.second, 2 * k + 1, l, mid);
  }
public:
  lazy_segment_tree(): N(0), M(0){}
  lazy_segment_tree(int n): N(n), M(ceil_pow2(N)), V(2 * M - 1, node()){}
  lazy_segment_tree(const std::vector<Val> &v): N(v.size()), M(ceil_pow2(N)), V(2 * M - 1){
    for(int i = 0; i < M; i++) V[i + M - 1] = node(i < N ? v[i] : id());
    for(int i = M - 2; i >= 0; i--) V[i].sum = merge(V[i * 2 + 1].sum, V[i * 2 + 2].sum);
  }
  int size(){return N;}
  // val[k] <- x
  Val set(int k, Val x){return set(k, x, 0, 0, M);}
  // val[k]
  Val get(int k){return get(k, 0, 0, M);}
  // sum[a, b)
  Val query(int a, int b){return query(a, b, 0, 0, M);}
  Val query_all(){return V[0].sum;}
  // apply([a, b), x)
  Val update(int a, int b, Lazy x){return update(a, b, x, 0, 0, M);}
  // f(sum[l, r))が初めてtrueになる
  // f(sum[l, i)), l <= i < r    = false
  // f(sum[l, j)), r <= j <= n   = true
  // となるような{r, sum[l, r)} 無い場合は r = -1
  template<typename F>
  std::pair<int, Val> bisect_from_left(int l, const F &f){
    return bisect_from_left(l, f, id(), 0, 0, M);
  }
  template<typename F>
  std::pair<int, Val> bisect_from_left2(int l, const F &f){
    return bisect_from_left2(l, f, id(), 0, 0, M);
  }
  // f(sum[l, r))が初めてtrueになる
  // f(sum[i, r)), l < i < r    = false
  // f(sum[j, r)), 0 <= j <= l  = true
  // となるような{l, sum[l, r)} 無い場合は l = -1
  template<typename F>
  std::pair<int, Val> bisect_from_right(int r, const F &f){
    return bisect_from_right(r, f, id(), 0, 0, M);
  }
  template<typename F>
  std::pair<int, Val> bisect_from_right2(int r, const F &f){
    return bisect_from_right2(r, f, id(), 0, 0, M);
  }
  std::vector<Val> to_list(){
    std::queue<std::tuple<int, int, int>> q;
    q.push({0, 0, M});
    std::vector<Val> res(N, id());
    while(!q.empty()){
      auto [k, l, r] = q.front();
      q.pop();
      if(r - l == 1){
        if(l < N) res[l] = V[k].sum;
      }else{
        push_down(k, l, r);
        int m = (l + r) / 2;
        q.push({k * 2 + 1, l, m});
        q.push({k * 2 + 2, m, r});
      }
    }
    return res;
  }
};

template<typename monoid>
struct lazy_segment_tree_reset{
  using Val = typename monoid::Val;
  using Lazy = typename monoid::Lazy;
  static constexpr auto id = monoid::id;
  static constexpr auto id_lazy = monoid::id_lazy;
  static constexpr auto merge = monoid::merge;
  static constexpr auto apply = monoid::apply;
  static constexpr auto propagate_lazy = monoid::propagate;
private:
  int N, M;
  struct node{
    Val sum;
    Lazy lazy;
    bool reset;
    node(Val val = id()): sum(val), lazy(id_lazy()), reset(false){}
  };
  static constexpr int ceil_pow2(int n){
    int m = 1;
    while(m < n) m <<= 1;
    return m;
  }
  std::vector<node> V;
  inline void push_down(int k, int l, int r){
    if(V[k].reset){
      int mid = (l + r) >> 1;
      __reset(k * 2 + 1, l, mid);
      __reset(k * 2 + 2, mid, r);
      V[k].reset = false;
    }
    if(V[k].lazy != id_lazy()){
      int mid = (l + r) >> 1;
      propagate(k * 2 + 1, l, mid, V[k].lazy);
      propagate(k * 2 + 2, mid, r, V[k].lazy);
      V[k].lazy = id_lazy();
    }
  }
  inline void propagate(int k, int l, int r, Lazy x){
    V[k].sum = apply(V[k].sum, x, l, r);
    V[k].lazy = propagate_lazy(V[k].lazy, x);
  }
  void __reset(int k, int l, int r){
    V[k].sum = id();
    V[k].lazy = id_lazy();
    V[k].reset = true;
  }
  Val set(int a, Val x, int k, int l, int r){
    if(r - l == 1) return V[k].sum = x;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid){
      return V[k].sum = merge(set(a, x, k * 2 + 1, l, mid), V[k * 2 + 2].sum);
    }else{
      return V[k].sum = merge(V[k * 2 + 1].sum, set(a, x, k * 2 + 2, mid, r));
    }
  }
  Val get(int a, int k, int l, int r){
    if(r - l == 1) return V[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    if(a < mid) return get(a, k * 2 + 1, l, mid);
    else return get(a, k * 2 + 2, mid, r);
  }
  Val update(int a, int b, Lazy x, int k, int l, int r){
    if(r <= a || b <= l) return V[k].sum;
    if(a <= l && r <= b){
      propagate(k, l, r, x);
      return V[k].sum;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return V[k].sum = merge(update(a, b, x, k * 2 + 1, l, mid), update(a, b, x, k * 2 + 2, mid, r));
  }
  void reset(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return;
    if(a <= l && r <= b){
      __reset(k, l, r);
      return;
    }
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    reset(a, b, k * 2 + 1, l, mid);
    reset(a, b, k * 2 + 2, mid, r);
    V[k].sum = merge(V[k * 2 + 1].sum, V[k * 2 + 2].sum);
  }
  Val query(int a, int b, int k, int l, int r){
    if(r <= a || b <= l) return id();
    if(a <= l && r <= b) return V[k].sum;
    push_down(k, l, r);
    int mid = (l + r) >> 1;
    return merge(query(a, b, k * 2 + 1, l, mid), query(a, b, k * 2 + 2, mid, r));
  }
  template<typename F>
  inline std::pair<int, Val> bisect_from_left(int a, const F &f, Val val, int k, int l, int r){
    if(r <= a) return {-1, val};
    if(k < M - 1) push_down(k, l, r);
    if(a <= l){
      Val merged = merge(val, V[k].sum);
      if(!f(merged)) return {-1, merged};
      if(k >= M - 1) return {l, merged};
    }
    int mid = (l + r) >> 1;
    auto left = bisect_from_left(a, f, val, 2 * k + 1, l, mid);
    if(left.first != -1) return left;
    return bisect_from_left(a, f, left.second, 2 * k + 2, mid, r);
  }
  template<typename F>
  inline std::pair<int, Val> bisect_from_left2(int a, const F &f, Val val, int k, int l, int r){
    if(r <= a) return {-1, val};
    if(k < M - 1) push_down(k, l, r);
    if(a <= l){
      Val merged = merge(val, V[k].sum);
      if(!f(merged, l, r)) return {-1, merged};
      if(k >= M - 1) return {l, merged};
    }
    int mid = (l + r) >> 1;
    auto left = bisect_from_left2(a, f, val, 2 * k + 1, l, mid);
    if(left.first != -1) return left;
    return bisect_from_left2(a, f, left.second, 2 * k + 2, mid, r);
  }
  template<typename F>
  inline std::pair<int, Val> bisect_from_right(int a, const F &f, Val val, int k, int l, int r){
    if(a <= l) return {-1, val};
    if(k < M - 1) push_down(k, l, r);
    if(r <= a){
      Val merged = merge(val, V[k].sum);
      if(!f(merged)) return {-1, merged};
      if(k >= M - 1) return {l, merged};
    }
    int mid = (l + r) >> 1;
    auto right = bisect_from_right(a, f, val, 2 * k + 2, mid, r);
    if(right.first != -1) return right;
    return bisect_from_right(a, f, right.second, 2 * k + 1, l, mid);
  }
  template<typename F>
  inline std::pair<int, Val> bisect_from_right2(int a, const F &f, Val val, int k, int l, int r){
    if(a <= l) return {-1, val};
    if(k < M - 1) push_down(k, l, r);
    if(r <= a){
      Val merged = merge(val, V[k].sum);
      if(!f(merged, l, r)) return {-1, merged};
      if(k >= M - 1) return {l, merged};
    }
    int mid = (l + r) >> 1;
    auto right = bisect_from_right2(a, f, val, 2 * k + 2, mid, r);
    if(right.first != -1) return right;
    return bisect_from_right2(a, f, right.second, 2 * k + 1, l, mid);
  }
public:
  lazy_segment_tree_reset(): N(0), M(0){}
  lazy_segment_tree_reset(int n): N(n), M(ceil_pow2(N)), V(2 * M - 1, node()){}
  lazy_segment_tree_reset(const std::vector<Val> &v): N(v.size()), M(ceil_pow2(N)), V(2 * M - 1){
    for(int i = 0; i < M; i++) V[i + M - 1] = node(i < N ? v[i] : id());
    for(int i = M - 2; i >= 0; i--) V[i].sum = merge(V[i * 2 + 1].sum, V[i * 2 + 2].sum);
  }
  int size(){return N;}
  // val[k] <- x
  Val set(int k, Val x){return set(k, x, 0, 0, M);}
  // val[k]
  Val get(int k){return get(k, 0, 0, M);}
  // sum[a, b)
  Val query(int a, int b){return query(a, b, 0, 0, M);}
  Val query_all(){return V[0].sum;}
  // apply([a, b), x)
  Val update(int a, int b, Lazy x){return update(a, b, x, 0, 0, M);}
  // [l, r)を単位元に戻す
  void reset(int a, int b){reset(a, b, 0, 0, M);}
  // f(sum[l, r))が初めてtrueになる
  // f(sum[l, i)), l <= i < r    = false
  // f(sum[l, j)), r <= j <= n   = true
  // となるような{r, sum[l, r)} 無い場合は r = -1
  template<typename F>
  std::pair<int, Val> bisect_from_left(int l, const F &f){
    return bisect_from_left(l, f, id(), 0, 0, M);
  }
  template<typename F>
  std::pair<int, Val> bisect_from_left2(int l, const F &f){
    return bisect_from_left2(l, f, id(), 0, 0, M);
  }
  // f(sum[l, r))が初めてtrueになる
  // f(sum[i, r)), l < i < r    = false
  // f(sum[j, r)), 0 <= j <= l  = true
  // となるような{l, sum[l, r)} 無い場合は l = -1
  template<typename F>
  std::pair<int, Val> bisect_from_right(int r, const F &f){
    return bisect_from_right(r, f, id(), 0, 0, M);
  }
  template<typename F>
  std::pair<int, Val> bisect_from_right2(int r, const F &f){
    return bisect_from_right2(r, f, id(), 0, 0, M);
  }
  std::vector<Val> to_list(){
    std::queue<std::tuple<int, int, int>> q;
    q.push({0, 0, M});
    std::vector<Val> res(N, id());
    while(!q.empty()){
      auto [k, l, r] = q.front();
      q.pop();
      if(r - l == 1){
        if(l < N) res[l] = V[k].sum;
      }else{
        push_down(k, l, r);
        int m = (l + r) / 2;
        q.push({k * 2 + 1, l, m});
        q.push({k * 2 + 2, m, r});
      }
    }
    return res;
  }
};
#endif