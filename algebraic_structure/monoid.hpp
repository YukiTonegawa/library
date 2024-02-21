#ifndef _MONOID_H_
#define _MONOID_H_
#include <limits>
#include <array>
#include <numeric>
#include <algorithm>
#include "../math/function.hpp"

template<typename T>
struct range_sum{
  using Val = T;
  static Val id(){return 0;}
  static Val merge(Val a, Val b){return a + b;}
  static Val flip(Val a){return a;}
};
template<typename T>
struct range_min{
  using Val = T;
  static Val id(){ return std::numeric_limits<Val>::max();}
  static Val merge(Val a, Val b){ return std::min(a, b);}
};
template<typename T>
struct range_max{
  using Val = T;
  static Val id(){return std::numeric_limits<Val>::min();}
  static Val merge(Val a, Val b){return std::max(a, b);}
};
template<typename T>
struct range_mul{
  using Val = T;
  static Val id(){return 1;}
  static Val merge(Val a, Val b){return a * b;}
};
template<typename T>
struct range_second_min{
  using Val = T;
  static Val id(){return {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};}
  static Val merge(Val a, Val b){
    // 最小値が2つある場合{最小値, 最小値}
    if(a.first <= b.first) return {a.first, std::min(a.second, b.first)};
    else return {b.first, std::min(a.first, b.second)};
    // 最小値が2つある場合{最小値, 2番目}
    /*
    if(a.first == b.first) return {a.first, std::min(a.second, b.second)};
    else if(a.first < b.first) return {a.first, std::min(a.second, b.first)};
    else return {b.first, std::min(a.first, b.second)};
    */
  }
};
template<typename T>
struct range_second_max{
  using Val = T;
  static Val id(){return {std::numeric_limits<T>::min(), std::numeric_limits<T>::min()};}
  static Val merge(Val a, Val b){
    // 最大値が2つある場合{最大値, 最大値}
    if(a.first <= b.first) return {a.first, std::max(a.second, b.first)};
    else return {b.first, std::max(a.first, b.second)};
    // 最大値が2つある場合{最大値, 2番目}
    /*
    if(a.first == b.first) return {a.first, std::max(a.second, b.second)};
    else if(a.first > b.first) return {a.first, std::max(a.second, b.first)};
    else return {b.first, std::max(a.first, b.second)};
    */
  }
};
template<typename mint>
struct range_composite{
  using Val = std::pair<mint, mint>;
  static Val id(){return {1, 0};}
  static Val merge(Val a, Val b){return {a.first * b.first, a.second * b.first + b.second};}
};
template<typename mint>
struct range_affine{
  using Val = std::pair<mint, mint>;
  static Val id(){return {1, 0};}
  static Val merge(Val a, Val b){return {a.first * b.first, a.second * b.first + b.second};}
};
template<typename mint>
struct range_composite_flip{
  using Val = std::tuple<mint, mint, mint>;
  static Val id(){return {1, 0, 0};}  
  static Val flip(Val a){
    auto [x, y, z] = a;
    return {x, z, y};
  }
  static Val merge(Val a, Val b){
    auto [x, y, z] = a;
    auto [X, Y, Z] = b;
    return {x * X, y * X + Y, Z * x + z};
  }
};
// f(x) = min(max(x + a, b), c)
template<typename T>
struct range_clamp{
  using Val = clamp_function<T>;
  static Val id(){return Val();}
  // g(f(x))
  static Val merge(Val f, Val g){return Val::merge(f, g);}
};
// f(x) = min(max(x + a, b), c), minで減少した値をスコアとする
template<typename T, typename Tsum>
struct range_clamp_score{
  using Val = clamp_function_score<T, Tsum>;
  static Val id(){return Val();}
  static Val merge(Val a, Val b){return Val::merge(a, b);}
};
template<typename T>
struct range_gcd{
  using Val = T;
  static Val __gcd(Val a, Val b){
    if(!a || !b) return !a ? b : a;
    if(a < b) std::swap(a, b);
    while(b){
      a %= b;
      std::swap(a, b);
    }
    return a;
  }
  static Val id(){return 0;}
  static Val merge(Val a, Val b){return __gcd(a, b);}
};
struct range_excess_value{
  using Val = excess_value;
  static Val id(){return Val();}
  static Val merge(Val a, Val b){return Val::merge(a, b);}
};
template<typename T>
struct range_prefixsum_min{
  using Val = prefixsum_min<T>;
  static Val id(){return Val::id();}
  static Val merge(Val a, Val b){return Val::merge(a, b);}
};
template<typename T>
struct range_prefixsum_max{
  using Val = prefixsum_max<T>;
  static Val id(){return Val::id();}
  static Val merge(Val a, Val b){return Val::merge(a, b);}
};
template<typename T>
struct range_substringsum_max{
  using Val = substringsum_max<T>;
  static Val id(){return Val::id();}
  static Val merge(Val a, Val b){return Val::merge(a, b);}
};
// 区間set, これまでにsetした物の中ならどれかを取得
template<typename T>
struct range_pointer_set_get_any{
  using Val = T*;
  using Lazy = T*;
  static Val id(){return nullptr;}
  static Lazy id_lazy(){return nullptr;}
  static Lazy propagate(Lazy a, Lazy b){return b ? b : a;}
  static Val apply(Val a, Lazy b, int l, int r){return b ? b : a;}
};
template<typename T>
struct range_add_range_sum{
  using Val = T;
  using Lazy = T;
  static Val id(){return 0;}
  static Lazy id_lazy(){return 0;}
  static Val merge(Val a, Val b){return a + b;}
  static Val apply(Val a, Lazy b, int l, int r){return a + (Val)b * (r - l);}
  static Lazy propagate(Lazy a, Lazy b){return a + b;}
  static Val flip(Val a){return a;}
};
template<typename T>
struct range_min_range_min{
  using Val = T;
  using Lazy = T;
  static Val id(){return std::numeric_limits<Val>::max();}
  static Lazy id_lazy(){return std::numeric_limits<Lazy>::max();}
  static Val merge(Val a, Val b){return std::min(a, b);}
  static Val apply(Val a, Lazy b, int l, int r){return std::min(a, b);}
  static Lazy propagate(Lazy a, Lazy b){return std::min(a, b);}
  static Val flip(Val a){return a;}
};
template<typename T>
struct range_max_range_max{
  using Val = T;
  using Lazy = T;
  static Val id(){return std::numeric_limits<Val>::min();}
  static Lazy id_lazy(){return std::numeric_limits<Lazy>::min();}
  static Val merge(Val a, Val b){return std::max(a, b);}
  static Val apply(Val a, Lazy b, int l, int r){return std::max(a, b);}
  static Lazy propagate(Lazy a, Lazy b){return std::max(a, b);}
  static Val flip(Val a){return a;}
};
template<typename T>
struct range_set_range_min{
  using Val = T;
  using Lazy = T;
  static Val id(){return std::numeric_limits<Val>::max();}
  static Lazy id_lazy(){return std::numeric_limits<Lazy>::max();}
  static Val merge(Val a, Val b){return std::min(a, b);}
  static Val apply(Val a, Lazy b, int l, int r){return b;}
  static Lazy propagate(Lazy a, Lazy b){return b;}
  static Val flip(Val a){return a;}
};
template<typename T>
struct range_set_range_max{
  using Val = T;
  using Lazy = T;
  static Val id(){return std::numeric_limits<Val>::min();}
  static Lazy id_lazy(){return std::numeric_limits<Lazy>::min();}
  static Val merge(Val a, Val b){return std::max(a, b);}
  static Val apply(Val a, Lazy b, int l, int r){return b;}
  static Lazy propagate(Lazy a, Lazy b){return b;}
  static Val flip(Val a){return a;}
};
template<typename T>
struct range_set_range_sum{
  using Val = T;
  using Lazy = T;
  static Val id(){return std::numeric_limits<Val>::min();}
  static Lazy id_lazy(){return std::numeric_limits<Lazy>::min();}
  static Val merge(Val a, Val b){return a + b;}
  static Val apply(Val a, Lazy b, int l, int r){return b * (r - l);}
  static Lazy propagate(Lazy a, Lazy b){return b;}
  static Val flip(Val a){return a;}
};
template<typename T>
struct range_add_range_min{
  using Val = T;
  using Lazy = T;
  static Val id(){return std::numeric_limits<T>::max();}
  static Lazy id_lazy(){return 0;}
  static Val merge(Val a, Val b){return std::min(a, b);}
  static Val apply(Val a, Lazy b, int l, int r){
    if(a == id()) return b;
    return a + b;
  }
  static Lazy propagate(Lazy a, Lazy b){return a + b;}
};
template<typename T>
struct range_add_range_max{
  using Val = T;
  using Lazy = T;
  static Val id(){return std::numeric_limits<T>::min();}
  static Lazy id_lazy(){return 0;}
  static Val merge(Val a, Val b){return std::max(a, b);}
  static Val apply(Val a, Lazy b, int l, int r){
    if(a == id()) return b;
    return a + b;
  }
  static Lazy propagate(Lazy a, Lazy b){return a + b;}
};
template<typename T>
struct range_add_range_argmin{
  using Val = std::pair<T, int>;
  using Lazy = T;
  static Val id(){return {std::numeric_limits<T>::max(), -1};}
  static Lazy id_lazy(){return 0;}
  static Val merge(Val a, Val b){return std::min(a, b);}
  static Val apply(Val a, Lazy b, int l, int r){
    if(a == id()) return a;
    return {a.first + b, a.second};
  }
  static Lazy propagate(Lazy a, Lazy b){return a + b;}
};

template<typename T>
struct range_add_range_mincount{
  using Val = std::pair<T, int>;
  using Lazy = T;
  static Val id(){return {std::numeric_limits<T>::max(), 0};}
  static Lazy id_lazy(){return 0;}
  static Val merge(Val a, Val b){
    if(a.first != b.first) return a.first < b.first ? a : b;
    return {a.first, a.second + b.second};
  }
  static Val apply(Val a, Lazy b, int l, int r){
    if(a == id()) return a;
    return {a.first + b, a.second};
  }
  static Lazy propagate(Lazy a, Lazy b){return a + b;}
};
template<typename mint>
struct range_affine_range_sum{
  using Val = mint;
  using Lazy = std::pair<mint, mint>;
  static Val id(){return 0;}
  static Lazy id_lazy(){return {1, 0};}
  static Val merge(Val a, Val b){return a + b;}
  static Val apply(Val a, Lazy b, int l, int r){return a * b.first + b.second * (r - l);}
  static Lazy propagate(Lazy a, Lazy b){return {a.first * b.first, a.second * b.first + b.second};}
};
#endif