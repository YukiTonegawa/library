#ifndef _SEMI_RING_H_
#define _SEMI_RING_H_
#include <limits>
#include <algorithm>
template<typename T>
struct add_mul{
  using Val = T;
  static Val add(Val a, Val b){return a + b;}
  static Val mul(Val a, Val b){return a * b;}
  static Val id_add(){return 0;}
  static Val id_mul(){return 1;}
};
template<typename T>
struct min_plus{
  using Val = T;
  static Val add(Val a, Val b){return std::min(a, b);}
  static Val mul(Val a, Val b){return a + b;}
  static Val id_add(){return std::numeric_limits<Val>::max() / 2;}
  static Val id_mul(){return 0;}
};
template<typename T>
struct max_plus{
  using Val = T;
  static Val add(Val a, Val b){return std::max(a, b);}
  static Val mul(Val a, Val b){return a + b;}
  static Val id_add(){return std::numeric_limits<Val>::min() / 2;}
  static Val id_mul(){return 0;}
};
// f(w) = max(w + x, y)というような関数を持たせる
template<typename T>
struct hawker_on_graph{
  using Val = std::pair<T, T>; // {sum, prefixsumのmax}
  static constexpr T minf = std::numeric_limits<T>::min() / 2;
  static Val add(Val a, Val b){
    return {std::max(a.first, b.first), std::max(a.second, b.second)};
  }
  static Val mul(Val a, Val b){
    if(a.first == minf || b.first == minf) return {minf, minf}; // 通行止め
    return {a.first + b.first, std::max(a.second + b.first, b.second)};
  }
  static Val id_add(){return {minf, minf};}
  static Val id_mul(){return {0, minf};}
};
template<typename T>
struct xor_and{
  using Val = T;
  static Val add(Val a, Val b){return a ^ b;}
  static Val mul(Val a, Val b){return a & b;}
  static Val id_add(){return 0;}
  static Val id_mul(){return std::numeric_limits<Val>::max();}
};
#endif