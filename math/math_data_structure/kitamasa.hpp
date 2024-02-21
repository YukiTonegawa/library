#ifndef _KITAMASA_H_
#define _KITAMASA_H_
#include <vector>
#include <algorithm>
#include "../../algebraic_structure/semi_ring.hpp"
// A[n] -> A[n + 1], O(s)
template<typename __func_st, typename Val>
std::vector<Val> _increment(const std::vector<Val> &x, const std::vector<Val> &c){
  int s = c.size();
  std::vector<Val> ret(s, __func_st::template add_id<Val>());
  for(int i = 0; i < s; i++){
    if(i) ret[i] = x[i - 1];
    ret[i] = __func_st::template add<Val>(ret[i], __func_st::template mul<Val>(x[s - 1], c[s - 1 - i]));
  }
  return ret;
}
// A[n] -> A[2n], O(s^2)
template<typename __func_st, typename Val>
std::vector<Val> _shift(const std::vector<Val> &x, const std::vector<Val> &c){
  int s = c.size();
  std::vector<Val> A = x;
  std::vector<Val> ret(s, __func_st::template add_id<Val>());
  for(int i = 0; i < s; i++){
    for(int j = 0; j < s; j++) ret[j] = __func_st::template add<Val>(ret[j], __func_st::template mul<Val>(x[i], A[j]));
    if(i != s - 1) A = _increment<__func_st, Val>(A, c);
  }
  return ret;
}
// O(s^2 log(k))
// A[i] = A[i - 1] * c[0] + A[i - 2] + c[1]...を満たす漸化式のk項目
template<typename __func_st, typename Val>
std::vector<Val> kitamasa(long long k, const std::vector<Val> &c){
  int s = c.size();
  std::vector<Val> ret(s, __func_st::template add_id<Val>());
  if(k < s){
    ret[k] = __func_st::template mul_id<Val>();
    return ret;
  }
  ret[0] = __func_st::template mul_id<Val>();
  for(int l = 63 - __builtin_clzll(k); l >= 0; l--){
    ret = _shift<__func_st, Val>(ret, c);
    if((k >> l) & 1) ret = _increment<__func_st, Val>(ret, c);
  }
  return ret;
}
#endif