#ifndef _ABELIAN_GROUP_H_
#define _ABELIAN_GROUP_H_

template<typename T>
struct range_sum_abelian{
  using Val = T;
  static Val id(){return 0;}
  static Val inv(Val a){return -a;}
  static Val merge(Val a, Val b){return a + b;}
};
template<typename T>
struct range_add_range_sum_abelian{
  using Val = T;
  using Lazy = T;
  static Val id(){return 0;}
  static Lazy id_lazy(){return 0;}
  static Val merge(Val a, Val b){return a + b;}
  static Val inv(Val a){return -a;}
  static Lazy propagate(Lazy a, Lazy b){return a + b;}
  static Lazy inv_lazy(Lazy a){return -a;}
  static Val apply(Val a, Lazy b, int l, int r){return a + b * (r - l);}
};
template<typename T>
using components_add_components_sum_abelian = range_add_range_sum_abelian<T>;
#endif