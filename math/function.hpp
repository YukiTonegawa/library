#ifndef _FUNCTION_H_
#define _FUNCTION_H_
#include <algorithm>
#include <limits>
#include <cassert>

template<typename T>
struct add_min_function{
  using F = add_min_function<T>;
  T add, upper;
  static constexpr T inf = std::numeric_limits<T>::max() / 2;
  add_min_function(): add(0), upper(inf){}
  add_min_function(T _add, T _upper): add(_add), upper(_upper){}
  // g(f(x))
  static F merge(const F &f, const F &g){return {f.add + g.add, std::min(g.upper, f.upper + g.add)};}
  T fx(T x){return std::min(x + add, upper);}
  bool operator == (const F &r)const{return add == r.add && upper == r.upper;}
};
template<typename T>
struct clamp_function{
  T add, lower, upper;
  static constexpr T inf = std::numeric_limits<T>::max() / 2;
  static constexpr T minf = std::numeric_limits<T>::min() / 2;
  clamp_function(): add(0), lower(minf), upper(inf){}
  clamp_function(T _add, T _lower, T _upper): add(_add), lower(_lower), upper(_upper){lower = std::min(_lower, _upper);}
  // g(f(x))
  static clamp_function<T> merge(const clamp_function<T> &f, const clamp_function<T> &g){
    return {f.add + g.add, std::max(g.lower, std::min(f.upper, f.lower) + g.add), std::min(g.upper, std::max(g.lower, f.upper + g.add))};
  }
  T fx(T x){
    return std::min(std::max(x + add, lower), upper);
  }
  bool operator == (const clamp_function<T> &r)const{return add == r.add && lower == r.lower && upper == r.upper;}
};
// min部分によって減少した値をスコアとする
template<typename T, typename Tsum>
struct clamp_function_score{
public:
  T add, lower, upper, score_upper;
  Tsum score_sum;
  static constexpr T inf = std::numeric_limits<T>::max() / 2;
  static constexpr T minf = std::numeric_limits<T>::min() / 2;
  clamp_function_score(): add(0), lower(minf), upper(inf), score_upper(inf), score_sum(0){}
  clamp_function_score(T _add, T _lower, T _upper): add(_add), lower(_lower), upper(_upper), score_upper(_upper), score_sum(0){lower = std::min(_lower, _upper);}
private:
  clamp_function_score(T _add, T _lower, T _upper, T _supper, Tsum _ssum): add(_add), lower(_lower), upper(_upper), score_upper(_supper), score_sum(_ssum){}
public:
  // g(f(x))
  static clamp_function_score<T, Tsum> merge(const clamp_function_score<T, Tsum> &f, const clamp_function_score<T, Tsum> &g){
    T _add = f.add + g.add;
    T _lower = std::max(g.lower, std::min(f.upper, f.lower) + g.add);
    T _upper = std::min(g.upper, std::max(g.lower, f.upper + g.add));
    _lower = std::min(_lower, _upper);
    Tsum _ssum = f.score_sum + g.score_sum;
    if(f.lower == f.upper){
      _ssum += std::max((T)0, f.lower + g.add - g.score_upper);
      return {_add, _lower, _upper, f.score_upper + g.add, _ssum};
    }else if(f.lower + g.add >= g.score_upper){
      assert(_lower == _upper);
      _ssum += std::max((T)0, f.lower + g.add - g.score_upper);
      return {_add, _lower, _upper, f.lower + g.add, _ssum};
    }else{
      return {_add, _lower, _upper, std::min(f.score_upper + g.add, g.score_upper), _ssum};
    }
  }
  T fx(T x){
    return std::min(std::max(x + add, lower), upper);
  }
  Tsum fx_score(T x){
    return score_sum + std::max((T)0, x + add - score_upper);
  }
  bool operator == (const clamp_function_score<T, Tsum> &r)const{return add == r.add && lower == r.lower && upper == r.upper && score_upper == r.score_upper && score_sum == r.score_sum;}
};
template<typename T>
struct prefixsum_min{
  T sum, pmin;
  static prefixsum_min<T> id(){
    return {0, std::numeric_limits<T>::max()};
  }
  static prefixsum_min<T> merge(prefixsum_min<T> a, prefixsum_min<T> b){
    if(b.pmin == std::numeric_limits<T>::max()) return a;
    return {a.sum + b.sum, std::max(a.pmin, a.sum + b.pmin)};
  }
};
template<typename T>
struct prefixsum_max{
  T sum, pmax;
  static prefixsum_max<T> id(){
    return {0, std::numeric_limits<T>::min()};
  }
  static prefixsum_max<T> merge(prefixsum_max<T> a, prefixsum_max<T> b){
    if(b.pmax == std::numeric_limits<T>::min()) return a;
    return {a.sum + b.sum, std::max(a.pmax, a.sum + b.pmax)};
  }
};
template<typename T>
struct substringsum_max{
  T sum, pmax, smax, ssmax; // {区間全体のsum, prefixsumのmax, suffixsumのmax, substringsumのmax}
  static substringsum_max<T> id(){
    return {0, std::numeric_limits<T>::min(), 0, 0};
  }
  static substringsum_max<T> merge(substringsum_max<T> a, substringsum_max<T> b){
    if(b.pmax == std::numeric_limits<T>::min()) return a;
    if(a.pmax == std::numeric_limits<T>::min()) return b;
    T sum = a.sum + b.sum;
    T pmax = std::max(a.pmax, a.sum + b.pmax);
    T smax = std::max(a.smax + b.sum, b.smax);
    return {sum, pmax, smax, std::max({a.ssmax, b.ssmax, pmax, smax})};
  }
};
// 01列の0を-1として扱う
// 1の数, 合計, 接頭辞のmin, 接頭辞のmax
struct excess_value{
public:
  int rank, sum, pmin, pmax;
  static constexpr int inf = 1 << 30;
  excess_value(): rank(inf), pmin(inf), pmax(-inf){}
  excess_value(bool f): rank(f), sum(f ? 1 : -1), pmin(sum), pmax(sum){}
private:
  excess_value(int a, int b, int c, int d): rank(a), sum(b), pmin(c), pmax(d){}
public:
  static excess_value merge(excess_value a, excess_value b){
    if(a.rank == inf) return b;
    if(b.rank == inf) return a;
    return {a.rank + b.rank, a.sum + b.sum, std::min(a.pmin, a.sum + b.pmin), std::max(a.pmax, a.sum + b.pmax)};
  }
};
// 01列の0を-1として扱う
// 1の数, 合計, 接頭辞のmin, 接頭辞のmax, 接尾辞のmin, 接尾辞のmax
struct excess_value2{
public:
  int rank, sum, pmin, pmax, smin, smax;
  static constexpr int inf = 1 << 30;
  excess_value2(): rank(inf), pmin(inf), pmax(-inf), smin(inf), smax(-inf){}
  excess_value2(bool f): rank(f), sum(f ? 1 : -1), pmin(sum), pmax(sum), smin(sum), smax(sum){}
private:
  excess_value2(int a, int b, int c, int d, int e, int f): rank(a), sum(b), pmin(c), pmax(d), smin(e), smax(f){}
public:
  static excess_value2 merge(excess_value2 a, excess_value2 b){
    if(a.rank == inf) return b;
    if(b.rank == inf) return a;
    return {a.rank + b.rank, a.sum + b.sum, std::min(a.pmin, a.sum + b.pmin), 
    std::max(a.pmax, a.sum + b.pmax), std::min(a.smin + b.sum, b.smin), std::max(a.smax + b.sum, b.smax)};
  }
};

#endif