#ifndef _SLOPE_CHMIN_H_
#define _SLOPE_CHMIN_H_
#include "cht.hpp"
#include "../../data_structure/range_query/sparse_table.hpp"
#include "../../data_structure/segment_tree/binary_indexed_tree.hpp"

// chtのenumerateを使って以下を列挙する
// {l, r, f} := [l, r)で直線or線分fが最小値をとる

// 以下のことができる
// f(x)をxにおける最小値とする
// range_min(l, r) : l <= x < rを満たす整数xの中でf(x)の最小値
// range_max(l, r) : l <= x < rを満たす整数xの中でf(x)の最大値
// range_sum(l, r) : l <= x < rを満たす整数xの中でf(x)の和
template<typename Val, typename ValSum>
struct static_slope_chmin{
private:
  bool is_built = false;
  Val x_low, x_high;
  lichao_tree<Val> cht;
  std::vector<Val> lx;
  std::vector<typename lichao_tree<Val>::line> ln;
  std::vector<Val> sum;

  static constexpr Val merge_min(Val a, Val b){
    return std::min(a, b);
  }
  static constexpr Val min_id(){return std::numeric_limits<Val>::max();}
  static constexpr Val merge_max(Val a, Val b){
    return std::max(a, b);
  }
  static constexpr Val max_id(){return std::numeric_limits<Val>::min();}
  
  // ak + b のl <= k < rにおける和
  // a / 2 * {r * (r - 1) - l * (l - 1)} + b * (r - l)
  static constexpr ValSum line_sum(Val a, Val b, Val l, Val r){
    return a * ((ValSum)r * (r - 1) - (ValSum)l * (l - 1)) / 2 + (ValSum)b * (r - l);
  }
  sparse_table<Val, merge_min, min_id> st_min;
  sparse_table<Val, merge_max, max_id> st_max;
  binary_indexed_tree<ValSum> bit;
public:
  static_slope_chmin(Val x_low, Val x_high): x_low(x_low), x_high(x_high), cht(x_low, x_high){}

  void add_line(Val a, Val b){
    assert(!is_built);
    cht.add_line(a, b);
  }
  void add_segment(Val l, Val r, Val a, Val b){
    assert(!is_built);
    cht.add_segment(l, r, a, b);
  }
  void build(){
    is_built = true;
    auto v = cht.enumerate();
    int n = v.size();
    std::vector<Val> _min(n), _max(n);
    std::vector<ValSum> _sum(n);
    for(int i = 0; i < n; i++){
      lx.push_back(v[i].first);
      ln.push_back(v[i].second);
      Val r = (i == n - 1 ? x_high : lx[i + 1] - 1); // [l, r]
      Val x = v[i].second.get(lx[i]);
      Val y = v[i].second.get(r);
      if(x > y) std::swap(x, y);
      _min[i] = x;
      _max[i] = y;
    }
    for(int i = 0; i < n; i++){
      Val l = lx[i];
      Val r = (i == n - 1 ? x_high + 1 : lx[i + 1]);
      _sum[i] = line_sum(ln[i].a, ln[i].b, l, r);
    }
    st_min = sparse_table<Val, merge_min, min_id>(_min);
    st_max = sparse_table<Val, merge_max, max_id>(_max);
    bit = binary_indexed_tree<ValSum>(_sum);
  }
  Val min(Val x){
    return cht.min(x);
  }
  Val range_min(Val l, Val r){
    assert(is_built);
    int l2 = std::lower_bound(lx.begin(), lx.end(), l) - lx.begin();
    int r2 = std::upper_bound(lx.begin(), lx.end(), r) - lx.begin();
    assert(l2 > 0 && r2 > 0);
    // [l2, r2 - 2]の線を完全に含む
    // l2 - 1, r2 - 1の線を部分的に含む
    Val res = cht.inf;
    if(l2 < r2 - 1) res = st_min.query(l2, r2 - 1);
    if(l2 != r2){
      // [l, lx[l2])
      if(l < lx[l2]) res = std::min(res, std::min(ln[l2 - 1].get(l), ln[l2 - 1].get(lx[l2] - 1)));
      // [lx[r2 - 1], r)
      if(lx[r2 - 1] < r) res = std::min(res, std::min(ln[r2 - 1].get(lx[r2 - 1]), ln[r2 - 1].get(r - 1)));
    }else{
      if(l < r) res = std::min(res, std::min(ln[l2 - 1].get(l), ln[l2 - 1].get(r - 1)));
    }
    return res;
  }
  Val range_max(Val l, Val r){
    assert(is_built);
    int l2 = std::lower_bound(lx.begin(), lx.end(), l) - lx.begin();
    int r2 = std::upper_bound(lx.begin(), lx.end(), r) - lx.begin();
    assert(l2 > 0 && r2 > 0);
    // [l2, r2 - 2]の線を完全に含む
    // l2 - 1, r2 - 1の線を部分的に含む
    Val res = -cht.inf;
    if(l2 < r2 - 1) res = st_max.query(l2, r2 - 1);
    if(l2 != r2){
      // [l, lx[l2])
      if(l < lx[l2]) res = std::max(res, std::max(ln[l2 - 1].get(l), ln[l2 - 1].get(lx[l2] - 1)));
      // [lx[r2 - 1], r)
      if(lx[r2 - 1] < r) res = std::max(res, std::max(ln[r2 - 1].get(lx[r2 - 1]), ln[r2 - 1].get(r - 1)));
    }else{
      if(l < r) res = std::max(res, std::max(ln[l2 - 1].get(l), ln[l2 - 1].get(r - 1)));
    }
    return res;
  }
  // https://atcoder.jp/contests/abc303/submissions/43402759
  ValSum range_sum(Val l, Val r){
    assert(is_built);
    int l2 = std::lower_bound(lx.begin(), lx.end(), l) - lx.begin();
    int r2 = std::upper_bound(lx.begin(), lx.end(), r) - lx.begin();
    assert(l2 > 0 && r2 > 0);
    // [l2, r2 - 2]の線を完全に含む
    // l2 - 1, r2 - 1の線を部分的に含む
    ValSum res = 0;
    if(l2 < r2 - 1) res = bit.query(l2, r2 - 1);
    if(l2 != r2){
      // [l, lx[l2])
      if(l < lx[l2]) res += line_sum(ln[l2 - 1].a, ln[l2 - 1].b, l, lx[l2]);
      // [lx[r2 - 1], r)
      if(lx[r2 - 1] < r) res += line_sum(ln[r2 - 1].a, ln[r2 - 1].b, lx[r2 - 1], r);
    }else{
      if(l < r) res += line_sum(ln[l2 - 1].a, ln[l2 - 1].b, l, r);
    }
    return res;
  }
};
#endif