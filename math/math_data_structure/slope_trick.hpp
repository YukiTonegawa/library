#ifndef _SLOPE_TRICK_H_
#define _SLOPE_TRICK_H_
#include "../../minior/multiset_sum.hpp"
#include <queue>
#include <vector>
#include <limits>

template<typename Val>
struct slope_trick{
private:
  Val min_fx;
  std::priority_queue<Val, std::vector<Val>> L;
  std::priority_queue<Val, std::vector<Val>, std::greater<Val>> R;
  static constexpr Val xlow = std::numeric_limits<Val>::min() / 2;
  static constexpr Val xhigh = std::numeric_limits<Val>::max() / 2;
  Val add_l, add_r;
  Val max_l(){
    return L.empty() ? xlow : L.top() + add_l;
  }
  Val min_r(){
    return R.empty() ? xhigh : R.top() + add_r;
  }
  std::vector<std::tuple<Val, Val, Val>> __to_func(){
    std::vector<std::tuple<Val, Val, Val>> res;
    auto tmpl = L;
    Val d = 0, y1 = min_fx;
    while(!tmpl.empty()){
      Val x0 = tmpl.top() + add_l;
      Val x1 = res.empty() ? x0 : std::get<0>(res.back());
      Val c = -d * x1 + y1;
      res.push_back({x0, d, c});
      y1 = x0 * d + c;
      while(!tmpl.empty() && tmpl.top() + add_l == x0){
        tmpl.pop();
        d--;
      }
    }
    res.push_back({xlow, d, -d * (res.empty() ? 0 : std::get<0>(res.back())) + y1});
    std::reverse(res.begin(), res.end());
    auto tmpr = R;
    d = 0;
    Val y0 = min_fx;
    while(!tmpr.empty()){
      Val x0 = tmpr.top() + add_r;
      while(!tmpr.empty() && tmpr.top() + add_r == x0){
        tmpr.pop();
        d++;
      }
      res.push_back({x0, d, -d * x0 + y0});
      if(!tmpr.empty()){
        Val x1 = tmpr.top() + add_r;
        y0 = std::get<1>(res.back()) * x1 + std::get<2>(res.back());
      }
    }
    return res;
  }
public:
  slope_trick(): min_fx(0), add_l(0), add_r(0){}

  // 傾きが変わる点が単調な場合の線形時間構築
  // {a, b}
  // b = 0: add_left(a)
  // b = 1: add_right(a)
  // b = 2: add_abs(a)
  slope_trick(const std::vector<std::pair<Val, int>> &v): min_fx(0), add_l(0), add_r(0){
    std::vector<Val> Ltmp;
    std::deque<Val> Rtmp;
    auto addl = [&](Val x) -> void {
      if(!Rtmp.empty()){
        min_fx += x - Rtmp.front();
        Rtmp.push_back(x);
        x = Rtmp.front();
        Rtmp.pop_front();
      }
      Ltmp.push_back(x);
    };
    auto addr = [&](Val x) -> void {Rtmp.push_back(x);};
    for(auto [a, b] : v){
      if(b == 0) addl(a);
      else if(b == 1) addr(a);
      else addl(a), addr(a);
    }
    std::reverse(Ltmp.begin(), Ltmp.end());
    L = std::priority_queue<Val, std::vector<Val>>(Ltmp.begin(), Ltmp.end());
    R = std::priority_queue<Val, std::vector<Val>, std::greater<Val>>(Rtmp.begin(), Rtmp.end());
  }
  // 最小値
  Val min(){
    return min_fx;
  }
  // 最小値を取る区間[l, r]
  std::pair<Val, Val> min2(){
    return std::make_pair(max_l(), min_r());
  }
  // 定数xを足す
  void add_const(Val x){
    min_fx += x;
  }
  // max(0, a - x)を足す
  void add_left(Val a){
    Val r0 = min_r();
    if(a > r0){
      min_fx += a - r0;
      R.push(a - add_r);
      a = R.top() + add_r;
      R.pop();
    }
    L.push(a - add_l);
  }
  // max(0, x - a)を足す
  void add_right(Val a){
    Val l0 = max_l();
    if(a < l0){
      min_fx += l0 - a;
      L.push(a - add_l);
      a = L.top() + add_l;
      L.pop();
    }
    R.push(a - add_r);
  }
  // |x - a|を足す
  void add_abs(Val a){
    add_left(a);
    add_right(a);
  }
  // 左側をクリア
  void accumulate_min_left(){
    L.clear();
  }
  // 右側をクリア
  void accumulate_min_right(){
    R.clear();
  }
  // 左側を右にa, 右側を右にb平行移動する
  // l0 + a > r0 + bだと左右で交差して壊れる
  void shift(Val a, Val b){
    Val l0 = max_l();
    Val r0 = min_r();
    assert(r0 - l0 >= a - b);
    add_l += a;
    add_r += b;
  }
  // 傾きが変わる点で分割(区間[l, r], ax + bで表される線分の集合になる)
  // {l, a, b}
  std::vector<std::tuple<Val, Val, Val>> to_func(){
    return __to_func();
  }
};

template<typename Val>
struct slope_trick_monotone{
private:
  Val min_fx;
  std::vector<Val> L;
  std::deque<Val> R;
  static constexpr Val xlow = std::numeric_limits<Val>::min() / 2;
  static constexpr Val xhigh = std::numeric_limits<Val>::max() / 2;
  Val add_l, add_r;
  Val max_l(){
    return L.empty() ? xlow : L.back() + add_l;
  }
  Val min_r(){
    return R.empty() ? xhigh : R.front() + add_r;
  }
  // これまでの最大値
  Val max_p(){
    return !R.empty() ? R.back() + add_r : !L.empty() ? L.back() + add_l : xlow;
  }
  std::vector<std::tuple<Val, Val, Val>> __to_func(){
    std::vector<std::tuple<Val, Val, Val>> res;
    auto tmpl = L;
    Val d = 0, y1 = min_fx;
    while(!tmpl.empty()){
      Val x0 = tmpl.back() + add_l;
      Val x1 = res.empty() ? x0 : std::get<0>(res.back());
      Val c = -d * x1 + y1;
      res.push_back({x0, d, c});
      y1 = x0 * d + c;
      while(!tmpl.empty() && tmpl.back() + add_l == x0){
        tmpl.pop_back();
        d--;
      }
    }
    res.push_back({xlow, d, -d * (res.empty() ? 0 : std::get<0>(res.back())) + y1});
    std::reverse(res.begin(), res.end());
    auto tmpr = R;
    d = 0;
    Val y0 = min_fx;
    while(!tmpr.empty()){
      Val x0 = tmpr.front() + add_r;
      while(!tmpr.empty() && tmpr.front() + add_r == x0){
        tmpr.pop_front();
        d++;
      }
      res.push_back({x0, d, -d * x0 + y0});
      if(!tmpr.empty()){
        Val x1 = tmpr.front() + add_r;
        y0 = std::get<1>(res.back()) * x1 + std::get<2>(res.back());
      }
    }
    return res;
  }
public:
  slope_trick_monotone(): min_fx(0), add_l(0), add_r(0){}

  // 最小値
  Val min(){
    return min_fx;
  }
  // 最小値を取る区間[l, r]
  std::pair<Val, Val> min2(){
    return std::make_pair(max_l(), min_r());
  }
  // 定数xを足す
  void add_const(Val x){
    min_fx += x;
  }
  // max(0, a - x)を足す
  void add_left(Val a){
    assert(max_p() <= a);
    if(!R.empty()){
      min_fx += a - min_r();
      R.push_back(a - add_r); 
      a = R.front() + add_r;
      R.pop_front();
    }
    L.push_back(a - add_l);
  }
  // max(0, x - a)を足す
  void add_right(Val a){
    assert(max_p() <= a);
    R.push_back(a - add_r);
  }
  // |x - a|を足す
  void add_abs(Val a){
    add_left(a);
    add_right(a);
  }
  // 左側をクリア
  void accumulate_min_left(){
    L.clear();
  }
  // 右側をクリア
  void accumulate_min_right(){
    R.clear();
  }
  // 左側を右にa, 右側を右にb平行移動する
  // l0 + a > r0 + bだと左右で交差して壊れる
  void shift(Val a, Val b){
    Val l0 = max_l();
    Val r0 = min_r();
    assert(r0 - l0 >= a - b);
    add_l += a;
    add_r += b;
  }
  // 傾きが変わる点で分割(区間[l, r], ax + bで表される線分の集合になる)
  // {l, a, b}
  std::vector<std::tuple<Val, Val, Val>> to_func(){
    return __to_func();
  }
};

template<typename Val>
struct slope_trick_bbst{
private:
  Val min_fx;
  static constexpr Val xlow = std::numeric_limits<Val>::min() / 2;
  static constexpr Val xhigh = std::numeric_limits<Val>::max() / 2;
  multiset_sum<Val, Val> L, R;
  Val max_l(){
    return L.empty() ? xlow : L.max();
  }
  Val min_r(){
    return R.empty() ? xhigh : R.min();
  }
  std::vector<std::tuple<Val, Val, Val>> __to_func(){
    std::vector<std::tuple<Val, Val, Val>> res;
    auto tmpl = L.enumerate();
    std::reverse(tmpl.begin(), tmpl.end());
    Val d = 0, y1 = min_fx;
    for(int i = 0; i < tmpl.size();){
      auto [x0, cnt] = tmpl[i];
      Val x1 = res.empty() ? x0 : std::get<0>(res.back());
      Val c = -d * x1 + y1;
      res.push_back({x0, d, c});
      y1 = x0 * d + c;
      while(i < tmpl.size() && tmpl[i].first == x0) d -= tmpl[i++].second;
    }
    res.push_back({xlow, d, -d * (res.empty() ? 0 : std::get<0>(res.back())) + y1});
    std::reverse(res.begin(), res.end());
    auto tmpr = R.enumerate();
    d = 0;
    Val y0 = min_fx;
    for(int i = 0; i < tmpr.size();){
      auto [x0, cnt] = tmpr[i];
      while(i < tmpr.size() && tmpr[i].first == x0) d += tmpr[i++].second;
      res.push_back({x0, d, -d * x0 + y0});
      if(i < tmpr.size()){
        Val x1 = tmpr[i].first;
        y0 = std::get<1>(res.back()) * x1 + std::get<2>(res.back());
      }
    }
    return res;
  }
public:
  slope_trick_bbst(): min_fx(0){}
  // f(x)
  Val get(Val x){
    Val l0 = max_l();
    Val r0 = min_r();
    if(l0 <= x && x <= r0) return min_fx;
    // R中のx未満の要素(傾きが変わる点)をa1, a2...とすると
    // f(x) = min_fx + (x - a1) + (x - a2) + ...
    if(r0 < x){
      Val asum = R.sum_left(x); // x未満の要素の和
      Val acnt = R.low_count(x); // x未満の要素の数
      return min_fx + acnt * x - asum;
    }else{
      Val asum = L.sum_right(x); // x以上の要素の和
      Val acnt = L.up_count(x); // x以上の要素の数
      return min_fx + asum - acnt * x;
    }
  }
  // 最小値
  Val min(){
    return min_fx;
  }
  // 最小値を取る区間[l, r]
  std::pair<Val, Val> min2(){
    return std::make_pair(max_l(), min_r());
  }
  // [l, r)のmin
  Val range_min(Val l, Val r){
    Val l0 = max_l();
    Val r0 = min_r();
    if(r0 < l){
      return get(l);
    }else if(r <= l0){
      return get(r - 1);
    }else{
      return min_fx;
    }
  }
  // 定数xを足す
  void add_const(Val x){
    min_fx += x;
  }
  // max(0, k * (a - x))を足す
  void add_left(Val a, Val k = 1){
    assert(k > 0);
    if(min_r() < a){
      R.insert(a, k);
      auto v = R.split_left(k);
      min_fx += a * k - v->sum;
      L.merge_right(v);
    }else{
      L.insert(a, k);
    }
  }
  // max(0, k * (x - a))を足す
  void add_right(Val a, Val k = 1){
    assert(k > 0);
    if(a < max_l()){
      L.insert(a, k);
      auto v = L.split_right(k);
      min_fx += v->sum - a * k;
      R.merge_left(v);
    }else{
      R.insert(a, k);
    }
  }
  // k * |x - a|を足す
  void add_abs(Val a, Val k = 1){
    add_left(a, k);
    add_right(a, k);
  }
  // 左側をクリア
  void accumulate_min_left(){
    L.clear();
  }
  // 右側をクリア
  void accumulate_min_right(){
    R.clear();
  }
  // 左側を右にa, 右側を右にb平行移動する
  // l0 + a > r0 + bだと左右で交差して壊れる
  void shift(Val a, Val b){
    Val l0 = max_l();
    Val r0 = min_r();
    assert(r0 - l0 >= a - b);
    L.add_all(a);
    R.add_all(b);
  }
  // 傾きが変わる点で分割(区間[l, r], ax + bで表される線分の集合になる)
  // {l, a, b}
  std::vector<std::tuple<Val, Val, Val>> to_func(){
    return __to_func();
  }
};
#endif