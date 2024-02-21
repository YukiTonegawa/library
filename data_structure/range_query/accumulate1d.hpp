#ifndef _ACCUMULATE1d_H_
#define _ACCUMULATE1d_H_
#include <vector>
#include <algorithm>

template<typename Val>
struct accumulate1d{
  std::vector<Val> sum;
  accumulate1d(){}
  accumulate1d(const std::vector<Val> &v): sum(v){
    for(int i = 1; i < v.size(); i++) sum[i] = sum[i - 1] + v[i];
  }
  // [0, r)の和, 範囲外の部分は全て単位元
  Val query(int r){
    r = std::min(r, (int)sum.size());
    if(r <= 0) return 0;
    return sum[r - 1];
  }
  // [l, r)の和, 範囲外の部分は全て単位元
  Val query(int l, int r){
    l = std::max(l, 0);
    r = std::min(r, (int)sum.size());
    if(r <= l) return 0;
    return (l == 0 ? sum[r - 1] : (sum[r - 1] - sum[l - 1]));
  }
  void push_back(Val x){
    Val y = (sum.empty() ? 0 : sum.back());
    sum.push_back(y + x);
  }
  void pop_back(){
    assert(!sum.empty());
    sum.pop_back();
  }
  // [0, k]がx以上になる最小インデックス, ない場合はn
  int lower_bound(Val x){
    return std::lower_bound(sum.begin(), sum.end(), x) - sum.begin();
  }
  // [0, k]がxより大きくなる最小インデックス, ない場合はn
  int upper_bound(Val x){
    return std::upper_bound(sum.begin(), sum.end(), x) - sum.begin();
  }
};
#endif
