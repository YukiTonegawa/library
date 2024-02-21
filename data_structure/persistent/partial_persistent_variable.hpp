#ifndef _PARTIAL_PERSISTENT_VARIABLE_H_
#define _PARTIAL_PERSISTENT_VARIABLE_H_
#include <vector>
#include <cassert>

template<typename Val>
struct partial_persistent_variable{
  static constexpr int brute_force_threshold = 6;
  std::vector<std::pair<int, Val>> v;
  partial_persistent_variable(){}
  partial_persistent_variable(Val x){v.push_back({0, x});}

  // t以前(t含む)に最後に追加された要素, 無いと壊れる
  Val get(int t){
    int n = v.size();
    if(n < brute_force_threshold){
      for(int i = 0; i < n - 1; i++) if(t < v[i + 1].first) return v[i].second;
      return v[n - 1].second;
    }else{
      int l = 0, m, r = n;
      while(r - l > 1){
        m = (l + r) >> 1;
        if(t < v[m].first) r = m;
        else l = m;
      }
      return v[l].second;
    }
  }
  // 最後に追加された時刻, 一度も追加されてない場合は-1
  int last_update(){
    return v.empty() ? -1 : v.back().first;
  }
  // 最後に追加された要素
  Val get_new(){
    assert(!v.empty());
    return v.back().second;
  }
  // 時刻tにxを追加, これまでに追加された時刻が非減少でなくてはならない
  void set(int t, Val x){
    assert(v.empty() || v.back().first <= t);
    if(v.empty() || v.back().first != t) v.push_back({t, x});
    else v.back().second = x;
  }
  // t時点に戻す(tより後に追加された要素を全て消す)
  void rollback(int t){
    while(!v.empty()){
      if(v.back().first <= t) return;
      v.pop_back();
    }
  }
  // 関数fが時刻t未満全てfalse, t以降全てtrueになるようなt, ない場合は-1
  template<typename F>
  int bisect(F &f){
    int n = v.size();
    if(n < brute_force_threshold){
      for(int i = 0; i < n; i++) if(f(v[i].second)) return v[i].first;
      return -1;
    }else{
      int l = -1, m, r = n;
      while(r - l > 1){
        m = (l + r) >> 1;
        if(f(v[m].second)) r = m;
        else l = m;
      }
      return (r == n ? -1 : r);
    }
  }
};
#endif