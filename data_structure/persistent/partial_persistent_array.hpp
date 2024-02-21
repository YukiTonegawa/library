#ifndef _PPARRAY_H_
#define _PPARRAY_H_
#include <vector>
#include <cassert>

template<typename Val>
struct partial_persistent_array{
private:
  // 要素を二分探索で探す時に, 要素数が少ない場合は全探索
  static constexpr int brute_force_threshold = 6;
  int T = 0;
  struct history{
    std::vector<std::pair<int, Val>> v;
    history(){}
    history(Val x){v.push_back({0, x});}
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
    Val get_new(){
      return v.back().second;
    }
    void set(int t, Val x){
      if(v.empty() || v.back().first != t) v.push_back({t, x});
      else v.back().second = x;
    }
    void rollback(int t){
      while(!v.empty()){
        if(v.back().first <= t) return;
        v.pop_back();
      }
    }
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
  std::vector<history> hist;
public:
  partial_persistent_array(){}
  partial_persistent_array(int n, Val x): hist(n){for(int i = 0; i < n; i++) hist[i] = history(x);};
  partial_persistent_array(const std::vector<Val> &v): hist(v.size()){for(int i = 0; i < v.size(); i++) hist[i] = history(v[i]);}

  int time(){
    return T;
  }
  void next(){
    T++;
  }
  Val get(int t, int i){
    return hist[i].get(t);
  }
  Val get_new(int i){
    return hist[i].get_new();
  }
  // f(elem(t未満)) = falseであり, f(elem(t以上)) = trueになるようなt
  // t = 0ですでに満たす場合は0
  // t = Tで満たさない場合は-1
  template<typename F>
  int bisect(int i, const F &f){
    return hist[i].bisect(f);
  }
  void set(int i, Val x){
    hist[i].set(T, x);
    g.push_back(i);// rollbackする場合必要
  }
  // t時点に戻す(tより後に追加された要素を全て消す)
  std::vector<int> g;
  void rollback(int t){
    assert(0 <= t);
    while(g.size()){
      hist[g.back()].rollback(t);
      g.pop_back();
    }
    T = t;
  }
};
#endif