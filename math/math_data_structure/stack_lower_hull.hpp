#ifndef _STACK_LOWER_HULL_H_
#define _STACK_LOWER_HULL_H_
#include <cassert>
#include <vector>
#include <array>

template<typename Val, typename Val2>
struct stack_lower_hull{
private:
  static constexpr int log_max_n = 20;
  std::vector<std::array<int, log_max_n>> par; // 凸包上の頂点を2^i個遡った時のインデックス
  std::vector<std::pair<Val, Val>> P;
  std::vector<int> sz;
  bool ccw(int a, int b, int c){
    return (Val2)(P[c].second - P[a].second) * (P[b].first - P[a].first) > (Val2)(P[b].second - P[a].second) * (P[c].first - P[a].first);
  }
  // n <= 2^kとなるような最小のk
  int ceillog2(int n){
    if(n == 0) return 0;
    int k = 63 - __builtin_clzll(n);
    return n > (1 << k) ? k + 1 : k;
  }
public:
  std::pair<Val, Val> get(int i){return P[i];}
  // O(logN)
  void push(std::pair<Val, Val> p){
    assert(P.empty() || P.back().first < p.first);
    int k = P.size() - 1;
    P.push_back(p);
    if(k != -1 && k != 0){
      if(!ccw(par[k][0], k, P.size() - 1)){
        for(int i = ceillog2(sz.back()); i >= 0; i--){
          int b = par[k][i];
          if(b == -1 || par[b][0] == -1) continue;
          if(!ccw(par[b][0], b, P.size() - 1)) k = b;
        }
        k = par[k][0];
      }
    }
    sz.push_back(k == -1 ? 1 : sz[k] + 1);
    par.push_back({});
    par.back().fill(-1);
    par.back()[0] = k;
    for(int i = 1; i < log_max_n; i++){
      if(par.back()[i - 1] == -1) break;
      par.back()[i] = par[par.back()[i - 1]][i - 1];
    }
  }
  // O(1)
  void pop(){
    P.pop_back();
    par.pop_back();
  }
  // O(1)
  int size_points(){
    return P.size();
  }
  // O(1)
  int size_hull(){
    if(sz.empty()) return 0;
    return sz.back();
  }
  // 全体の下側凸包でx座標の小さい順にk個目, O(logN)
  int kth_point(int k){
    assert(0 <= k && k < size_hull());
    k = size_hull() - 1 - k;
    int pos = (int)P.size() - 1;
    for(int i = ceillog2(sz.back()); i >= 0; i--){
      if((1 << i) > k) continue;
      pos = par[pos][i];
      k -= 1 << i;
    }
    return pos;
  }
  // 全体の下側凸包にi番目の点が含まれるか(含まれない場合は-1, 含まれる場合はその順番)
  // O(logN)
  int order(int i){
    aassert(0 <= i && i < P.size());
    int pos = P.size() - 1;
    int res = sz.back() - 1;
    for(int j = ceillog2(sz.back()); j >= 0; j--){
      if(par[pos][j] < i) continue;
      res -= 1 << j;
      if(par[pos][j] == i) return res;
      pos = par[pos][j];
    }
    return (pos == i ? res : -1);
  }
  // [0, i]の点だけで作った下側凸包でiの左の点(ない場合は-1)
  // iが全体の下側凸包に含まれる場合この左の点も含まれ, 隣り合っている
  // O(1)
  int left(int i){
    return par[i][0];
  }
  // [0, i]の点だけで作った下側凸包でiのk個左の点(ない場合は-1)
  // O(logN)
  int left_k(int i, int k){
    for(int j = ceillog2(sz[i]); j >= 0; j--){
      if((1 << j) > k) continue;
      k -= (1 << j);
      i = par[i][j];
      if(i == -1 || k == 0) return i;
    }
    return -1;
  }
  // [0, r)の点にax+bを足したと仮定した時の[0, r)の{最小値を取る点, 最小値}
  // O(logN)
  std::pair<int, Val> linear_add_min(int r, Val a, Val b){
    assert(0 < r);
    int pos = r - 1;
    for(int i = ceillog2(sz[r - 1]); i >= 0; i--){
      int mid = par[pos][i];
      if(mid == -1) continue;
      if((P[mid].first - P[pos].first) * a < (P[pos].second - P[mid].second)){
        pos = mid;
      }
    }
    return {pos, P[pos].first * a + P[pos].second};
  }
};
#endif