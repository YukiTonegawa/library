#ifndef _RANGE_MEX_H_
#define _RANGE_MEX_H_
#include <set>
#include "../../segment_tree/segment_tree2d.hpp"
#include "../rectangle_min.hpp"

// 隣接する同じ要素同士に辺を張る
// 区間を2次元平面上の点と見なすと, [l, r)のmex は [-1, l) × [r, n]の点の値のminを求める問題に言い換えられる
struct offline_static_range_mex{
private:
  int n;
  offline_static_upper_rectangle_min<int, chmin_struct<int>> rect;
  void build(const std::vector<int> &v){
    n = v.size();
    // mexはn以下なので, それより大きい要素は無視していい
    std::vector<int> last(n + 1, -1); // 最後に現れるインデックス
    for(int i = 0; i < n; i++){
      assert(v[i] >= 0);
      if(v[i] > n) continue;
      // {a_last, a_i}
      rect.update(last[v[i]] + 1, i, v[i]);
      last[v[i]] = i;
    }
    for(int i = 0; i <= n; i++) rect.update(last[i] + 1, n, i);
  }
public:
  offline_static_range_mex(){}
  offline_static_range_mex(const std::vector<int> &v){build(v);}
  // mex(A[l, r))
  void query(int l, int r){
    rect.query(l + 1, r, n + 1);
  }
  std::vector<int> solve(){
    return rect.solve();
  }
};

struct static_range_mex{
private:
  int n;
  compressed_segment_tree2d<int, range_min<int>> rect;
  void build(const std::vector<int> &v){
    n = v.size();
    std::vector<std::tuple<int, int, int>> P;
    // mexはn以下なので, それより大きい要素は無視していい
    std::vector<int> last(n + 1, -1); // 最後に現れるインデックス
    for(int i = 0; i < n; i++){
      assert(v[i] >= 0);
      if(v[i] > n) continue;
      // {a_last, a_i}
      P.push_back({last[v[i]], i, v[i]});
      last[v[i]] = i;
    }
    for(int i = 0; i <= n; i++) P.push_back({last[i], n, i});
    rect = compressed_segment_tree2d<int, range_min<int>>(P);
  }
public:
  static_range_mex(){}
  static_range_mex(const std::vector<int> &v){build(v);}
  // mex(A[l, r)), O(log^2N)
  int query(int l, int r){
    return rect.query(-1, l, r, n + 1);
  }
};


struct incremental_mex{
  int mex, maxQ;
  std::vector<bool> flag;
  // maxQ := 追加クエリの最大回数
  incremental_mex(int maxQ): mex(0), maxQ(maxQ), flag(maxQ + 1, false){}
  // 追加 O(1)
  void push(long long x){
    assert(0 <= x);
    if(x > maxQ) return;
    flag[x] = 1;
  }
  // 全体のmex ならしO(1)
  int query(){
    while(mex < flag.size() && flag[mex]) mex++;
    assert(mex != flag.size());
    return mex;
  }
};
#endif