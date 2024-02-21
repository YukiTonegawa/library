#ifndef _RANGE_TYPES_H_
#define _RANGE_TYPES_H_
#include "../rectangle_sum.hpp"
#include "../../basic/set_cache.hpp"
#include <unordered_map>

// O((N + Q)logN)
template<typename Val>
struct offline_static_range_types{
  int n;
  std::vector<std::pair<int, int>> E;
  std::vector<std::tuple<int, int, int>> Q;
  offline_static_range_types(const std::vector<Val> &val){
    std::unordered_map<Val, int> mp; // {値, 最後に現れたインデックス}
    n = val.size();
    for(int i = 0; i < n; i++){
      auto itr = mp.find(val[i]);
      int prev;
      if(itr == mp.end()){
        prev = -1;
        mp.emplace(val[i], i);
      }else{
        prev = itr->second;
        itr->second = i;
      }
      // lがprev+1以上になったら追加する
      E.push_back({prev + 1, i});
    }
  }
  // A[l, r)の種類数
  void query(int l, int r){
    int id = Q.size();
    Q.push_back({l, r, id});
  }
  std::vector<int> solve(){
    std::vector<int> ans(Q.size());
    std::sort(Q.begin(), Q.end());
    std::sort(E.begin(), E.end());
    binary_indexed_tree<int> bit(n);
    int j = 0;
    for(auto [l, r, id] : Q){
      while(j < E.size() && E[j].first <= l){
        bit.update(E[j].second, 1);
        if(E[j].first != 0) bit.update(E[j].first - 1, -1);
        j++;
      }
      ans[id] = bit.query(l, r);
    }
    return ans;
  }
};
/*
// O((N + Q)logN)
template<typename Val>
struct offline_static_range_types{
  int n;
  std::vector<std::pair<int, int>> E;
  std::vector<long long> Q;
  offline_static_range_types(const std::vector<Val> &val){
    n = val.size();
    E.resize(n);
    int j = 0;
    uhash_map<unsigned, int> mp(2 * n);
    for(int i = n - 1; i >= 0; i--){
      auto [f, x] = mp.at(val[i]);
      if(f) E[j++] = {i + 1, x};
      mp.emplace_replace(val[i], i);
    }
    for(auto [x, i] : mp.enumerate()) E[j++] = {0, i};
    std::reverse(E.begin(), E.end());
  }
  // A[l, r)の種類数
  void query(int l, int r){
    int id = Q.size();
    Q.push_back(((long long)l << 40) + ((long long)r << 20) + id);
  }
  std::vector<int> solve(){
    std::vector<int> ans(Q.size());
    std::sort(Q.begin(), Q.end());
    binary_indexed_tree<int> bit(n);
    int j = 0;
    for(long long x : Q){
      int l = x >> 40, r = (x >> 20) - ((long long)l << 20), id = x & ((1 << 20) - 1);
      while(j < E.size() && E[j].first <= l){
        bit.update(E[j].second, 1);
        if(E[j].first != 0) bit.update(E[j].first - 1, -1);
        j++;
      }
      ans[id] = bit.query(l, r);
    }
    return ans;
  }
};
*/

// O((N + Q)log^2 N)
template<typename Val>
struct offline_static_range_types2{
  int n;
  std::vector<std::tuple<int, int, Val>> E;
  std::vector<std::tuple<int, int, Val, Val, int>> Q;

  offline_static_range_types2(const std::vector<Val> &val){
    std::unordered_map<Val, int> mp; // {値, 最後に現れたインデックス}
    n = val.size();
    for(int i = 0; i < n; i++){
      auto itr = mp.find(val[i]);
      int prev;
      if(itr == mp.end()){
        prev = -1;
        mp.emplace(val[i], i);
      }else{
        prev = itr->second;
        itr->second = i;
      }
      // lがprev+1以上になったら追加する
      E.push_back({prev + 1, i, val[i]});
    }
  }
  // A[l, r)の値が[s, t)の要素の種類数
  void query(int l, int r, Val s, Val t){
    int id = Q.size();
    Q.push_back({l, r, s, t, id});
  }
  std::vector<int> solve(){
    std::vector<int> ans(Q.size());
    std::sort(Q.begin(), Q.end());
    std::sort(E.begin(), E.end());
    offline_point_add_rectangle_sum<Val, int> rect;
    int j = 0;
    for(auto [l, r, s, t, id] : Q){
      while(j < E.size() && std::get<0>(E[j]) <= l){
        auto [prev, cur, v] = E[j];
        rect.update(cur, v, 1);
        if(prev != 0) rect.update(prev - 1, v, -1);
        j++;
      }
      ans[id] = rect.query(l, r, s, t);
    }
    return ans;
  }
};

// 十分大きな値INFとして, 初期値は A = {INF, INF + 1, INF + 2, INF + 3...}, (つまり初期状態はN種類)
// O((N + Q)log^2(N + Q))
template<typename Val>
struct offline_point_set_range_types{
private:
  struct query_type{
    int t, id, l, r;
    Val val;
    query_type(int a, int b, int c, int d, Val e): t(a), id(b), l(c), r(d), val(e){}
  };
  std::vector<query_type> q;
  int qs = 0;
  std::vector<int> solve_inner(){
    std::vector<int> ans(qs, 0);
    offline_point_add_rectangle_sum<int, int> rect;
    std::vector<set_avl_cache<int>> idx_table;
    int max_idx = -1;
    // 追加される値で座圧
    std::vector<Val> zrev;
    for(int i = 0; i < q.size(); i++){
      if(q[i].t == 0){
        zrev.push_back(q[i].val);
        max_idx = std::max(max_idx, q[i].id);
      }
    }
    std::sort(zrev.begin(), zrev.end());
    zrev.erase(std::unique(zrev.begin(), zrev.end()), zrev.end());
    std::vector<int> v(max_idx + 1, -1);
    idx_table.resize(zrev.size());
    int inf = set_avl_cache<int>::inf;
    for(auto qi : q){
      if(qi.t == 0){
        int i = qi.id, y = v[i];
        int z = std::lower_bound(zrev.begin(), zrev.end(), qi.val) - zrev.begin();
        if(z == y) continue;
        if(y != -1){
          int prev = idx_table[y].lower_bound_rev(i - 1);
          int next = idx_table[y].lower_bound(i + 1);
          if(prev != inf && next != inf) rect.update(prev, next, 1);
          if(prev != inf) rect.update(prev, i, -1);
          if(next != inf) rect.update(qi.id, i, -1);
          idx_table[y].erase(i);
        }
        int prev = idx_table[z].lower_bound_rev(i - 1);
        int next = idx_table[z].lower_bound(i + 1);
        if(prev != inf && next != inf) rect.update(prev, next, -1);
        if(prev != inf) rect.update(prev, i, 1);
        if(next != inf) rect.update(i, next, 1);
        idx_table[z].insert(i);
        v[i] = z;
      }else{
        rect.query(qi.l, qi.r, qi.l, qi.r);
      }
    }
    auto rect_sum = rect.solve();
    int j = 0;
    for(auto qi : q){
      if(qi.t == 1){
        ans[qi.id] = qi.r - qi.l - rect_sum[j++];
      }
    }
    return ans;
  }
public:
  offline_point_set_range_types(){}
  // A[i] <- val
  void update(int i, Val val){
    q.push_back(query_type(0, i, 0, 0, val));
  }
  // A[l, r)の種類数
  void query(int l, int r){
    q.push_back(query_type(1, qs++, l, r, 0));
  }
  std::vector<int> solve(){
    return solve_inner();
  }
};

#endif