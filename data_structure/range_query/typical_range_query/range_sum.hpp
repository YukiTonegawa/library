#ifndef _RANGE_SUM_H_
#define _RANGE_SUM_H_
#include <vector>
#include "../../segment_tree/binary_indexed_tree_compressed.hpp"

template<typename Idx, typename Val>
struct offline_static_range_add_range_sum{
  std::vector<std::tuple<Idx, Idx, Val>> upd;
  std::vector<std::pair<Idx, Idx>> q;
  bool built = false;
  void update(Idx l, Idx r, Val x){
    assert(!built);
    upd.push_back({l, r, x});
  }
  void query(Idx l, Idx r){
    assert(!built);
    q.push_back({l, r});
  }
  std::vector<Val> solve(){
    assert(!built);
    built = true;
    std::vector<Idx> I;
    for(auto [l, r, x] : upd){
      I.push_back(l);
      I.push_back(r);
    }
    std::sort(allof(I));
    I.erase(std::unique(I.begin(), I.end()), I.end());
    auto lb = [&](Idx k){
      return std::lower_bound(I.begin(), I.end(), k) - I.begin();
    };
    std::vector<Val> dx(I.size(), 0);
    std::vector<Val> sum(I.size(), 0); // sum[r] := [0, I[r])の和
    for(auto [l, r, x] : upd){
      dx[lb(l)] += x;
      dx[lb(r)] -= x;
    }
    for(int i = 1; i < I.size(); i++){
      dx[i] += dx[i - 1];
      sum[i] = sum[i - 1] + dx[i - 1] * (I[i] - I[i - 1]);
    }
    std::vector<Val> ans;
    auto get = [&](Idx r){
      int i = lb(r);
      if(i == 0) return 0;
      i--;
      return sum[i] + dx[i] * (r - I[i]);
    };
    for(auto [l, r] : q) ans.push_back(get(r) - get(l));
    return ans;
  }
};

template<typename Idx, typename Val>
struct offline_range_add_range_sum{
  std::vector<bool> is_query;
  std::vector<std::tuple<Idx, Idx, Val>> upd;
  std::vector<std::pair<Idx, Idx>> q;
  bool built = false;
  void update(Idx l, Idx r, Val x){
    assert(!built);
    upd.push_back({l, r, x});
    is_query.push_back(0);
  }
  void query(Idx l, Idx r){
    assert(!built);
    q.push_back({l, r});
    is_query.push_back(1);
  }
  std::vector<Val> solve(){
    assert(!built);
    built = true;
    std::vector<Idx> I;
    for(auto [l, r, x] : upd){
      I.push_back(l);
      I.push_back(r);
    }
    for(auto [l, r] : q){
      I.push_back(l);
      I.push_back(r);
    }
    std::sort(allof(I));
    I.erase(std::unique(I.begin(), I.end()), I.end());
    auto lb = [&](Idx k){
      return std::lower_bound(I.begin(), I.end(), k) - I.begin();
    };
    binary_indexed_tree_compressed<Idx, Val> bit(I.size());
    std::vector<Val> ans;
    int i = 0, j = 0;
    for(bool f : is_query){
      if(!f){
        auto [l, r, x] = upd[i++];
        bit.update(l, r, lb(l), lb(r), x);
      }else{
        auto [l, r] = q[j++];
        ans.push_back(bit.query(l, r, lb(l), lb(r)));
      }
    }
    return ans;
  }
};


#endif