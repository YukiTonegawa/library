#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

template<typename S>
struct mo_dynamic{
private:
  const int block_size;
  using Val = typename S::Val;
  using Ret = typename S::Ret;
  struct query{
    int i, l, r, lb, rb, t;
    query(int i, int l, int r, int lb, int rb, int t)
    : i(i), l(l), r(r), lb(lb), rb(rb), t(t){}
  };
  S s;
  int pq, pl, pr, pt;
  std::vector<query> q;
public:
  mo_dynamic(std::vector<Val> &v)
  : block_size(round(pow(int(v.size()), 2.0 / 3.0))), s(v), pq(0), pl(0), pr(0), pt(0){}
  inline int insert_query(int l, int r){
    int qsz = q.size();
    q.emplace_back(query(qsz, l, r, l / block_size, r / block_size, pt));
    return qsz;
  }
  inline void insert_update(int i, Val x){
    s.update_query.push_back({i, {x, x}});
    pt++;
  }
  std::vector<Ret> solve(){
    pt = 0;
    std::sort(q.begin(), q.end(), [&](query &a, query &b){
      if(a.lb != b.lb) return a.lb < b.lb;
      if(a.rb != b.rb) return a.lb & 1 ? a.rb > b.rb : a.rb < b.rb;
      return a.rb & 1 ? a.t > b.t : a.t < b.t;
    });
    s.simulate();
    std::vector<Ret> ret(q.size());
    for(query &qi:q){
      while(pt < qi.t) s.time_add(pl, pr, pt++);
      while(pt > qi.t) s.time_sub(pl, pr, --pt);
      while(pl > qi.l) s.add(--pl);
      while(pr < qi.r) s.add(pr++);
      while(pl < qi.l) s.del(pl++);
      while(pr > qi.r) s.del(--pr);
      ret[qi.i] = s.get();
    }
    return ret;
  }
};

struct S{
  using Val = int;
  using Ret = long long;
  std::vector<Val> A;
  std::vector<std::pair<int, std::pair<Val, Val>>> update_query;
  S(std::vector<Val> &A):A(A){}
  void simulate(){
    std::vector<Val> sim = A;
    for(int i=0;i<update_query.size();i++){
      auto [idx, p] = update_query[i];
      Val pv = sim[idx], nw = update_func(sim[idx], p.first);
      sim[idx] = update_func(sim[idx], p.first);
      update_query[i].second = {pv, nw};
    }
  }
  Val update_func(Val old, Val x){
    return old + x;
  }
  void time_add(int pl, int pr, int qi){
    auto [idx, p] = update_query[qi];
    if(pl <= idx && idx < pr){
      del(idx);
      A[idx] = p.second;
      add(idx);
    }else{
      A[idx] = p.second;
    }
  }
  void time_sub(int pl, int pr, int qi){
    auto [idx, p] = update_query[qi];
    if(pl <= idx && idx < pr){
      del(idx);
      A[idx] = p.first;
      add(idx);
    }else{
      A[idx] = p.first;
    }
  }
  long long sum = 0;
  void add(int idx){
    sum += A[idx];
  }
  void del(int idx){
    sum -= A[idx];
  }
  Ret get(){
    return sum;
  }
};
