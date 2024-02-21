#ifndef _OFFLINE_DYNAMIC_BRIDGE_CONNECTIVITY_H_
#define _OFFLINE_DYNAMIC_BRIDGE_CONNECTIVITY_H_
#include "rollback_bridge_connectivity.hpp"
#include "../../data_structure/range_query/pseudo_tree.hpp"
#include <map>

struct offline_dynamic_bridge_connectivity{
private:
  int bridge = 0;
  pseudo_segment_tree<int> segment;
  rollback_bridge_connectivity bc;
  enum query_type: int{_bridge, _relation, _is_tree, _size, _is_tec, _link, _cut};
  struct query{
    int tm;
    query_type type;
    int u, v, idx;
    query(int tm, query_type t, int u, int v, int idx) : tm(tm), type(t), u(u), v(v), idx(idx){}
  };
  struct update{
    query_type type;
    int a, b, idx;
    update(query_type t, int a, int b, int idx): type(t), a(a), b(b), idx(idx){}
  };
  struct answer{
    query_type type;
    int a, b, c, re;
    answer(query_type t): type(t), a(-1), b(-1), c(-1), re(0){}
  };
  std::vector<update> U;
  std::vector<query> Q;
  std::vector<answer> ans;
  using node = rollback_bridge_connectivity::node;
  std::vector<node*> V;

  void _solve(){
    segment = pseudo_segment_tree<int>(U.size());
    std::multimap<std::pair<int, int>, int> E;// (s, t), time
    std::vector<std::vector<std::pair<int, int>>> e(segment.M * 2 - 1);
    for(int i = 0; i < U.size(); i++){
      if(U[i].a > U[i].b) std::swap(U[i].a, U[i].b);
      if(U[i].type == _link){
        E.emplace(std::make_pair(U[i].a, U[i].b), i);
      }else{
        auto itr = E.find(std::make_pair(U[i].a, U[i].b));
        // 辺が無かった
        if(itr == E.end()){
          ans[U[i].idx].a = -1;
          ans[U[i].idx].b = U[i].idx;
          continue;
        }
        int link_time = U[itr->second].idx;
        int cut_time = U[i].idx;
        ans[link_time].a = ans[cut_time].a = link_time,
        ans[link_time].b = ans[cut_time].b = cut_time;
        auto seg_idx = segment.range_to_index(itr->second, i);
        for(int idx: seg_idx) e[idx].push_back(std::make_pair(U[i].a, U[i].b));
        E.erase(itr);
      }
    }
    for(auto itr = E.begin(); itr != E.end(); itr++){
      int link_time = itr->second;
      int cut_time = segment.N;
      auto seg_idx = segment.range_to_index(link_time, cut_time);
      std::pair<int, int> p = itr->first;
      for(int idx: seg_idx) e[idx].push_back(p);
      ans[link_time].a = link_time;
      ans[link_time].b = U.size();
    }
    int qidx = 0;
    while(qidx < Q.size() && Q[qidx].tm == -1){
      query &q = Q[qidx];
      if(q.type == _bridge){
        ans[q.idx].re = bridge;
      }else if(q.type == _relation){
        ans[q.idx].re = bc.relation(V[q.u], V[q.v]);
      }else if(q.type == _is_tree){
        ans[q.idx].re = bc.is_tree(V[q.u]);
      }else if(q.type == _size){
        ans[q.idx].re = bc.components_size(V[q.u]);
      }else if(q.type == _is_tec){
        ans[q.idx].re = bc.is_two_edge_connected(V[q.u]);
      }
      qidx++;
    }
    auto dfs = [&](auto &&dfs, int k, int l, int r)->void{
      int mid = (l + r) / 2;
      if(r - l == 1){
        for(std::pair<int, int> _e : e[k]) bridge += bc.link(V[_e.first], V[_e.second]).second;
        while(qidx < Q.size() && Q[qidx].tm <= l){
          query &q = Q[qidx];
          if(q.type == _bridge){
            ans[q.idx].c = bridge;
          }else if(q.type == _relation){
            ans[q.idx].re = bc.relation(V[q.u], V[q.v]);
          }else if(q.type == _is_tree){
            ans[q.idx].re = bc.is_tree(V[q.u]);
          }else if(q.type == _size){
            ans[q.idx].re = bc.components_size(V[q.u]);
          }else if(q.type == _is_tec){
            ans[q.idx].re = bc.is_two_edge_connected(V[q.u]);
          }
          qidx++;
        }
        for(int i = 0; i < e[k].size(); i++) bridge += bc.rollback();
      }else{
        for(std::pair<int, int> _e : e[k]) bridge += bc.link(V[_e.first], V[_e.second]).second;
        dfs(dfs, k * 2 + 1, l, mid);
        dfs(dfs, k * 2 + 2, mid, r);
        for(int i = 0; i < e[k].size(); i++) bridge += bc.rollback();
      }
    };
    dfs(dfs, 0, 0, segment.M);
  }
public:
  int N;
  offline_dynamic_bridge_connectivity(int n): V(n), N(n){
    for(int i = 0; i < n; i++) V[i] = bc.make_node();
  }
  // 橋の数
  void count_bridge(){
    Q.push_back(query(int(U.size()) - 1, _bridge, 0, 0, ans.size()));
    ans.push_back(answer(_bridge));
  }
  // aとbの関係　
  // 0: 非連結, 1: 連結, 2: 2辺連結
  void relation(int a, int b){
    assert(a < N && b < N);
    Q.push_back(query(int(U.size()) - 1, _relation, a, b, ans.size()));
    ans.push_back(answer(_relation));
  }
  // aを含む連結成分が木か
  void is_tree(int a){
    assert(a < N);
    Q.push_back(query(int(U.size()) - 1, _is_tree, a, a, ans.size()));
    ans.push_back(answer(_is_tree));
  }
  // aを含む連結成分のサイズ
  void size(int a){
    assert(a < N);
    Q.push_back(query(int(U.size()) - 1, _size, a, a, ans.size()));
    ans.push_back(answer(_size));
  }
  // aが2辺連結成分に含まれるか
  void is_two_edge_connected(int a){
    assert(a < N);
    Q.push_back(query(int(U.size()) - 1, _is_tec, a, a, ans.size()));
    ans.push_back(answer(_is_tec));
  }
  // 辺a-bを追加する, すでにあっても追加する
  void link(int a, int b){
    assert(a < N && b < N);
    if(a > b) std::swap(a, b);
    U.push_back(update(_link, a, b, ans.size()));
    Q.push_back(query(int(U.size()) - 1, _bridge, a, b, ans.size()));
    ans.push_back(answer(_link));
  }
  // 辺a-bを消す, ない場合は何もしない
  void cut(int a, int b){
    assert(a < N && b < N);
    if(a > b) std::swap(a, b);
    U.push_back(update(_cut, a, b, ans.size()));
    Q.push_back(query(int(U.size()) - 1, _bridge, a, b, ans.size()));
    ans.push_back(answer(_cut));
  }
  std::vector<answer> solve(){
    _solve();
    return ans;
  }
};
// ansの種類
// enum query_type: int{_bridge, _relation, _is_tree, _size, _is_tec, _link, _cut};

//_bridge: {c} = 橋の数
//_relation: {re} 0: 非連結, 1: 連結, 2: 2辺連結
// _is_tree: {re} aを含む連結成分が木か
// _size: {re} aを含む連結成分のサイズ
// _is_tec: {re} aが2辺連結成分に含まれるか
//_link: {a, b, re} = 辺のlink, cutされたクエリ番号, 削除されなかった場合はb = query.size(), re : 操作後の橋の数　
//_cut: {a, b, re} = 辺のlink, cutされたクエリ番号, 辺がなかった場合はa = -1,  re : 操作後の橋の数
#endif