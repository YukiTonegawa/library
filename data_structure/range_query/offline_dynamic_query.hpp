#ifndef _OFFLINE_DYNAMIC_QUERY_H_
#define _OFFLINE_DYNAMIC_QUERY_H_
#include <vector>
#include <cassert>
#include <algorithm>
#include <limits>
#include <functional>
#include <map>
// F: (query_id) -> void  取得
// G: 更新用の構造体 operator <, == が必要
// コンストラクタ (function<void(G&)>, function<void()>)  Gを入力とする更新用関数, rollback用関数

template<typename F, typename G>
struct offline_dynamic_query{
private:
  static constexpr int inf_time = 1 << 30;
  const std::function<void(G&)> update_func;
  const std::function<void()> rollback_func;

  std::multimap<G, int> E; // {g, inesrt_time}
  std::vector<std::pair<int, F>> Q;
  struct Update{ G g; int link_time, cut_time;};
  std::vector<Update> U{Update{G(), 0, 0}};

  void solve(int l, int r, int &j, int &k){
    assert(l != r);
    if(r - l == 1){
      while(k < Q.size() && Q[k].first == l) Q[k++].second(j++);
      return;
    }
    int mid = (l + r) / 2;
    int cnt = 0;
    for(int i = mid; i < r; i++) if(U[i].link_time <= l) cnt++, update_func(U[i].g);
    solve(l, mid, j, k);
    while(cnt) cnt--, rollback_func();
    for(int i = l; i <= mid; i++) if(r <= U[i].cut_time) cnt++, update_func(U[i].g);
    solve(mid, r, j, k);
    while(cnt) cnt--, rollback_func();
  }
public:
  offline_dynamic_query(std::function<void(G&)> ufunc, std::function<void()> rfunc): update_func(ufunc), rollback_func(rfunc){}

  void update(G g){
    int insert_time = U.size();
    E.emplace(g, insert_time);
    U.push_back(Update{g, insert_time, inf_time});
  }
  void update_unique(G g){
    int insert_time = U.size();
    auto itr = E.find(g);
    if(itr != E.end()) return;
    E.emplace(g, insert_time);
    U.push_back(Update{g, insert_time, inf_time});
  }
  // hを消す, ない場合は何もしない
  void erase(G g){
    auto itr = E.find(g);
    if(itr == E.end()) return;
    int link_time = itr->second;
    int cut_time = U.size();
    U[link_time].cut_time = cut_time;
    U.push_back({g, link_time, cut_time});
    E.erase(itr);
  }
  void query(F f){
    int query_time = (int)U.size() - 1;
    Q.push_back(std::make_pair(query_time, f));
  }
  void solve(){
    int j = 0, k = 0;
    solve(0, U.size(), j, k);
  }
};

/*
https://judge.yosupo.jp/problem/dynamic_graph_vertex_add_component_sum
int main(){
  using UF = union_find_commutative<range_add_range_sum_commutative, ll, ll>;
  io_init();
  int n, q;
  std::cin >> n >> q;
  auto b = read_vec<ll>(n);

  UF uf(b);
  offline_dynamic_query<function<void(int)>, ll> dc(
    [&](ll x){
      int f = x >> 60; // 1なら辺, 0なら加算
      int a = x >> 30;
      int b = x - ((ll)a << 30);
      if(!f) uf.unite(a, b);
      else uf.update(a, b);
    }, 
    [&](){
      uf.rollback();
    }
  );
  vector<ll> ans;
  queue<pair<int, int>> Q;
  range(i, 0, q){
    int x, y, z;
    std::cin >> x >> y;
    if(x != 3) std::cin >> z;
    if((x == 0 || x == 1) && y > z) swap(y, z);
    if(x == 0){
      dc.update(((ll)y << 30) + z);
    }else if(x == 1){
      dc.erase(((ll)y << 30) + z);
    }else if(x == 2){
      Q.push({y, z});
      dc.query([&](int qid){
        auto [cc, addval] = Q.front();
        Q.pop();
        uf.update(cc, addval);
      });
    }else{
      Q.push({y, -1});
      dc.query([&](int qid){
        int cc = Q.front().first;
        Q.pop();
        ans.push_back(uf.query_cc(cc));
      });
    }
  }
  dc.solve();
  range(i, 0, ans.size()) std::cout << ans[i] << '\n';
}
*/
#endif