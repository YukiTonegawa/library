#ifndef _RECTANGLE_PLANE_ADD_H_
#define _RECTANGLE_PLANE_ADD_H_
#include <vector>
#include <cassert>
#include <algorithm>
template<typename Val>
struct binary_indexed_tree3{
  int N;
  std::vector<std::array<Val, 3>> sum;
  binary_indexed_tree3(){}
  binary_indexed_tree3(int N): N(N), sum(N + 1 , {0, 0, 0}){}
  void update(int k, Val a, Val b, Val c){
    for(int i = k + 1; i <= N; i += (i & (-i))) sum[i][0] += a, sum[i][1] += b, sum[i][2] += c;
  }
  std::tuple<Val, Val, Val> query(int r){
    Val a = 0, b = 0, c = 0;
    for(int k = r; k > 0; k -= (k & (-k))) a += sum[k][0], b += sum[k][1], c += sum[k][2];
    return {a, b, c};
  }
};

template<typename Idx, typename Val>
struct offline_static_plane_add{
private:
  struct Query{
    Idx x, y;
    int id;
  };
  struct Update{
    Idx lx, rx, ly, ry;
    Val a, b, c;
  };
  std::vector<Update> U;
  std::vector<Query> Q;
public:
  // lx <= x < rx, ly <= y < ryを満たす全ての(x, y)にax + by + cを足す
  void update(Idx lx, Idx rx, Idx ly, Idx ry, Val a, Val b, Val c){U.push_back(Update{lx, rx, ly, ry, a, b, c});}
  // get(x, y)
  void query(Idx x, Idx y){Q.push_back(Query{x, y, (int)Q.size()});}
  std::vector<Val> solve(){
    int N = Q.size();
    std::vector<Val> ans(N, 0);
    if(U.empty() || Q.empty()) return ans;
    std::vector<Idx> Y;
    for(const Update &u : U){
      Y.push_back(u.ly);
      Y.push_back(u.ry);
    }
    std::sort(Y.begin(), Y.end());
    Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
    int M = U.size();
    for(int i = 0; i < M; i++){
      Update &u = U[i];
      int ly = std::lower_bound(Y.begin(), Y.end(), u.ly) - Y.begin();
      int ry = std::lower_bound(Y.begin(), Y.end(), u.ry) - Y.begin();
      u.ly = ly;
      u.ry = ry;
      U.push_back(Update{u.rx, u.rx, ly, ry, -u.a, -u.b, -u.c});
    }
    std::sort(U.begin(), U.end(), [](const Update &a, const Update &b){return a.lx < b.lx;});
    std::sort(Q.begin(), Q.end(), [](const Query &a, const Query &b){return a.x < b.x;});
    binary_indexed_tree3<Val> bit(Y.size());
    int u = 0, q = 0;
    while(q < N){
      if(u == U.size() || Q[q].x < U[u].lx){
        int y = std::upper_bound(Y.begin(), Y.end(), Q[q].y) - Y.begin();
        auto [a, b, c] = bit.query(y);
        ans[Q[q].id] = a * Q[q].x + b * Q[q].y + c;
        q++;
      }else{
        Val a = U[u].a, b = U[u].b, c = U[u].c;
        bit.update(U[u].ly, a, b, c);
        bit.update(U[u].ry, -a, -b, -c);
        u++;
      }
    }
    return ans;
  }
};

template<typename Idx, typename Val>
struct offline_dynamic_plane_add{
private:
  static constexpr int qlim = 1e8;
  struct Query{
    Idx x, y;
    int id;
  };
  struct Update{
    Idx lx, rx, ly, ry;
    Val a, b, c;
  };
  struct Event{
    Idx x;
    int lyc, ryc; // ryc >= qlim  -> 取得クエリ
    Val a, b, c;
  };
  std::vector<Update> U;
  std::vector<Query> Q;
  std::vector<std::pair<bool, int>> T;
  void solve(int l, int r, std::vector<Val> &ans){
    if(r - l < 2) return;
    int mid = (l + r) / 2;
    solve(l, mid, ans);
    solve(mid, r, ans);
    std::vector<Idx> Y;
    for(int i = l; i < mid; i++){
      if(!T[i].first) continue;
      int id = T[i].second;
      Y.push_back(U[id].ly);
      Y.push_back(U[id].ry);
    }
    if(Y.size() == 0) return;
    std::sort(Y.begin(), Y.end());
    Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
    std::vector<Event> E;
    for(int i = l; i < mid; i++){
      if(!T[i].first) continue;
      int id = T[i].second;
      int lyc = std::lower_bound(Y.begin(), Y.end(), U[id].ly) - Y.begin();
      int ryc = std::lower_bound(Y.begin(), Y.end(), U[id].ry) - Y.begin();
      E.push_back(Event{U[id].lx, lyc, ryc, U[id].a, U[id].b, U[id].c});
      E.push_back(Event{U[id].rx, lyc, ryc, -U[id].a, -U[id].b, -U[id].c});
    }
    for(int i = mid; i < r; i++){
      if(T[i].first) continue;
      int id = T[i].second;
      int y = std::upper_bound(Y.begin(), Y.end(), Q[id].y) - Y.begin();
      E.push_back(Event{Q[id].x, y, Q[id].id + qlim, 0, Q[id].y, 0});
    }
    std::sort(E.begin(), E.end(), [](const Event &a, const Event &b){
      if(a.x == b.x) return a.ryc < b.ryc;
      return a.x < b.x;
    });
    binary_indexed_tree3<Val> bit(Y.size());
    for(const Event &e : E){
      if(e.ryc < qlim){
        bit.update(e.lyc, e.a, e.b, e.c);
        bit.update(e.ryc, -e.a, -e.b, -e.c);
      }else{
        int id = e.ryc - qlim;
        auto [a, b, c] = bit.query(e.lyc);
        ans[id] += a * e.x + b * e.b + c;
      }
    }
  }
public:
  // lx <= x < rx, ly <= y < ryを満たす全ての(x, y)にax + by + cを足す
  void update(Idx lx, Idx rx, Idx ly, Idx ry, Val a, Val b, Val c){
    T.push_back({1, U.size()});
    U.push_back(Update{lx, rx, ly, ry, a, b, c});
  }
  // get(x, y)
  void query(Idx x, Idx y){
    T.push_back({0, Q.size()});
    Q.push_back(Query{x, y, (int)Q.size()});
  }
  std::vector<Val> solve(){
    std::vector<Val> ans(Q.size(), 0);
    solve(0, T.size(), ans);
    return ans;
  }
};
#endif
