#ifndef _RECTANGLE_SUM_H_
#define _RECTANGLE_SUM_H_
#include <vector>
#include <tuple>
#include <array>
#include <algorithm>
#include "../segment_tree/binary_indexed_tree.hpp"
#include "../segment_tree/binary_indexed_tree_compressed.hpp"
#include "../../minior/map_sum_cache.hpp"

template<typename Idx, typename Val>
struct offline_static_rectangle_sum{
private:
  struct Point{
    Idx x, y;
    Val z;
  };
  struct Query{
    Idx lx, rx, ly, ry;
  };
  std::vector<Point> P;
  std::vector<Query> Q;
public:
  void update(Idx x, Idx y, Val z){P.push_back(Point{x, y, z});}
  void query(Idx lx, Idx rx, Idx ly, Idx ry){Q.push_back(Query{lx, rx, ly, ry});}
  std::vector<Val> solve(){
    struct Event{
      Idx x;
      int ly, ry;
      int id;
    };
    int N = Q.size();
    if(P.empty() || Q.empty()) return std::vector<Val>(N, 0);
    std::vector<Event> Q2(2 * N);
    std::vector<Val> ans(N, 0);
    std::vector<Idx> Y;
    std::sort(P.begin(), P.end(), [](const Point &a, const Point &b){return a.y < b.y;});
    for(Point &t : P){
      if(Y.empty() || Y.back() != t.y) Y.push_back(t.y);
      t.y = int(Y.size()) - 1;
    }
    for(int i = 0; i < N; i++){
      int ly = std::lower_bound(Y.begin(), Y.end(), Q[i].ly) - Y.begin();
      int ry = std::lower_bound(Y.begin(), Y.end(), Q[i].ry) - Y.begin();
      Q2[2 * i] = Event{Q[i].lx, ly, ry, i};
      Q2[2 * i + 1] = Event{Q[i].rx, ly, ry, i + N};
    }
    std::sort(P.begin(), P.end(), [](const Point &a, const Point &b){return a.x < b.x;});
    std::sort(Q2.begin(), Q2.end(), [](const Event &a, const Event &b){return a.x < b.x;});
    binary_indexed_tree<Val> bit(Y.size());
    int p = 0, q = 0;
    while(q < 2 * N){
      if(p == P.size() || Q2[q].x <= P[p].x){
        Val sum = (!p ? Val(0) : bit.query(Q2[q].ly, Q2[q].ry));
        if(Q2[q].id >= N) ans[Q2[q].id - N] += sum;
        else ans[Q2[q].id] -= sum;
        q++;
      }else{
        bit.update(P[p].y, P[p].z);
        p++;
      }
    }
    return ans;
  }
};

template<typename Idx, typename Val>
struct offline_point_add_rectangle_sum{
private:
  struct query_type{
    int id;
    Idx lx, rx, ly, ry;
    Val z;
    query_type(int a, Idx b, Idx c, Idx d, Idx e, Val f): id(a), lx(b), rx(c), ly(d), ry(e), z(f){}
  };
  struct event_type{
    Idx x;
    int q, id;
    int lyc, ryc;
    Val z;
    event_type(Idx a, int b, int c, int d, int e, Val f): x(a), q(b), id(c), lyc(d), ryc(e), z(f){}
  };
  std::vector<query_type> q;
  std::vector<int> qcount{0};
  int qs = 0;
  void solve(int l, int r, std::vector<Val> &ans){
    if(r - l < 2) return;
    int mid = (l + r) >> 1;
    solve(l, mid, ans);
    solve(mid, r, ans);
    int left_update = (mid - l) - (qcount[mid] - qcount[l]);
    int right_query= qcount[r] - qcount[mid];
    if(left_update == 0 || right_query == 0) return;

    // compress y
    std::vector<Idx> y;
    for(int i = l; i < mid; i++) if(q[i].id == -1) y.push_back(q[i].ly);
    std::sort(y.begin(), y.end());
    y.erase(std::unique(y.begin(), y.end()), y.end());

    binary_indexed_tree<Val> bit(y.size());
    std::vector<event_type> e;
    for(int i = l; i < mid; i++){
      if(q[i].id == -1){
        int y_idx = std::lower_bound(y.begin(), y.end(), q[i].ly) - y.begin();
        e.push_back(event_type(q[i].lx, 2, -1, y_idx, 0, q[i].z));
      }
    }
    for(int i = mid; i < r; i++){
      if(q[i].id != -1){
        int ly_idx = std::lower_bound(y.begin(), y.end(), q[i].ly) - y.begin();
        int ry_idx = std::lower_bound(y.begin(), y.end(), q[i].ry) - y.begin();
        e.push_back(event_type(q[i].lx, 0, q[i].id, ly_idx, ry_idx, 0));
        e.push_back(event_type(q[i].rx, 1, q[i].id, ly_idx, ry_idx, 0));
      }
    }
    std::sort(e.begin(), e.end(), [](event_type &a, event_type &b){
      if(a.x != b.x) return a.x < b.x;
      return a.q < b.q;
    });
    bool updated = false;
    for(event_type &ei : e){
      if(ei.q == 0){
        if(updated) ans[ei.id] -= bit.query(ei.lyc, ei.ryc);
      }else if(ei.q == 1){
        if(updated) ans[ei.id] += bit.query(ei.lyc, ei.ryc);
      }else{
        updated = true;
        bit.update(ei.lyc, ei.z);
      }
    }
  }
public:
  offline_point_add_rectangle_sum(){}
  void update(Idx x, Idx y, Val z){
    q.push_back(query_type(-1, x, 0, y, 0, z));
    qcount.push_back(0);
  }
  void query(Idx lx, Idx rx, Idx ly, Idx ry){
    q.push_back(query_type(qs++, lx, rx, ly, ry, 0));
    qcount.push_back(1);
  }
  std::vector<Val> solve(){
    std::vector<Val> ret(qs, 0);
    for(int i = 1; i < qcount.size(); i++) qcount[i] += qcount[i - 1];
    solve(0, q.size(), ret);
    return ret;
  }
};

template<typename Idx, typename Val>
struct partial_offline_point_add_rectangle_sum{
private:
  bool is_built = false;
  int M = 1;
  std::vector<Idx> X;
  std::vector<std::vector<Idx>> Y;
  std::vector<binary_indexed_tree<Val>> BITs;
  using point = std::tuple<Idx, Idx, Val>;
  int lbx(Idx x){
    return std::lower_bound(X.begin(), X.end(), x) - X.begin();
  }
public:
  std::vector<point> P;
  void set_initial_point(Idx x, Idx y, Val z){
    assert(!is_built);
    P.push_back({x, y, z});
  }
  void build(){
    assert(!is_built);
    is_built = true;
    int n = P.size();
    std::sort(P.begin(), P.end(), [](const point &a, const point &b){
      return std::get<1>(a) < std::get<1>(b);
    });
    for(int i = 0; i < n; i++) X.push_back(std::get<0>(P[i]));
    std::sort(X.begin(), X.end());
    X.erase(std::unique(X.begin(), X.end()), X.end());
    M = X.size();
    std::vector<std::vector<Val>> tmp(M + 1);
    BITs.resize(M + 1);
    Y.resize(M + 1);
    for(int i = 0; i < n; i++){
      auto [x, y, z] = P[i];
      int k = lbx(x);
      for(int j = k + 1; j <= M; j += (j & (-j))){
        if(Y[j].empty()||Y[j].back()!=y){
          Y[j].push_back(y);
          tmp[j].push_back(z);
        }else{
          tmp[j].back() += z;
        }
      }
    }
    for(int i = 0; i <= M; i++) BITs[i] = binary_indexed_tree<Val>(tmp[i]);
  }
  partial_offline_point_add_rectangle_sum(){}
  partial_offline_point_add_rectangle_sum(const std::vector<point> &v){
    P = v;
    build();
  }
  void update(Idx x, Idx y, Val z){
    assert(is_built);
    int xidx = lbx(x);
    for(int i = xidx + 1; i <= M; i += (i & (-i))){
      auto yidx = std::lower_bound(Y[i].begin(), Y[i].end(), y) - Y[i].begin();
      BITs[i].update(yidx, z);
    }
  }
  Val query(int rx, Idx ly, Idx ry){
    assert(is_built);
    Val ret = 0;
    for(int i = rx; i > 0; i -= (i & (-i))){
      int ridx = std::lower_bound(Y[i].begin(), Y[i].end(), ry) - Y[i].begin();
      int lidx = std::lower_bound(Y[i].begin(), Y[i].end(), ly) - Y[i].begin();
      ret += BITs[i].query(lidx, ridx);
    }
    return ret;
  }
  Val query(Idx lx, Idx rx, Idx ly, Idx ry){
    return query(lbx(rx), ly, ry) - query(lbx(lx), ly, ry);
  }
};

template<typename Idx, typename Val, int max_n_log = 30>
struct online_point_add_rectangle_sum{
private:
  struct node{
    map_sum_cache<Idx, Val> mp;
    std::vector<std::pair<Idx, Val>> initial_points; // 初期化用　
    node *l, *r;
    node(): l(nullptr), r(nullptr){}
  };
  node *root;
  void update_inner(Idx x, Idx y, Val z){
    node *v = root;
    Idx lx = 0, rx = (Idx)1 << max_n_log;
    while(true){
      Idx mx = (lx + rx) / 2;
      if(rx <= x){
        if(!v->r) v->r = new node();
        v = v->r;
        lx = rx, rx += rx - mx;
      }else{
        v->mp.update(y, z);
        if(rx - 1 == x) return;
        rx = mx;
        if(!v->l) v->l = new node();
        v = v->l;
      }
    }
  }
  Val query_inner(Idx x, Idx ly, Idx ry){
    Idx lx = 0, rx = (Idx)1 << max_n_log;
    Val ret = 0;
    node *v = root;
    while(v){
      Idx mx = (lx + rx) / 2;
      if(rx <= x){
        ret += v->mp.query(ly, ry);
        if(rx == x) return ret;
        v = v->r;
        lx = rx;
        rx += rx - mx;
      }else{
        v = v->l;
        rx = mx;
      }
    }
    return ret;
  }
  using point = std::tuple<Idx, Idx, Val>;
public:
  online_point_add_rectangle_sum(): root(new node()){}
  online_point_add_rectangle_sum(std::vector<point> v): root(new node()){
    sort(v.begin(), v.end(), [](const point &a, const point &b){
      return std::get<1>(a) < std::get<1>(b);
    });
    auto push = [&](Idx x, Idx y, Val z){
      node *v = root;
      Idx lx = 0, rx = (Idx)1 << max_n_log;
      while(true){
        Idx mx = (lx + rx) / 2;
        if(rx <= x){
          if(!v->r) v->r = new node();
          v = v->r;
          lx = rx, rx += rx - mx;
        }else{
          if(v->initial_points.empty() || v->initial_points.back().first != y){
            v->initial_points.push_back({y, z});
          }else{
            v->initial_points.back().second += z;
          }
          if(rx - 1 == x) return;
          rx = mx;
          if(!v->l) v->l = new node();
          v = v->l;
        }
      }
    };
    for(auto [x, y, z] : v) push(x, y, z);
    auto init = [&](auto &&init, node *v) -> void {
      v->mp.init_sorted(v->initial_points);
      v->initial_points.clear();
      if(v->l) init(init, v->l);
      if(v->r) init(init, v->r);
    };
    init(init, root);
  }
  void update(Idx x, Idx y, Val z){
    update_inner(x, y, z);
  }
  Val query(Idx lx, Idx rx, Idx ly, Idx ry){
    return query_inner(rx, ly, ry) - query_inner(lx, ly, ry);
  }
};

template<typename Idx, typename Val>
struct offline_static_rectangle_add_rectangle_sum{
private:
  struct bit4{
    int N;
    std::vector<std::array<Val, 4>> sum;
    bit4(int N): N(N), sum(N + 1, {0, 0, 0, 0}){}
    void update(Idx l, Idx r, int lc, int rc, Val z1, Val z2){
      Val a = l * z1, b = r * z1, c = l * z2, d = r * z2;
      for(int i = lc + 1; i <= N; i += (i & (-i))){
        sum[i][0] -= a;
        sum[i][1] += z1;
        sum[i][2] -= c;
        sum[i][3] += z2;
      }
      for(int i = rc + 1; i <= N; i += (i & (-i))){
        sum[i][0] += b;
        sum[i][1] -= z1;
        sum[i][2] += d;
        sum[i][3] -= z2;
      }
    }
    std::pair<Val, Val> query(Idx r, int rc){
      Val a = 0, b = 0, c = 0, d = 0;
      for(int i = rc; i > 0; i -= (i & (-i))){
        a += sum[i][0];
        b += sum[i][1];
        c += sum[i][2];
        d += sum[i][3];
      }
      return {a + (b * r), c + (d * r)};
    }
    std::pair<Val, Val> query(Idx l, Idx r, int lc, int rc){
      auto [cr, dxr] = query(r, rc);
      auto [cl, dxl] = query(l, lc);
      return {cr - cl, dxr - dxl};
    }
  };
  struct Update{
    Idx lx, rx, ly, ry;
    Val z;
    int lyc, ryc;
    Update(Idx lx, Idx rx, Idx ly, Idx ry, Val z, int lyc = 0, int ryc = 0): lx(lx), rx(rx), ly(ly), ry(ry), z(z), lyc(lyc), ryc(ryc){}
  };
  struct Query{
    Idx lx, rx, ly, ry;
    int id, lyc, ryc;
    Query(Idx lx, Idx rx, Idx ly, Idx ry, int id, int lyc = 0, int ryc = 0): lx(lx), rx(rx), ly(ly), ry(ry), id(id), lyc(lyc), ryc(ryc){}
  };
  std::vector<Query> Q;
  std::vector<Update> U;

  void solve(std::vector<Val> &ans){
    int N = U.size(), M = Q.size();
    std::vector<Idx> Y;
    for(int i = 0; i < N; i++){
      Y.push_back(U[i].ly);
      Y.push_back(U[i].ry);
    }
    std::sort(Y.begin(), Y.end());
    Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
    auto lb = [&](Idx y) -> int { return std::lower_bound(Y.begin(), Y.end(), y) - Y.begin();};
    for(int i = 0; i < N; i++){      
      int lyc = lb(U[i].ly), ryc = lb(U[i].ry);
      U[i].lyc = lyc, U[i].ryc = ryc;
      U.push_back(Update(U[i].rx, 0, U[i].ly, U[i].ry, -U[i].z, lyc, ryc));
    }
    for(int i = 0; i < M; i++){
      int lyc = lb(Q[i].ly), ryc = lb(Q[i].ry);
      Q[i].lyc = lyc, Q[i].ryc = ryc;
      Q.push_back(Query(Q[i].rx, 0, Q[i].ly, Q[i].ry, Q[i].id + M, lyc, ryc));
    }
    std::sort(U.begin(), U.end(), [](const Update &a, const Update &b){return a.lx < b.lx;});
    std::sort(Q.begin(), Q.end(), [](const Query &a, const Query &b){return a.lx < b.lx;});
    assert(U.size() == 2 * N && Q.size() == 2 * M);
    bit4 bit(Y.size());
    int uid = 0, qid = 0;
    while(qid < 2 * M){
      if(uid < 2 * N && U[uid].lx < Q[qid].lx){
        bit.update(U[uid].ly, U[uid].ry, U[uid].lyc, U[uid].ryc, -U[uid].z * U[uid].lx, U[uid].z);
        uid++;
      }else{
        auto [a, b] = bit.query(Q[qid].ly, Q[qid].ry, Q[qid].lyc, Q[qid].ryc);
        int id = Q[qid].id;
        if(id >= M){
          ans[id - M] += a + Q[qid].lx * b;
        }else{
          ans[id    ] -= a + Q[qid].lx * b;
        }
        qid++;
      }
    }
  }
public:
  offline_static_rectangle_add_rectangle_sum(){}
  void update(Idx lx, Idx rx, Idx ly, Idx ry, Val z){
    U.push_back(Update(lx, rx, ly, ry, z));
  }
  void query(Idx lx, Idx rx, Idx ly, Idx ry){
    Q.push_back(Query(lx, rx, ly, ry, Q.size()));
  }
  std::vector<Val> solve(){
    std::vector<Val> ans(Q.size(), 0);
    solve(ans);
    return ans;
  }
};

template<typename Val>
struct offline_dynamic_rectangle_add_rectangle_sum{
private:
  struct query_type{
    int lx, rx, ly, ry;
    Val z;
    int type, lyc, ryc;
    query_type(int _lx, int _rx, int _ly, int _ry, Val _z, int _type, int _lyc = 0, int _ryc = 0):
    lx(_lx), rx(_rx), ly(_ly), ry(_ry), z(_z), type(_type), lyc(_lyc), ryc(_ryc){}
  };
  int q = 0;
  std::vector<int> qcnt{0};
  std::vector<query_type> Q;
  void solve(int l, int r, std::vector<Val> &ans){
    if(r - l < 2) return;
    int mid = (l + r) >> 1;
    int left_add = (mid - l) - (qcnt[mid] - qcnt[l]);
    int right_query = qcnt[r] - qcnt[mid];
    if(left_add) solve(l, mid, ans);
    if(right_query) solve(mid, r, ans);
    if(!left_add || !right_query) return;
    std::vector<query_type> Q_tmp;
    // do naive
    if(left_add <= 6 || right_query <= 6){
      for(int i = l; i < mid; i++) if(Q[i].type == -1) Q_tmp.push_back(Q[i]);
      for(int i = mid; i < r; i++){
        if(Q[i].type == -1) continue;
        for(query_type qi: Q_tmp){
          int lx = std::max(Q[i].lx, qi.lx);
          int rx = std::min(Q[i].rx, qi.rx);
          if(lx >= rx) continue;
          int ly = std::max(Q[i].ly, qi.ly);
          int ry = std::min(Q[i].ry, qi.ry);
          if(ly >= ry) continue;
          ans[Q[i].type] += qi.z * Val(rx - lx) * Val(ry - ly);
        }
      }
      return;
    }
    std::vector<int> Y;
    for(int i = l; i < mid; i++){
      if(Q[i].type != -1) continue;
      Y.push_back(Q[i].ly);
      Y.push_back(Q[i].ry);
    }
    for(int i = mid; i < r; i++){
      if(Q[i].type == -1) continue;
      Y.push_back(Q[i].ly);
      Y.push_back(Q[i].ry);
    }
    std::sort(Y.begin(), Y.end());
    Y.erase(std::unique(Y.begin(), Y.end()), Y.end());

    auto lb = [&](int y) -> int {
      return std::lower_bound(Y.begin(), Y.end(), y) - Y.begin();
    };
    for(int i = l; i < mid; i++){
      if(Q[i].type != -1) continue;
      query_type qi = Q[i];
      int lyc = lb(qi.ly), ryc = lb(qi.ry);
      Q_tmp.push_back(query_type(qi.lx, 0, qi.ly, qi.ry, qi.z, -1, lyc, ryc));
      Q_tmp.push_back(query_type(qi.rx, 0, qi.ly, qi.ry, qi.z, -2, lyc, ryc));
    }
    for(int i = mid; i < r; i++){
      if(Q[i].type == -1) continue;
      query_type qi = Q[i];
      int lyc = lb(qi.ly), ryc = lb(qi.ry);
      Q_tmp.push_back(query_type(qi.lx, 0, qi.ly, qi.ry, qi.z, qi.type, lyc, ryc));
      Q_tmp.push_back(query_type(qi.rx, 0, qi.ly, qi.ry, qi.z, qi.type + q, lyc, ryc));
    }
    std::sort(Q_tmp.begin(), Q_tmp.end(), [](const query_type &a, const query_type &b){
      if(a.lx == b.lx) return a.type > b.type;
      return a.lx < b.lx;
    });
    binary_indexed_tree_compressed<int, Val> slope(Y.size()), intercept(Y.size());
    for(query_type &qi: Q_tmp){
      if(qi.type == -1){
        // 傾き +z, x分切片を減らす
        slope.update(qi.ly, qi.ry, qi.lyc, qi.ryc, qi.z);
        intercept.update(qi.ly, qi.ry, qi.lyc, qi.ryc, Val(qi.z) * Val(-qi.lx));
      }else if(qi.type == -2){
        // 傾き -z, x分切片を増やす
        slope.update(qi.ly, qi.ry, qi.lyc, qi.ryc, Val(-qi.z));
        intercept.update(qi.ly, qi.ry, qi.lyc, qi.ryc, Val(qi.z) * Val(qi.lx));
      }else if(qi.type < q){
        // [0, lx) × [ly, ry)を減らす
        Val a = slope.query(qi.ly, qi.ry, qi.lyc, qi.ryc);
        Val b = intercept.query(qi.ly, qi.ry, qi.lyc, qi.ryc);
        ans[qi.type] -= Val(a) * Val(qi.lx) + Val(b);
      }else{
        // [0, rx) × [ly, ry)を増やす
        qi.type -= q;
        Val a = slope.query(qi.ly, qi.ry, qi.lyc, qi.ryc);
        Val b = intercept.query(qi.ly, qi.ry, qi.lyc, qi.ryc);
        ans[qi.type] += Val(a) * Val(qi.lx) + Val(b);
      }
    }
  }
public:
  offline_dynamic_rectangle_add_rectangle_sum(){}
  void update(int lx, int rx, int ly, int ry, Val z){
    Q.push_back(query_type(lx, rx, ly, ry, z, -1));
    qcnt.push_back(0);
  }
  void query(int lx, int rx, int ly, int ry){
    Q.push_back(query_type(lx, rx, ly, ry, 0, q++));
    qcnt.push_back(1);
  }
  std::vector<Val> solve(){
    std::vector<Val> ans(q, 0);
    for(int i = 1; i < qcnt.size(); i++) qcnt[i] += qcnt[i - 1];
    solve(0, Q.size(), ans);
    return ans;
  }
};

template<typename Idx, typename Val>
struct offline_static_cuboid_sum{
private:
  offline_point_add_rectangle_sum<Idx, Val> rect;
  std::vector<std::tuple<Idx, Idx, Idx, Val>> P;
  struct qst{
    Idx rx, ly, ry, lz, rz;
    int id;
    bool add;
  };
  std::vector<qst> Q;
  std::vector<Val> __solve(){
    std::sort(P.begin(), P.end());
    std::sort(Q.begin(), Q.end(), [&](qst &a, qst &b){return a.rx < b.rx;});
    int j = 0;
    for(auto &q : Q){
      while(j < P.size() && std::get<0>(P[j]) < q.rx){
        auto [x, y, z, w] = P[j++];
        rect.update(y, z, w);
      }
      rect.query(q.ly, q.ry, q.lz, q.rz);
    }
    std::vector<Val> ans(Q.size() / 2);
    auto s = rect.solve();
    j = 0;
    for(auto &q : Q){
      if(q.add) ans[q.id] += s[j++];
      else ans[q.id] -= s[j++];
    }
    return ans;
  }
public:
  offline_static_cuboid_sum(){}
  // (x, y, z)に重みwを追加
  void update(Idx x, Idx y, Idx z, Val w){
    P.push_back({x, y, z, w});
  }
  // [lx, rx) × [ly, ry) × [lz, rz)の点の重みの和
  void query(Idx lx, Idx rx, Idx ly, Idx ry, Idx lz, Idx rz){
    int id = Q.size() / 2;
    Q.push_back({rx, ly, ry, lz, rz, id, true});
    Q.push_back({lx, ly, ry, lz, rz, id, false});
  }
  // O((N + Q)log^2(N + Q))
  std::vector<Val> solve(){
    return __solve();
  }
};
template<typename Idx, typename Val>
struct offline_rectangle_add_point_get{
private:
  static constexpr int qlim = 1e8;
  struct Query{
    Idx x, y;
    int id;
  };
  struct Update{
    Idx lx, rx, ly, ry;
    Val z;
  };
  struct Event{
    Idx x;
    int lyc, ryc;
    Val z;
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
    for(int i = mid; i < r; i++){
      if(T[i].first) continue;
      int id = T[i].second;
      Y.push_back(Q[id].y);
    }
    if(Y.empty()) return;
    std::sort(Y.begin(), Y.end());
    Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
    std::vector<Event> E;
    for(int i = l; i < mid; i++){
      if(!T[i].first) continue;
      int id = T[i].second;
      int lyc = std::lower_bound(Y.begin(), Y.end(), U[id].ly) - Y.begin();
      int ryc = std::lower_bound(Y.begin(), Y.end(), U[id].ry) - Y.begin();
      E.push_back(Event{U[id].lx, lyc, ryc, U[id].z});
      E.push_back(Event{U[id].rx, lyc, ryc, -U[id].z});
    }
    for(int i = mid; i < r; i++){
      if(T[i].first) continue;
      int id = T[i].second;
      int y = std::lower_bound(Y.begin(), Y.end(), Q[id].y) - Y.begin();
      E.push_back(Event{Q[id].x, y, Q[id].id + qlim, 0});
    }
    std::sort(E.begin(), E.end(), [](const Event &a, const Event &b){
      if(a.x == b.x) return a.ryc < b.ryc;
      return a.x < b.x;
    });
    binary_indexed_tree<Val> bit(Y.size());
    for(const Event &e : E){
      if(e.ryc < qlim){
        if(e.lyc < Y.size()) bit.update(e.lyc, e.z);
        if(e.ryc < Y.size()) bit.update(e.ryc, -e.z);
      }else{
        int id = e.ryc - qlim;
        ans[id] += bit.query(e.lyc + 1);
      }
    }
  }
public:
  // [lx, rx) × [ly, ry)にzを足す
  void update(Idx lx, Idx rx, Idx ly, Idx ry, Val z){
    T.push_back({1, U.size()});
    U.push_back(Update{lx, rx, ly, ry, z});
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
