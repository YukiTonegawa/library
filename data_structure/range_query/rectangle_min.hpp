#ifndef _RECTANGLE_MIN_H_
#define _RECTANGLE_MIN_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include "kd_tree.hpp"
#include "../segment_tree/segment_tree2d.hpp"

template<typename T>
struct chmin_struct{
  using Val = T;
  static Val id(){return std::numeric_limits<T>::max() / 2;}
  static bool compare(Val a, Val b){return a < b;}
  static Val update(Val a, Val b){return std::min(a, b);}
  static Val merge(Val a, Val b){return std::min(a, b);}
};
template<typename T>
struct chmax_struct{
  using Val = T;
  static Val id(){return std::numeric_limits<T>::min() / 2;}
  static bool compare(Val a, Val b){return a > b;}
  static Val update(Val a, Val b){return std::max(a, b);}
  static Val merge(Val a, Val b){return std::max(a, b);}
};

// 長方形がN個ありそれぞれ値が割り振られている, 各クエリでは以下の操作を行う
// q(x, y) := (x, y)を含む長方形に書かれている数字のうち最小のものを答える
template<typename Idx, typename compare_struct>
struct offline_static_rectangle_chmin_point_get{
  using Val = typename compare_struct::Val;
private:
  struct Rect{
    Idx lx, rx, ly, ry;
    Val z;
  };
  std::vector<Rect> R;
  kd_tree<Idx, int> kdt;
  using Point = typename kd_tree<Idx, int>::point;
  std::vector<Point> P;
public:
  void update(Idx lx, Idx rx, Idx ly, Idx ry, Val z){
    R.push_back({lx, rx, ly, ry, z});
  }
  void query(Idx x, Idx y){
    int qcnt = P.size();
    P.push_back(Point(x, y, qcnt));
  }
  // 追加される長方形を値の小さい順に見る
  // これがある質問点を含む場合, その質問点を消していい(それ以降の長方形はより値が大きいため)
  std::vector<Val> solve(){
    kdt = kd_tree<Idx, int>(P);
    std::sort(R.begin(), R.end(), [&](const Rect &a, const Rect &b){
      return compare_struct::compare(a.z, b.z);
    });
    std::vector<Val> ans(P.size(), compare_struct::id());
    for(Rect &r : R){
      auto V = kdt.range_find(r.lx, r.rx, r.ly, r.ry);
      for(auto v : V){
        ans[v->p.z] = r.z;
        kdt.erase(v);
      }
    }
    return ans;
  }
};

template<typename Idx, typename compare_struct>
struct offline_dynamic_rectangle_chmin_point_get{
  using Val = typename compare_struct::Val;
private:
  struct Rect{
    Idx lx, rx, ly, ry;
    int t;
    Val z;
  };
  std::vector<Rect> R;
  kd_tree_3d<Idx, int> kdt;
  using Point = typename kd_tree_3d<Idx, int>::point;
  std::vector<Point> P;
public:
  void update(Idx lx, Idx rx, Idx ly, Idx ry, Val z){
    int qcnt = P.size();
    R.push_back({lx, rx, ly, ry, qcnt, z});
  }
  void query(Idx x, Idx y){
    int qcnt = P.size();
    P.push_back(Point(x, y, qcnt, qcnt));
  }
  // 追加される長方形を値の小さい順に見る
  // これがある質問点を含む場合, その質問点を消していい(それ以降の長方形はより値が大きいため)
  std::vector<Val> solve(){
    kdt = kd_tree_3d<Idx, int>(P);
    std::sort(R.begin(), R.end(), [&](const Rect &a, const Rect &b){
      return compare_struct::compare(a.z, b.z);
    });
    std::vector<Val> ans(P.size(), compare_struct::id());
    for(Rect &r : R){
      auto V = kdt.range_find(r.lx, r.rx, r.ly, r.ry, r.t, P.size());
      for(auto v : V){
        ans[v->p.w] = r.z;
        kdt.erase(v);
      }
    }
    return ans;
  }
};


template<typename Idx, typename compare_struct>
struct offline_static_rectangle_min{
  using Val = typename compare_struct::Val;
private:
  kd_tree_rectangle<Idx, int> kdt;
  std::vector<std::tuple<Idx, Idx, Idx, Idx, int>> R;
  std::vector<std::tuple<Idx, Idx, Val>> P;
public:
  void query(Idx lx, Idx rx, Idx ly, Idx ry){
    int qcnt = R.size();
    R.push_back({lx, rx, ly, ry, qcnt});
  }
  void update(Idx x, Idx y, Val z){
    P.push_back({x, y, z});
  }
  std::vector<Val> solve(){
    kdt = kd_tree_rectangle<Idx, int>(R);
    std::sort(P.begin(), P.end(), [&](const std::tuple<Idx, Idx, Val> &a, const std::tuple<Idx, Idx, Val> &b){
      return compare_struct::compare(std::get<2>(a), std::get<2>(b));
    });
    std::vector<Val> ans(R.size(), compare_struct::id());
    for(auto [x, y, z] : P){
      auto V = kdt.fp(x, y);
      for(auto v : V){
        ans[v->p.z] = z;
        kdt.erase(v);
      }
    }
    return ans;
  }
};

template<typename Idx, typename compare_struct>
struct offline_static_upper_rectangle_min{
  using Val = typename compare_struct::Val;
private:
  struct Point{
    Idx x, y;
    Val z;
  };
  struct Query{
    Idx rx, ly, ry;
  };
  std::vector<Point> P;
  std::vector<Query> Q;
public:
  void update(Idx x, Idx y, Val z){P.push_back(Point{x, y, z});}
  // [0, rx) × [ly, ry)のmin
  void query(Idx rx, Idx ly, Idx ry){Q.push_back(Query{rx, ly, ry});}
  
  std::vector<Val> solve(){
    struct Event{
      Idx x;
      int ly, ry;
      int id;
    };
    int N = Q.size();
    std::vector<Val> ans(N, compare_struct::id());
    if(P.empty() || Q.empty()) return ans;
    std::vector<Event> Q2(N);
    std::vector<Idx> Y;
    std::sort(P.begin(), P.end(), [](const Point &a, const Point &b){return a.y < b.y;});
    for(Point &t : P){
      if(Y.empty() || Y.back() != t.y) Y.push_back(t.y);
      t.y = int(Y.size()) - 1;
    }
    for(int i = 0; i < N; i++){
      int ly = std::lower_bound(Y.begin(), Y.end(), Q[i].ly) - Y.begin();
      int ry = std::lower_bound(Y.begin(), Y.end(), Q[i].ry) - Y.begin();
      Q2[i] = Event{Q[i].rx, ly, ry, i};
    }
    std::sort(P.begin(), P.end(), [](const Point &a, const Point &b){return a.x < b.x;});
    std::sort(Q2.begin(), Q2.end(), [](const Event &a, const Event &b){return a.x < b.x;});
    segment_tree<compare_struct> seg(Y.size());
    
    int p = 0, q = 0;
    while(q < N){
      if(p == P.size() || Q2[q].x <= P[p].x){
        Val x = (!p ? compare_struct::id() : seg.query(Q2[q].ly, Q2[q].ry));
        ans[Q2[q].id] = compare_struct::merge(ans[Q2[q].id], x);
        q++;
      }else{
        Val tmp = seg.get(P[p].y);
        seg.set(P[p].y, std::min(tmp, P[p].z));
        p++;
      }
    }
    return ans;
  }
};

template<typename Idx, typename compare_struct>
struct offline_point_min_upper_rectangle_min{
  using Val = typename compare_struct::Val;
private:
  struct query_type{
    int id;
    Idx rx, ly, ry;
    Val z;
    query_type(int a, Idx c, Idx d, Idx e, Val f): id(a), rx(c), ly(d), ry(e), z(f){}
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

    segment_tree<compare_struct> seg(y.size());

    std::vector<event_type> e;
    for(int i = l; i < mid; i++){
      if(q[i].id == -1){
        int y_idx = std::lower_bound(y.begin(), y.end(), q[i].ly) - y.begin();
        e.push_back(event_type(q[i].rx, 1, -1, y_idx, 0, q[i].z));
      }
    }
    for(int i = mid; i < r; i++){
      if(q[i].id != -1){
        int ly_idx = std::lower_bound(y.begin(), y.end(), q[i].ly) - y.begin();
        int ry_idx = std::lower_bound(y.begin(), y.end(), q[i].ry) - y.begin();
        e.push_back(event_type(q[i].rx, 0, q[i].id, ly_idx, ry_idx, 0));
      }
    }
    std::sort(e.begin(), e.end(), [](event_type &a, event_type &b){
      if(a.x != b.x) return a.x < b.x;
      return a.q < b.q;
    });
    bool updated = false;
    for(event_type &ei : e){
      if(ei.q == 0){
        if(updated){
          Val x = seg.query(ei.lyc, ei.ryc);
          ans[ei.id] = compare_struct::merge(ans[ei.id], x);
        }
      }else{
        updated = true;
        seg.update(ei.lyc, ei.z);
      }
    }
  }
public:
  offline_point_min_upper_rectangle_min(){}
  void update(Idx x, Idx y, Val z){
    q.push_back(query_type(-1, x, y, 0, z));
    qcount.push_back(0);
  }
  // [0, rx)　× [ly, ry)
  void query(Idx rx, Idx ly, Idx ry){
    q.push_back(query_type(qs++, rx, ly, ry, 0));
    qcount.push_back(1);
  }
  std::vector<Val> solve(){
    std::vector<Val> res(qs, compare_struct::id());
    for(int i = 1; i < qcount.size(); i++) qcount[i] += qcount[i - 1];
    solve(0, q.size(), res);
    return res;
  }
};
#endif
