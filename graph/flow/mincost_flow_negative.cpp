#ifndef _MINCOST_FLOW_NEGATIVE_H_
#define _MINCOST_FLOW_NEGATIVE_H_

#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>
#include "../../minior/simple_queue.hpp"

/*
始点、終点、フロー固定、負の辺を扱える

flow():
  pair(達成できたフロー、最小コスト)を返す
  指定した流量と異なる場合例外処理する必要がある

正の辺の場合edge.flow流れている
負の辺の場合edge.cap - edge.flow流れている
*/
/*
template<typename Cap, typename Cost>
struct mcf_graph{
public:
  ll negative_flow = 0;//負辺の容量
  ll negative_cost_sum = 0;// < 0
  int s, t, S, T;
  Cap target_flow;
  mcf_graph() {}
  mcf_graph(int n, int s, int t, Cap f): s(s), t(t), S(n), T(n+1), target_flow(f), _n(n+2), g(n+2){
    add_edge(S, s, target_flow, 0);
    add_edge(t, T, target_flow, 0);
  } //s:n, t:n+1

  int add_edge(int from, int to, Cap cap, Cost cost) {
    assert(0 <= from && from < _n);
    assert(0 <= to && to < _n);

    if(cost<0){
      negative_flow += cap;
      negative_cost_sum += cap * cost;
      add_edge(S, to, cap, 0);
      add_edge(to, from, cap, -cost);
      add_edge(from, T, cap, 0);
      return 0;
    }

    int m = int(pos.size());
    pos.push_back({from, int(g[from].size())});
    g[from].push_back(_edge{to, int(g[to].size()), cap, cost});
    g[to].push_back(_edge{from, int(g[from].size()) - 1, 0, -cost});
    return m;
  }

  struct edge {
    int from, to;
    Cap cap, flow;
    Cost cost;
  };

  edge get_edge(int i) {
    int m = int(pos.size());
    assert(0 <= i && i < m);
    auto _e = g[pos[i].first][pos[i].second];
    auto _re = g[_e.to][_e.rev];
    return edge{pos[i].first, _e.to, _e.cap + _re.cap, _re.cap, _e.cost,};
  }

  vector<edge> edges() {
    int m = int(pos.size());
    vector<edge> result(m);
    for(int i=0;i<m;i++) result[i] = get_edge(i);
    return result;
  }

  pair<Cap, Cost> flow() {
    pair<Cap, Cost> ret = slope().back();
    ret.first -= negative_flow;
    ret.second += negative_cost_sum;
    return ret;
  }

  vector<pair<Cap, Cost>> slope() {
    assert(0 <= S && S < _n);
    assert(0 <= T && T < _n);
    assert(S != T);
    // variants (C = maxcost):
    // -(n-1)C <= dual[s] <= dual[i] <= dual[t] = 0
    // reduced cost (= e.cost + dual[e.from] - dual[e.to]) >= 0 for all edge
    vector<Cost> dual(_n, 0), dist(_n);
    vector<int> pv(_n), pe(_n);
    vector<bool> vis(_n);
    auto dual_ref = [&]() {
      fill(dist.begin(), dist.end(), numeric_limits<Cost>::max());
      fill(pv.begin(), pv.end(), -1);
      fill(pe.begin(), pe.end(), -1);
      fill(vis.begin(), vis.end(), false);
      struct Q {
        Cost key;
        int to;
        bool operator<(Q r) const { return key > r.key; }
      };
      priority_queue<Q> que;
      dist[S] = 0;
      que.push(Q{0, S});
      while (!que.empty()) {
        int v = que.top().to;
        que.pop();
        if (vis[v]) continue;
        vis[v] = true;
        if (v == T) break;
        // dist[v] = shortest(s, v) + dual[s] - dual[v]
        // dist[v] >= 0 (all reduced cost are positive)
        // dist[v] <= (n-1)C
        for(int i=0;i<int(g[v].size());i++) {
          auto e = g[v][i];
          if (vis[e.to] || !e.cap) continue;
          // |-dual[e.to] + dual[v]| <= (n-1)C
          // cost <= C - -(n-1)C + 0 = nC
          Cost cost = e.cost - dual[e.to] + dual[v];
          if(dist[e.to] - dist[v] > cost) {
            dist[e.to] = dist[v] + cost;
            pv[e.to] = v;
            pe[e.to] = i;
            que.push(Q{dist[e.to], e.to});
          }
        }
      }
      if(!vis[T]) return false;
      for(int v=0;v<_n;v++) {
        if(!vis[v]) continue;
        // dual[v] = dual[v] - dist[t] + dist[v]
        //         = dual[v] - (shortest(s, t) + dual[s] - dual[t]) + (shortest(s, v) + dual[s] - dual[v])
        //         = - shortest(s, t) + dual[t] + shortest(s, v)
        //         = shortest(s, v) - shortest(s, t) >= 0 - (n-1)C
        dual[v] -= dist[T] - dist[v];
      }
      return true;
    };
    Cap flow = 0;
    Cost cost = 0, prev_cost = -1;
    vector<pair<Cap, Cost>> result;
    result.push_back({flow, cost});
    while(flow < target_flow + negative_flow) {
      if(!dual_ref()) break;
      Cap c = target_flow + negative_flow - flow;
      for(int v=T;v!=S;v=pv[v]) c = min(c, g[pv[v]][pe[v]].cap);
      for(int v=T;v!=S;v=pv[v]) {
        auto& e = g[pv[v]][pe[v]];
        e.cap -= c;
        g[v][e.rev].cap += c;
      }
      Cost d = -dual[S];
      flow += c;
      cost += c * d;
      if(prev_cost == d) result.pop_back();
      result.push_back({flow, cost});
      prev_cost = cost;
    }
    return result;
  }
private:
  int _n;
  struct _edge {
    int to, rev;
    Cap cap;
    Cost cost;
  };
  vector<pair<int, int>> pos;
  vector<vector<_edge>> g;
};
*/
#endif