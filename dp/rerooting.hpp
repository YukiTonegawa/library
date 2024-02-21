#ifndef _REROOTING_H_
#define _REROOTING_H_
#include <vector>
#include <cassert>

// 全方位木dp

// 注意
// 辺{s, t}がある時, 辺{t, s}が必要
// 実行すると辺の順番が入れ替わる(親への辺が各頂点のリストの最後になる)
// 辺属性と頂点属性を同時に持てるが, 同じ辺について向きによって値を変えることはできない
// 頂点vでマージするとき, vの頂点属性および子部分木の値を好きな順番でマージできる必要がある
template<typename rerooting_structure, typename edge, typename Val>
struct rerooting{
  using Tree = std::vector<std::vector<edge>>;
  // sum_left / sum_right[i][j] := iの子部分木の値(merge_up済)の左 / 右からの累積和
  std::vector<std::vector<Val>> sum_left, sum_right;
  // sum[i] := iの頂点属性と子部分木の値をマージした値
  std::vector<Val> sum;

  void build(Tree &g, const std::vector<Val> &v){
    auto dfs_up = [&](auto &&dfs_up, int cur, int par)->Val{
      int s = g[cur].size();
      sum_left[cur].resize(s, rerooting_structure::template id<Val>());
      Val res = rerooting_structure::template id<Val>();
      int par_idx = -1;
      for(int i = 0; i < s; i++){
        if(g[cur][i].t == par){
          par_idx = i;
          sum_left[cur][i] = rerooting_structure::template id<Val>();
        }else sum_left[cur][i] = rerooting_structure::template merge_up<edge, Val>(g[cur][i], dfs_up(dfs_up, g[cur][i].t, cur));
        res = rerooting_structure::template merge_sibling<Val>(res, sum_left[cur][i]);
      }
      assert(par == -1 || par_idx != -1);
      // 親への辺を末尾にする
      if(par_idx != -1){
        std::swap(sum_left[cur][par_idx], sum_left[cur][s - 1]);
        std::swap(g[cur][par_idx], g[cur][s - 1]);
      }
      return rerooting_structure::template merge_parent<Val>(v[cur], res);
    };
    auto dfs_down = [&](auto &&dfs_down, int cur, int par)->void{
      int s = g[cur].size();
      sum_right[cur].resize(s, rerooting_structure::template id<Val>());
      for(int i = 0; i < s; i++){
        if(!i) sum_right[cur][i] = sum_left[cur][s - 1];
        else sum_right[cur][i] = rerooting_structure::template merge_sibling<Val>(sum_left[cur][s - 1 - i], sum_right[cur][i - 1]);
      }
      for(int i = 1; i < s; i++) sum_left[cur][i] = rerooting_structure::template merge_sibling<Val>(sum_left[cur][i], sum_left[cur][i - 1]);

      for(int i = 0; i < s; i++){
        if(g[cur][i].t == par) continue;
        Val left = (!i ? rerooting_structure::template id<Val>() : sum_left[cur][i - 1]);
        Val right = (i == s - 1 ? rerooting_structure::template id<Val>() : sum_right[cur][s - i - 2]);
        Val mid = rerooting_structure::template merge_parent<Val>(v[cur], rerooting_structure::template merge_sibling<Val>(left, right));
        sum_left[g[cur][i].t].back() = rerooting_structure::template merge_up<edge, Val>(g[cur][i], mid);
        dfs_down(dfs_down, g[cur][i].t, cur);
      }
    };
    int root = 0;
    int n = g.size();
    sum_left.resize(n);
    sum_right.resize(n);
    dfs_up(dfs_up, root, -1);
    dfs_down(dfs_down, root, -1);
    sum = v;
    for(int i = 0; i < n; i++){
      if(!sum_left[i].empty()){
        sum[i] = rerooting_structure::template merge_parent<Val>(sum[i], sum_left[i].back());
      }
    }
  }
  rerooting(){}
  rerooting(Tree &g, const std::vector<Val> &v){
    build(g, v);
  }
};

struct edpc_subtree{
  template<typename Val>
  static Val id(){
    return 1;
  }
  template<typename edge, typename Val>
  static Val merge_up(edge &e, Val a){
    return a + 1;
  }
  template<typename Val>
  static Val merge_sibling(Val a, Val b){
    return a * b;
  }
  template<typename Val>
  static Val merge_parent(Val p, Val c){
    return merge_sibling<Val>(p, c);
  }
};

struct maxdist{
  template<typename Val>
  static Val id(){
    return {0, 0}; // 距離, 頂点番号
  }
  template<typename edge, typename Val>
  static Val merge_up(edge &e, Val a){
    return {a.first + 1, a.second};
  }
  template<typename Val>
  static Val merge_sibling(Val a, Val b){
    return std::max(a, b);
  }
  template<typename Val>
  static Val merge_parent(Val p, Val c){
    return merge_sibling<Val>(p, c);
  }
};
struct mindist_to_marked_vertex{
  static constexpr int inf = 1 << 25;
  template<typename Val>
  static Val id(){
    return inf; // マークされた頂点までの距離(マークされていない頂点はinf)
  }
  template<typename edge, typename Val>
  static Val merge_up(edge &e, Val a){
    return min(a + 1, inf);
  }
  template<typename Val>
  static Val merge_sibling(Val a, Val b){
    return std::min(a, b);
  }
  template<typename Val>
  static Val merge_parent(Val p, Val c){
    return merge_sibling<Val>(p, c);
  }
};
struct tree_path_composite_sum{
  template<typename Val>
  static Val id(){
    return {0, 0};
  }
  template<typename edge, typename Val>
  static Val merge_up(edge &e, Val a){
    return {a.first, e.w.first * a.second + e.w.second * a.first};
  }
  template<typename Val>
  static Val merge_sibling(Val a, Val b){
    return {a.first + b.first, a.second + b.second};
  }
  template<typename Val>
  static Val merge_parent(Val p, Val c){
    return merge_sibling<Val>(p, c);
  }
};

struct expensive_expense{
  template<typename Val>
  static Val id(){
    return 0;
  }
  template<typename edge, typename Val>
  static Val merge_up(edge &e, Val a){
    return a + e.w;
  }
  template<typename Val>
  static Val merge_sibling(Val a, Val b){
    return std::max(a, b);
  }
  template<typename Val>
  static Val merge_parent(Val p, Val c){
    return merge_sibling<Val>(p, c);
  }
};
struct vertex_deletion{
  template<typename Val>
  static Val id(){
    return {0, -1};
  }
  template<typename edge, typename Val>
  static Val merge_up(edge &e, Val a){
    return a;
  }
  template<typename Val>
  static Val merge_sibling(Val a, Val b){
    if(a.second == -1 || b.second == -1) return a.second == -1 ? b : a;
    return {std::max(a.first + b.second, a.second + b.first), a.second + b.second};
  }
  template<typename Val>
  static Val merge_parent(Val p, Val c){
    if(c.second == -1) return p;
    return {std::max(c.first, c.second), c.first + 1};
  }
};

struct cf899d{
  template<typename Val>
  static Val id(){
    return {0, 0}; // コストの和, サイズ
  }
  template<typename edge, typename Val>
  static Val merge_up(edge &e, Val a){
    return {a.first + (long long)e.w * a.second, a.second};
  }
  template<typename Val>
  static Val merge_sibling(Val a, Val b){
    return {a.first + b.first, a.second + b.second};
  }
  template<typename Val>
  static Val merge_parent(Val p, Val c){
    return merge_sibling<Val>(p, c);
  }
};

/*
int main(){
  io_init();
  int n;
  std::cin >> n;
  using e = simple_edge<int>;
  tree<e> t(n);
  range(i, 0, n - 1){
    int a, b;
    std::cin >> a >> b;
    a--, b--;
    t.add_edge(a, {a, b});
    t.add_edge(b, {b, a});
  }
  std::vector<pair<int, int>> v(n, {0, 0});
  rerooting<vertex_deletion, e, std::pair<int, int>> r(t.g, v);
  int max_matching = max(r.sum[0].first, r.sum[0].second);
  int ans = 0;
  range(i, 0, n){
    int tmp = std::max(r.sum_left[i].back().first, r.sum_left[i].back().second);
    if(tmp == max_matching) ans++;
  }
  std::cout << ans << '\n';
}
*/
#endif