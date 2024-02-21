#ifndef _GRAPH_ALGORITHM_EXTRA_H_
#define _GRAPH_ALGORITHM_EXTRA_H_
#include <vector>
#include <set>
#include <algorithm>
#include "../misc/random_number.hpp"
#include "../data_structure/bit_sequence/dynamic_bitset.hpp"

template<typename T>
using vec = std::vector<T>;

// i-bit目が1 -> 頂点iを使う
long long maximum_independent_set(const vec<long long> &g2, long long rem = -1){
  int n = g2.size();
  if(rem == -1) rem = (1LL << n) - 1;
  long long ret = 0;
  int k = -1, m = -1;
  while(true){
    bool update = false;
    for(int i = 0; i < n; i++){
      if(!((rem >> i) & 1)) continue;
      int s = __builtin_popcountll(rem & g2[i]); //次数
      if(s > m) k = i, m = s;
      if(s <= 1){
        rem &= ~(g2[i] | (1LL << i));
        ret |= (1LL << i), update = true;
      }
    }
    if(!update) break;
    k = -1, m = -1;
  }
  if(rem > 0){
    rem &= ~(1LL << k);
    long long p = maximum_independent_set(g2, rem); //kを使わない
    long long q = maximum_independent_set(g2, rem & ~g2[k]); //kを使う
    if(__builtin_popcountll(p) > __builtin_popcountll(q)) ret |= p;
    else ret |= ((1LL << k) | q);
  }
  return ret;
}

int chromatic_number(const vec<long long> &g){
  using int_type = unsigned int;
  using mat = std::vector<int_type>;
  int n = g.size();
  std::vector<int_type> is_ok(1 << n, 0);
  for(int i = 0; i < (1 << n); i++){
    bool f = true;
    for(int j = 0; j < n; j++){
      if(!((i >> j) & 1)) continue;
      if(i & g[j]){
        f = false;
        break;
      }
    }
    is_ok[i] = f;
  }
  auto subset_zeta = [](int n, const std::vector<int_type> &A){
    assert(A.size() == (1 << n));
    int m = 1 << n;
    mat C(m);
    for(int i = 0; i < m; i++) C[i] = A[i];
    for(int i = 1; i < m; i <<= 1){
      for(int j = 0; j < m; j++){
        if(!(j & i)){
          C[j | i] += C[j];
        }
      }
    }
    return C;
  };
  auto multiply = [](int n, mat &C, const mat &D){
    int m = 1 << n;
    for(int i = 0; i < m; i++){
      C[i] *= D[i];
    }
  };
  // メビウス変換をした後に(1 << n) - 1が1以上になるか
  auto check = [](int n, const mat &C){
    int m = 1 << n;
    auto D = C;
    for(int i = 1; i < m; i <<= 1){
      for(int j = 0; j < m; j++){
        if(!(j & i)){
          D[j | i] -= D[j];
        }
      }
    }
    return D[(1 << n) - 1];
  };
  auto x = subset_zeta(n, is_ok);
  decltype(x) y;
  for(int i = 0; i < n; i++){
    if(i == 0) y = x;
    else multiply(n, y, x);
    if(check(n, y)) return i + 1;
  }
  return -1;
}
/*
void enumerate_cliques(){
  int n, m;
  std::cin >> n >> m;
  std::vector<long long> x(n);
  constexpr static int mod = 998244353;
  for(int i = 0; i < n; i++){
    int t;
    std::cin >> t;
    x[i] = t;
  }
  std::vector<__uint128_t> g(n, 0);
  __uint128_t is_alive = ((__uint128_t)1 << n) - 1;

  auto pop_count_128 = [](__uint128_t x)->int{
    constexpr static uint64_t mask_split = ((__uint128_t)1 << 64) - 1;
    return __builtin_popcountll(x & mask_split) + __builtin_popcountll(x >> 64);
  };

  for(int i = 0; i < m; i++){
    int a, b;
    std::cin >> a >> b;
    g[a] |= (__uint128_t)1 << b;
    g[b] |= (__uint128_t)1 << a;
  }

  auto dfs = [&](auto &&dfs, int k, long long mul, __uint128_t s, __uint128_t s2, std::vector<int> &V, int must)->long long{
    if(k == V.size()){
      if((s & s2) == s2) return mul;
      else return 0;
    }
    int cur = V[k];
    __uint128_t cur_bit = (__uint128_t)1 << cur;
    long long res = 0;

    // v[i]を使う
    res += dfs(dfs, k + 1, (mul * x[cur]) % mod, s & (g[cur] | cur_bit), s2 | cur_bit, V, must);
    // v[i]を使わない
    if(cur != must) res += dfs(dfs, k + 1, mul, s, s2, V, must);
    return res >= mod ? res - mod : res;
  };

  long long ans = 0;
  for(int iter = 0; iter < n; iter++){
    int upper = 1;
    while(upper * upper <= 2 * m) upper++;

    int lower_vertex = -1;
    for(int i = 0; i < n; i++){
      if(!((is_alive >> i) & 1)) continue;
      if(pop_count_128(g[i]) < upper){
        lower_vertex = i;
        break;
      }
    }
    // n <= √2m
    if(lower_vertex == -1){
      std::vector<int> v;
      for(int i = 0; i < n; i++) if((is_alive >> i) & 1) v.push_back(i);
      ans += dfs(dfs, 0, 1, ((__uint128_t)1 << n) - 1, 0, v, -1);
      ans--;
      if(ans < 0) ans += mod;
      if(ans >= mod) ans -= mod;
      break;
    }else{
      std::vector<int> v{lower_vertex};
      for(int i = 0; i < n; i++){
        if(((g[lower_vertex] >> i) & 1) && i != lower_vertex && (is_alive >> i) & 1){
          v.push_back(i);
        }
      }
      ans += dfs(dfs, 0, 1, ((__uint128_t)1 << n) - 1, 0, v, lower_vertex);
      if(ans >= mod) ans -= mod;

      // erase lower_vertex
      m -= pop_count_128(g[lower_vertex]);
      is_alive ^= (__uint128_t)1 << lower_vertex;
      for(int i = 0; i < n; i++){
        if((g[i] >> lower_vertex) & 1){
          g[i] ^= (__uint128_t)1 << lower_vertex;
        }
      }
    }
  }
  std::cout << ans << '\n';
}
*/

// 自己辺はあっても無視, 多重辺があると時間計算量が壊れる(enumerateの場合)
template<typename edge>
void undirected_enumerate_triangles(vec<vec<edge>> &g){
  int n = g.size();
  vec<vec<edge>> g2(n);
  for(int i = 0; i < n; i++){
    int deg_i = g[i].size();
    for(edge &j : g[i]){
      int deg_j = g[j.to()].size();
      if(deg_i > deg_j || (deg_i == deg_j && i > j.to())) g2[i].push_back(j.to());
    }
  }
  std::vector<bool> use(n, false);
  for(int i = 0; i < n; i++){
    for(edge &j : g2[i]) use[j.to()] = true;
    for(edge &j : g2[i]){
      for(edge &k : g2[j.to()]){
        if(use[k.to()]){
          // {i, j, k} is triangle
          //
        }
      }
    }
    for(edge &j : g2[i]) use[j.to()] = false;
  }
}

struct cutset_structure{
private:
  static constexpr int C = 12;
  static constexpr long long mask = (1LL << 32) - 1;
  static constexpr long long encode(int a, int b){
    if(a > b) std::swap(a, b);
    return ((long long)a << 32) + b;
  }
  static constexpr std::pair<int, int> decode(long long x){
    return {x >> 32, x & mask};
  }
  // probability of success is at least 1/9.
  // if make C lg_2(n) versions, The probability that all versions fail to find the edge is at most
  // (1 - 1/9)^C lg_2(n) = n^-0.17C
  int Clgn;
  struct cutset_structure_independent{
    int N, maxlv;
    std::vector<std::vector<long long>> XORed;
    bool is_alive_random(int lv, long long x){
      return random_number() & 1;
    }
    void make_new_level(){
      XORed.push_back(std::vector<long long>(N, 0));
      maxlv++;
    }
    cutset_structure_independent(){}
    cutset_structure_independent(int N): N(N), maxlv(1), XORed(1, std::vector<long long>(N, 0)){}
    // 操作の時点で(a, b)がないことが保証されている
    void insert(int a, int b){
      long long x = encode(a, b);
      int lv = 0;
      while(true){
        XORed[lv][a] ^= x;
        XORed[lv][b] ^= x;
        if(!is_alive_random(lv, x)) return;
        if(++lv == maxlv) make_new_level();
      }
    }
    // 操作の時点で(a, b)があることが保証されている
    void erase(int a, int b){
      insert(a, b);
    }
  };
  bool find(int a, const std::vector<int> &T){
    int idx = std::lower_bound(T.begin(), T.end(), a) - T.begin();
    return idx != T.size() && T[idx] == a;
  }
  std::set<long long> st;
  std::vector<cutset_structure_independent> G;
public:
  cutset_structure(){}
  cutset_structure(int N): Clgn(C * log2(N)){
    for(int i = 0; i < Clgn; i++){
      G.push_back(cutset_structure_independent(N));
    }
  }
  void insert(int a, int b){
    for(int i = 0; i < Clgn; i++){
      G[i].insert(a, b);
    }
    st.insert(encode(a, b));
  }
  void erase(int a, int b){
    insert(a, b);
    st.erase(encode(a, b));
  }
  std::pair<int, int> find_edge(std::vector<int> T){
    std::sort(T.begin(), T.end());
    for(int i = 0; i < Clgn; i++){
      for(int j = 0; j < G[i].maxlv; j++){
        long long x = 0;
        for(int v: T) x ^= G[i].XORed[j][v];
        if(!x) continue;
        if(st.find(x) != st.end()){
          auto [a, b] = decode(x);
          if(find(a, T) ^ find(b, T)) return {a, b};
        }
      }
    }
    return {-1, -1};
  }
};
struct cutset_structure_bitset{
private:
  int N;
  std::vector<dynamic_bitset> bit;
public:
  cutset_structure_bitset(){}
  cutset_structure_bitset(int N): N(N), bit(N){
    for(int i = 0; i < N; i++) bit[i] = dynamic_bitset(N, 0);
  }
  void insert(int a, int b){
    bit[a].set(b, 1);
    bit[b].set(a, 1);
  }
  void erase(int a, int b){
    bit[a].set(b, 0);
    bit[b].set(a, 0);
  }
  std::pair<int, int> find_edge(std::vector<int> T){
    dynamic_bitset mask(N, 1);
    for(int v : T) mask.set(v, 0);
    for(int v : T){
      int f = (mask & bit[v]).find_first();
      if(f != -1) return {v, f};
    }
    return {-1, -1};
  }
};
// 有向完全グラフから辺集合Eを取り除いた補グラフでsから各頂点への最短距離, 到達不可能な頂点はinf
std::vector<int> dense_graph_bfs_directed(int n, const std::vector<std::pair<int, int>> &E, int s){
  static constexpr int inf = std::numeric_limits<int>::max();
  std::set<int> unused;
  for(int i = 0; i < n; i++){
    if(i != s){
      unused.insert(i);
    }
  }
  std::vector<int> ans(n, inf);
  ans[s] = 0;
  std::vector<std::vector<int>> ban(n);
  for(auto [a, b] : E) ban[a].push_back(b);
  for(int i = 0; i < n; i++){
    std::sort(ban[i].begin(), ban[i].end());
    ban[i].erase(std::unique(ban[i].begin(), ban[i].end()), ban[i].end());
  }
  std::vector<std::vector<int>> D(n); // 各距離の頂点
  D[0].push_back(s);
  for(int i = 0; ; i++){
    if(D[i].empty() || unused.empty()) break;
    for(int v : D[i]){
      int j = 0;
      for(auto itr = unused.begin(); itr != unused.end();){
        int u = *itr;
        while(j < ban[v].size() && ban[v][j] < u) j++;
        if(j < ban[v].size() && ban[v][j] == u){
          itr++;
        }else{
          D[i + 1].push_back(u);
          ans[u] = i + 1;
          itr = unused.erase(itr);
        }
      }
    }
  }
  return ans;
}
std::vector<int> dense_graph_bfs_undirected(int n, const std::vector<std::pair<int, int>> &E, int s){
  std::vector<std::pair<int, int>> tmp;
  for(auto [a, b] : E){
    tmp.push_back({a, b});
    tmp.push_back({b, a});
  }
  return dense_graph_bfs_directed(n, tmp, s);
}
#endif