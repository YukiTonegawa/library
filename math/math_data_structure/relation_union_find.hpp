#ifndef _RELATION_UNION_FIND_H_
#define _RELATION_UNION_FIND_H_
#include <vector>
#include <numeric>

// 与えられる関係式 x_i = x_j * a + bにおいて常にa = ±1が成り立つ場合, 確定する答えは整数/2になる
// @param 2 * {(add_range, add_constで追加した値) + (重み * n)} がオーバーフローしない
template<typename Val>
struct relation_union_find_integer{
private:
  int num_cc;
  std::vector<int> par, sz;
  std::vector<std::pair<Val, Val>> relation_with_parent; // その時点での親との関係式
  std::vector<std::pair<Val, Val>> ans; // 答えとしてあり得る値[l, r]
  static constexpr Val inf = std::numeric_limits<Val>::max();
  static constexpr Val minf = std::numeric_limits<Val>::min();

  // x_1 = x_2 * a + b (a != 0)のとき
  // x_2 = x_1 * A + Bを満たす{A, B}を返す
  std::pair<Val, Val> inv(const std::pair<Val, Val> &a){
    if(a.first == 1) return {1, -a.second};
    return {-1, a.second};
  }
  std::pair<Val, Val> composite(const std::pair<Val, Val> &a, const std::pair<Val, Val> &b){
    return {a.first * b.first, a.second * b.first + b.second};
  }
  // inf, minfをaとして使うと壊れる
  Val composite(Val a, const std::pair<Val, Val> &b){
    if(a == inf) return b.first == 1 ? inf : minf;
    if(a == minf) return b.first == 1 ? minf : inf;
    return a * b.first + 2 * b.second;
  }
  int find(int u){
    if(par[u] == u) return u;
    int p = find(par[u]);
    if(par[u] != p){
      relation_with_parent[u] = composite(relation_with_parent[par[u]], relation_with_parent[u]);
      par[u] = p;
    }
    return par[u];
  }
  int size(int u){
    return sz[find(u)];
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
public:
  relation_union_find_integer(int n): num_cc(n), par(n), sz(n, 1), relation_with_parent(n, {1, 0}), ans(n, {minf, inf}){
    std::iota(par.begin(), par.end(), 0);
  }
  // @param |a| <= 1
  bool check_relation(int u, int v, Val a, Val b){
    assert(abs(a) <= 1);
    if(a == 0) return check_const(u, b);
    int ur = find(u), vr = find(v);
    if(ur == vr){
      auto [A, B] = composite(relation_with_parent[v], {a, b});
      auto [X, Y] = relation_with_parent[u];
      if(A - X == 0) return (Y - B) == 0;
      Val x = 2 * (Y - B) / (A - X);
      return ans[ur].first <= x && x <= ans[ur].second;
    }
    std::pair<Val, Val> c = composite(relation_with_parent[v], composite({a, b}, inv(relation_with_parent[u])));
    Val lur = composite(ans[vr].first, c), rur = composite(ans[vr].second, c);
    if(c.first == -1) std::swap(lur, rur);
    return std::max(ans[ur].first, lur) <= std::min(ans[ur].second, rur);
  }
  // x_u = x_v * a + b
  bool add_relation(int u, int v, Val a, Val b){
    assert(abs(a) <= 1);
    if(a == 0) return add_const(u, b);
    int ur = find(u), vr = find(v);
    if(ur == vr){
      auto [A, B] = composite(relation_with_parent[v], {a, b});
      auto [X, Y] = relation_with_parent[u];
      if(A - X == 0) return (Y - B) == 0;
      Val x = 2 * (Y - B) / (A - X);
      if(x < ans[ur].first || ans[ur].second < x) return false;
      ans[ur] = {x, x};
      return true;
    }
    if(sz[ur] > sz[vr]){
      std::swap(u, v);
      std::swap(ur, vr);
      std::tie(a, b) = inv({a, b});
    }
    std::pair<Val, Val> c = composite(relation_with_parent[v], composite({a, b}, inv(relation_with_parent[u])));
    Val lur = composite(ans[vr].first, c), rur = composite(ans[vr].second, c);
    if(c.first == -1) std::swap(lur, rur);
    lur = std::max(ans[ur].first, lur);
    rur = std::min(ans[ur].second, rur);
    if(lur > rur) return false;
    num_cc--;
    ans[ur] = {lur, rur};
    sz[vr] += sz[ur];
    par[ur] = par[vr];
    relation_with_parent[ur] = c;
    std::pair<Val, Val> i = inv(c);
    Val lvr = composite(ans[ur].first, i), rvr = composite(ans[ur].second, i);
    if(i.first == -1) std::swap(lvr, rvr);
    ans[vr].first = std::max(ans[vr].first, lvr);
    ans[vr].second = std::min(ans[vr].second, rvr);
    return true;
  }
  // l <= a < rという情報が矛盾するか確認, 1つでも条件を満たす整数/2が存在するならtrue
  bool check_range(int u, Val l, Val r){
    l *= 2, r = 2 * (r - 1);
    int ur = find(u);
    std::pair<Val, Val> f = inv(relation_with_parent[u]);
    Val lur = composite(l, f), rur = composite(r, f);
    if(f.first == -1) std::swap(lur, rur);
    return std::max(ans[ur].first, lur) <= std::min(ans[ur].second, rur);
  }
  // l <= a < rという情報が矛盾するか確認, 1つでも条件を満たす整数が存在するならok
  // 矛盾しない場合は採用してtrue
  // 矛盾する場合は何もせずfalse
  bool add_range(int u, Val l, Val r){
    l *= 2, r = 2 * (r - 1);
    int ur = find(u);
    std::pair<Val, Val> f = inv(relation_with_parent[u]);
    Val lur = composite(l, f), rur = composite(r, f);
    if(f.first == -1) std::swap(lur, rur);
    lur = std::max(ans[ur].first, lur), rur = std::min(ans[ur].second, rur);
    if(lur > rur) return false;
    ans[ur] = {lur, rur};
    return true;
  }
  // u = aという情報が矛盾するか確認
  bool check_const(int u, Val a){
    a *= 2;
    int ur = find(u);
    Val aur = composite(a, inv(relation_with_parent[u]));
    return ans[ur].first <= aur && aur <= ans[ur].second;
  }
  // u = a
  // 矛盾しない場合は採用してtrue
  // 矛盾する場合は何もせずfalse
  bool add_const(int u, Val a){
    a *= 2;
    int ur = find(u);
    Val aur = composite(a, inv(relation_with_parent[u]));
    if(aur < ans[ur].first || ans[ur].second < aur) return false;
    ans[ur] = {aur, aur};
    return true;
  }
  // あり得るuの値 * 2の区間[l, r)を返す(l <= 2 * u < r)
  // 区間長が1かつそれが偶数でない場合, 整数解を持たない(整数/2の解はある)
  std::pair<Val, Val> get2x(int u){
    int ur = find(u);
    Val lu = composite(ans[ur].first, relation_with_parent[u]);
    Val ru = composite(ans[ur].second, relation_with_parent[u]);
    if(relation_with_parent[u].first == -1) std::swap(lu, ru);
    return {lu, (ru == inf ? inf : ru + 1)};
  }
  // u = v * a + bという式が存在するか, 存在する場合その値
  std::tuple<bool, Val, Val> relation(int u, int v){
    int ur = find(u), vr = find(v);
    if(ur != vr) return {false, 0, 0};
    auto [x, y] = composite(inv(relation_with_parent[v]), relation_with_parent[u]);
    return {true, x, y};
  }
  // 連結成分の数
  int count_cc(){
    return num_cc;
  }
};

struct relation_union_find_mod2{
private:
  int num_cc;
  std::vector<int> par, sz;
  std::vector<int> relation_with_parent; // その時点での親との関係式
  std::vector<bool> decided;
  std::vector<bool> ans;
  int find(int u){
    if(par[u] == u) return u;
    int p = find(par[u]);
    if(par[u] != p){
      relation_with_parent[u] ^= relation_with_parent[par[u]];
      par[u] = p;
    }
    return par[u];
  }
  int size(int u){
    return sz[find(u)];
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
public:
  relation_union_find_mod2(int n): num_cc(n), par(n), sz(n, 1), relation_with_parent(n, 0), decided(n, 0), ans(n, 0){
    std::iota(par.begin(), par.end(), 0);
  }
  // u = v ^ a
  bool check_relation(int u, int v, bool a){
    int ur = find(u), vr = find(v);
    if(ur == vr) return (relation_with_parent[v] ^ a) == relation_with_parent[u];
    // ur = vr ^ c
    int c = relation_with_parent[v] ^ relation_with_parent[u] ^ a;
    if(decided[ur] && decided[vr] && (ans[vr] ^ c) != ans[ur]) return false;
    return true;
  }
  // x_u = x_v ^ a
  bool add_relation(int u, int v, bool a){
    int ur = find(u), vr = find(v);
    if(ur == vr) return (relation_with_parent[v] ^ a) == relation_with_parent[u];
    if(sz[ur] > sz[vr]){
      std::swap(u, v);
      std::swap(ur, vr);
    }
    // ur = vr ^ c
    int c = relation_with_parent[v] ^ relation_with_parent[u] ^ a;
    if(decided[ur] && decided[vr] && (ans[vr] ^ c) != ans[ur]) return false;
    num_cc--;
    sz[vr] += sz[ur];
    par[ur] = par[vr];
    relation_with_parent[ur] = c;
    if(!decided[vr] && decided[ur]){
      decided[vr] = true;
      ans[vr] = ans[ur] ^ c;
    }
    return true;
  }
  // u = aという情報が矛盾するか確認
  bool check_const(int u, bool a){
    int ur = find(u);
    if(!decided[ur]) return true;
    return (ans[ur] ^ relation_with_parent[u]) == a;
  }
  // u = a
  // 矛盾しない場合は採用してtrue
  // 矛盾する場合は何もせずfalse
  bool add_const(int u, bool a){
    int ur = find(u);
    if(!decided[ur]){
      ans[ur] = a;
      decided[ur] = true;
    }
    return (ans[ur] ^ relation_with_parent[u]) == a;
  }
  // {確定しているか, そうでない場合, その値}
  std::pair<bool, bool> get(int u){
    int ur = find(u);
    if(!decided[ur]) return {false, 0};
    return {true, (ans[ur] ^ relation_with_parent[u])};
  }
  // u = v ^ aという式が存在するか, 存在する場合そのa
  std::pair<bool, bool> relation(int u, int v){
    int ur = find(u), vr = find(v);
    if(ur != vr) return {false, 0};
    return {true, relation_with_parent[u] ^ relation_with_parent[v]};
  }
  // 連結成分の数
  int count_cc(){
    return num_cc;
  }
};

template<typename mint>
struct relation_union_find_mod{
private:
  int num_cc;
  std::vector<int> par, sz;
  std::vector<std::pair<mint, mint>> relation_with_parent; // その時点での親との関係式
  std::vector<bool> decided;
  std::vector<mint> ans;
  // x_1 = x_2 * a + b (a != 0)のとき
  // x_2 = x_1 * A + Bを満たす{A, B}を返す
  std::pair<mint, mint> inv(const std::pair<mint, mint> &a){
    mint A = a.first.inv();
    return {A, -a.second * A};
  }
  std::pair<mint, mint> composite(const std::pair<mint, mint> &a, const std::pair<mint, mint> &b){
    return {a.first * b.first, a.second * b.first + b.second};
  }
  int find(int u){
    if(par[u] == u) return u;
    int p = find(par[u]);
    if(par[u] != p){
      relation_with_parent[u] = composite(relation_with_parent[par[u]], relation_with_parent[u]);
      par[u] = p;
    }
    return par[u];
  }
  int size(int u){
    return sz[find(u)];
  }
  bool same(int u, int v){
    return find(u) == find(v);
  }
public:
  relation_union_find_mod(int n): num_cc(n), par(n), sz(n, 1), relation_with_parent(n, {1, 0}), decided(n, false), ans(n, 0){
    std::iota(par.begin(), par.end(), 0);
  }
  bool check_relation(int u, int v, mint a, mint b){
    if(a == 0) return check_const(u, b);
    int ur = find(u), vr = find(v);
    if(ur == vr){
      auto [A, B] = composite(relation_with_parent[v], {a, b});
      auto [X, Y] = relation_with_parent[u];
      // R * (A - X) = Y - B
      if(decided[ur]) return ans[ur] * (A - X) == Y - B;
      if(A - X == 0) return (Y - B) == 0;
      return true; // ans[ur] = (Y - B) / (A - X)
    }
    if(!decided[ur] || !decided[vr]) return true;
    auto [c, d] = composite(relation_with_parent[v], composite({a, b}, inv(relation_with_parent[u])));
    // ur = vr * c + d
    return ans[vr] * c + d == ans[ur];
  }
  // x_u = x_v * a + b
  bool add_relation(int u, int v, mint a, mint b){
    if(a == 0) return add_const(u, b);
    int ur = find(u), vr = find(v);
    if(ur == vr){
      auto [A, B] = composite(relation_with_parent[v], {a, b});
      auto [X, Y] = relation_with_parent[u];
      // R * (a - X) = Y - b
      if(decided[ur]) return ans[ur] * (A - X) == Y - B;
      if(A - X == 0) return (Y - B) == 0;
      decided[ur] = true;
      ans[ur] = (Y - B) / (A - X);
      return true;
    }
    if(sz[ur] > sz[vr]){
      std::swap(u, v);
      std::swap(ur, vr);
      std::tie(a, b) = inv({a, b});
    }
    // ur = vr * c + d
    auto [c, d] = composite(relation_with_parent[v], composite({a, b}, inv(relation_with_parent[u])));
    if(decided[ur] && decided[vr] && ans[vr] * c + d != ans[ur]) return false;
    num_cc--;
    sz[vr] += sz[ur];
    par[ur] = par[vr];
    relation_with_parent[ur] = {c, d};
    if(!decided[vr] && decided[ur]){
      decided[vr] = true;
      auto [e, f] = inv({c, d});
      ans[vr] = ans[ur] * e + f;
    }
    return true;
  }
  // u = aという情報が矛盾するか確認
  bool check_const(int u, mint a){
    int ur = find(u);
    if(!decided[ur]) return true;
    return (ans[ur] * relation_with_parent[u].first + relation_with_parent[u].second) == a;
  }
  // u = a
  // 矛盾しない場合は採用してtrue
  // 矛盾する場合は何もせずfalse
  bool add_const(int u, mint a){
    int ur = find(u);
    if(!decided[ur]){
      auto [c, d] = inv(relation_with_parent[u]);
      ans[ur] = a * c + d;
      decided[ur] = true;
    }
    return (ans[ur] * relation_with_parent[u].first + relation_with_parent[u].second) == a;
  }
  // {確定しているか, そうでない場合, その値}
  std::pair<bool, mint> get(int u){
    int ur = find(u);
    if(!decided[ur]) return {false, 0};
    return {true, (ans[ur] * relation_with_parent[u].first + relation_with_parent[u].second)};
  }
  // u = v * a + bという式が存在するか, 存在する場合その値
  std::tuple<bool, mint, mint> relation(int u, int v){
    int ur = find(u), vr = find(v);
    if(ur != vr) return {false, 0, 0};
    auto [x, y] = composite(inv(relation_with_parent[v]), relation_with_parent[u]);
    return {true, x, y};
  }
  // 連結成分の数
  int count_cc(){
    return num_cc;
  }
};
#endif