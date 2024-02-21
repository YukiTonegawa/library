#ifndef _FPS_EXTRA_H_
#define _FPS_EXTRA_H_
#include "fps.hpp"
template<typename mint>
std::vector<mint> multipoint_evaluation(const formal_power_series<mint> &f, const std::vector<mint> &v){
  using fps = formal_power_series<mint>;
  int m = v.size();
  int N = 1;
  while(N < m) N *= 2;
  std::vector<fps> t(2 * N - 1, fps{1});
  for(int i = 0; i < m; i++) t[N - 1 + i] = fps{-v[i], 1};
  for(int i = N - 2; i >= 0; i--) t[i] = t[i * 2 + 1] * t[i * 2 + 2];
  t[0] = f % t[0];
  for(int i = 1; i < 2 * N - 1; i++){
    t[i] = t[(i - 1) / 2] % t[i];
  }
  std::vector<mint> res(m);
  for(int i = 0; i < m; i++){
    res[i] = t[N - 1 + i][0];
  }
  return res;
}
template<typename mint, typename fps>
fps interpolation_divide_and_conquer(const std::vector<mint> &fx, const fps &f, const std::vector<fps> &tree, int k, int l, int r){
  if(r - l == 1) return fps{fx[l] * f[0].inv()};
  int mid = (l + r) / 2;
  if(tree[k * 2 + 2].empty()) return interpolation_divide_and_conquer<mint, fps>(fx, f, tree, k * 2 + 1, l, mid);
  fps left = interpolation_divide_and_conquer<mint, fps>(fx, f % tree[k * 2 + 1], tree, k * 2 + 1, l, mid);
  fps right = interpolation_divide_and_conquer<mint, fps>(fx, f % tree[k * 2 + 2], tree, k * 2 + 2, mid, r);
  return left * tree[k * 2 + 2] + right * tree[k * 2 + 1];
}
template<typename mint>
formal_power_series<mint> polynomial_interpolation(const std::vector<mint> &xi, const std::vector<mint> &fx){
  using fps = formal_power_series<mint>;
  int n = xi.size();
  int N = 1;
  while(N < n) N *= 2;
  std::vector<fps> tree(2 * N - 1, fps{});
  for(int i = 0; i < n; i++) tree[N - 1 + i] = fps{-xi[i], 1};
  for(int i = N - 2; i >= 0; i--){
    if(tree[i * 2 + 2].empty()) tree[i] = tree[i * 2 + 1];
    else tree[i] = tree[i * 2 + 1] * tree[i * 2 + 2];
  }
  for(int i = 0; i < n; i++) tree[0][i] = tree[0][i + 1] * (i + 1);
  tree[0].pop_back();
  return interpolation_divide_and_conquer<mint, fps>(fx, tree[0], tree, 0, 0, N).prefix(n);
}
// a / bx + c を大量に足す
template<typename mint>
formal_power_series<mint> fractions_sum_binomial(const std::vector<std::tuple<mint, mint, mint>> &v, int max_size = -1){
  using fps = formal_power_series<mint>;
  if(max_size == -1) max_size = std::numeric_limits<int>::max();
  int n = v.size();
  std::vector<std::pair<fps, fps>> res(n);
  for(int i = 0; i < n; i++){
    auto [a, b, c] = v[i];
    assert(b.val() || c.val());
    res[i].first = fps{a};
    res[i].second = fps{c, b};
  }
  for(int d = 1; d < n; d <<= 1){
    for(int j = 0; j < n; j += 2 * d){
      if(j + d < n){
        v[j].first *= v[j + d].second;
        v[j + d].first *= v[j].second;
        v[j].second *= v[j + d].second;
        v[j].first += v[j + d].first;
        if(v[j].first.size() > max_size) v[j].first.resize(max_size);
        if(v[j].second.size() > max_size) v[j].second.resize(max_size);
      }
    }
  }
  return res[0].first * res[0].second.inv();
}
/*
template<typename mint>
std::vector<mint> shiftof_sampling_points(mint x0, mint d, const std::vector<mint> &y, const std::vector<mint> &p){
  using fps = formal_power_series<mint>;
  static modcomb<mint> mcb;
  int n = y.size();
  mcb.recalc(n);
  mint pi_xi_x0(1);
  // xi == 0 or xi == xjで壊れる
  // xi == 0は調べる
  // 素数modだと仮定して, n >= mod or dがmodの倍数で壊れる
  assert(n < mint::mod() && d.val()) ;
  for(int i = 0; i < n; i++){
    mint xi = x0 + d * i;
    assert(xi.val());
    pi_xi_x0 *= d * i;
  }
  std::vector<fps> tmp;
  std::vector<std::tuple<mint, mint, mint>> tmp_frac;
  for(int i = 0; i < n; i++){
    mint xi = x0 + d * i;
    tmp.push_back({-xi, 1});
  }
  fps x = multiply_constant_degree<mint>(tmp);
  x /= pi_xi_x0;
  for(int i = 0; i < n; i++){
    mint iQ = mcb.finv(i) * mcb.finv(n - 1 - i);
    mint xi = x0 + d * i;
    if((n - i - 1) & 1) iQ = -iQ;
    mint iC = y[i] * iQ;
    tmp_frac.push_back({iC, 1, -xi});
  }
  x *= fractions_sum_binomial(tmp_frac);
  x.print();
  exit(0);
}
*/

/*
//f_i(x) := 1 ± a_i x^b_i を大量に掛ける
//O(N^2/th + thMlogM)
// N = M = 1e5 なら50が早かった
template<ll MOD>
StaticModFPS<MOD> multiply_binomial(vector<vector<ll>> v, ll m, ll th){
  using fps = StaticModFPS<MOD>;
  fps ret_big(m, 0);
  vvl big;
  vector<vector<fps>> small(th+1, vector<fps>());
  vector<ll> inv((m/th)+2);
  inv[0] = inv[1] = 1;
  for(ll i=2;i<=(m/th)+1;i++){
    ll u = MOD - (MOD/i);
    inv[i] = (u*inv[MOD%i])%MOD;
  }
  for(auto e:v){
    e[0] %= MOD;
    if(e[0] < 0) e[0] += MOD;
    if(e[1]>th) big.push_back(e);
    else small[e[1]].push_back({1, e[0]});
  }
  for(int i=0;i<big.size();i++){
    ll j = 1;
    ll mul = big[i][0];
    while(j * big[i][1] < m){
      ll y;
      if(j%2==0) y = ((MOD - abs(mul)) * inv[j])%MOD;
      else y = (abs(mul) * inv[j])%MOD;
      ret_big[j * big[i][1]] = (ret_big[j * big[i][1]] + y)%MOD;
      mul = (mul * big[i][0])%MOD;
      j++;
    }
  }
  ret_big = ret_big.exp(m);
  for(ll i=1;i<=th;i++){
    if(small[i].size()==0) continue;
    fps tmp = multiply_lower_degree<MOD>(small[i], m/i+1);
    bool f = true;

    for(int j=m-1;j>=1;j--){
      if(f&&j%i==0&&tmp.size()>(j/i)&&tmp[j/i]){
        tmp.resize(j+1);
        break;
      }
    }
    for(int j=tmp.size()-1;j>=1;j--){
      if(j%i==0) tmp[j] = tmp[j/i];
      else tmp[j] = 0;
    }
    ret_big *= tmp;
    if(ret_big.size()>m) ret_big.resize(m);
  }
  return ret_big;
}
*/
/*
//1 ± x^a_i を大量に掛ける
template<ll MOD>
StaticModFPS<MOD> multiply_binomial_one(const vector<vector<ll>> &v, ll m){
  using fps = StaticModFPS<MOD>;
  vector<ll> pl(m, 0), mi(m, 0);
  for(int i=0;i<v.size();i++){
    if(v[i][0]==-1) mi[v[i][1]]++;
    else pl[v[i][1]]++;
  }
  vector<ll> inv(m+1);
  inv[0] = inv[1] = 1;
  for(ll i=2;i<=m;i++){
    ll u = MOD - (MOD/i);
    inv[i] = (u*inv[MOD%i])%MOD;
  }
  fps ret(m, 0);
  for(ll i=1;i<m;i++){
    ll j = 1;
    while(j*i<m){
      ll mlt = (j%2==0?(-pl[i]-mi[i]):pl[i]-mi[i]);
      mlt = (mlt * inv[j])%MOD;
      mlt = (mlt + MOD)%MOD;
      ret[i*j] += mlt;
      if(ret[i*j]>=MOD) ret[i*j] -= MOD;
      j++;
    }
  }
  return ret.exp(m);
}
*/

template<typename mint>
mint _kth_coefficient(long long k, formal_power_series<mint> &p, formal_power_series<mint> &q){
  using fps = formal_power_series<mint>;
  if(k <= 100000) return (p.prefix(k + 1) * q.inv(k + 1))[k];
  int n = p.deg_fix() + 1, m = q.deg_fix() + 1;
  fps _q(m, 0);
  for(int i = 0; i < m; i++) _q[i] = i & 1 ? -q[i] : q[i];
  // p *= _q
  // q *= _q
  // pとqがほとんど同じサイズである場合は早い
  int N = std::max(n, m), M = _q.size();
  int z = 1 << ceil_pow2(N + M - 1);
  p.resize(z), q.resize(z), _q.resize(z);
  butterfly(p), butterfly(q), butterfly(_q);
  mint iz = mint(z).inv();
  for(int i = 0; i < z; i++){
    p[i] *= _q[i];
    q[i] *= _q[i];
  }
  butterfly_inv(p), butterfly_inv(q);
  p.resize(n + M - 1), q.resize(m + M - 1);
  //
  n = p.size(), m = q.size();
  for(int i = k & 1; i < n; i += 2) p[i / 2] = p[i] * iz;
  for(int i = 0; i < m; i += 2) q[i / 2] = q[i] * iz;
  p.resize((n + 1) / 2);
  q.resize((m + 1) / 2);
  return _kth_coefficient<mint>(k / 2, p, q);
}
template<typename mint>
mint kth_coefficient(long long k, formal_power_series<mint> p, formal_power_series<mint> q){
  return _kth_coefficient<mint>(k, p, q);
}
// _a := 初期値, a[0, d)
// c := (a[i] = ∑{0 <= j < d} a[i - 1 - j] * c[j]を満たす)
// 0-indexed
template<typename mint>
mint kth_term_of_linear_reccurrence(long long k, const std::vector<mint> &_a, const std::vector<mint> &c){
  using fps = formal_power_series<mint>;
  int d = c.size();
  assert(_a.size() == d);
  if(k < d) return _a[k];
  fps Q(d + 1);
  Q[0] = 1;
  for(int i = 0; i < d; i++) Q[i + 1] = -c[i];
  fps P = Q * fps(_a);
  P.resize(d);
  return kth_coefficient<mint>(k, P, Q);
}
#endif