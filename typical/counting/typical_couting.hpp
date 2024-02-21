#include "../../math/integer_extra.hpp"
namespace count_equation{
  using ll = long long;
  // A([0, ra)) + B([0, rb)) = C([0, rc)) の解の数
  ll eq1(ll ra, ll rb, ll rc){
    if(ra <= 0 || rb <= 0 || rc <= 0) return 0;
    ra = std::min(ra, rc);
    rb = std::min(rb, rc);
    // 0 <= A <= (rc - rb)のとき, 対応するBはrb種類
    ll ans = rb * std::min(ra, rc - rb + 1);
    // (rc - rb) < A < raのとき, 対応するBが1ずつ減っていく
    ll len = std::max(0LL, ra - (rc - rb + 1));
    ans += len * rb - (len * (len + 1)) / 2;
    return ans;
  }
  // A([la, ra)) + B([lb, rb)) = C([lc, rc)) の解の数
  ll eq2(ll la, ll ra, ll lb, ll rb, ll lc, ll rc){
    ll low = la + lb;
    ra -= la, rb -= lb;
    lc -= low, rc -= low;
    return eq1(ra, rb, rc) - eq1(ra, rb, lc);
  }
  // A([la, ra)) * B([lb, rb)) = C([0, rc)) の解の数
  // @param 変数 >= 0
  ll __eq3(ll la, ll ra, ll lb, ll rb, ll rc){
    if(la >= ra || ra <= 0 || lb >= rb || rb <= 0) return 0;
    assert(la >= 0 && lb >= 0 && rc >= 0);
    if(rc == 0) return 0;
    ll ans = 0;
    for(auto [bmax, amin, amax] : enumerate_quotients3(rc - 1)){
      // [amin, amax)のとき, [lb, bmax)のbが可能
      bmax = std::max(std::min(bmax + 1, rb), lb);
      amin = std::max(la, amin);
      amax = std::min(amax, ra);
      if(amin < amax) ans += (amax - amin) * (bmax - lb);
    }
    return ans;
  }
  // A([la, ra)) * B([lb, rb)) = C([lc, rc)) の解の数
  ll eq3(ll la, ll ra, ll lb, ll rb, ll lc, ll rc){
    if(la >= ra || lb >= rb || lc >= rc) return 0;
    ll ans = 0;
    // 0以上 × 0以上 = 0以上
    if(ra > 0 && rb > 0 && rc > 0){
      ans += __eq3(std::max(0LL, la), ra, std::max(0LL, lb), rb, rc) - __eq3(std::max(0LL, la), ra, std::max(0LL, lb), rb, std::max(0LL, lc));
    }
    // 負 × 負 = 正
    if(la < 0 && lb < 0 && rc > 1){
      int La = std::max(1LL, -ra + 1), Ra = -la + 1;
      int Lb = std::max(1LL, -rb + 1), Rb = -lb + 1;
      ans += __eq3(La, Ra, Lb, Rb, rc) - __eq3(La, Ra, Lb, Rb, std::max(0LL, lc));
    }
    // 0以上 × 負 = 0以下
    if(ra > 0 && lb < 0 && lc <= 0){
      int Lb = std::max(1LL, -rb + 1), Rb = -lb + 1;
      int Lc = std::max(0LL, -rc + 1), Rc = -lc + 1;
      ans += __eq3(std::max(0LL, la), ra, Lb, Rb, Rc) - __eq3(std::max(0LL, la), ra, Lb, Rb, Lc);
    }
    // 負 × 0以上 = 0以下
    if(la < 0 && rb > 0 && lc <= 0){
      int La = std::max(1LL, -ra + 1), Ra = -la + 1;
      int Lc = std::max(0LL, -rc + 1), Rc = -lc + 1;
      ans += __eq3(La, Ra, std::max(0LL, lb), rb, Rc) - __eq3(La, Ra, std::max(0LL, lb), rb, Lc);
    }
    return ans;
  }
  // gcd(A([la, ra)), B([lb, rb))) = C([lc, rc)) の解の数
  ll eq4(int la, int ra, int lb, int rb, int lc, int rc){
    static std::vector<std::vector<int>> low;
    if(low.empty()) low = coprime_sum_table(2000, 2000);
    std::unordered_map<ll, ll> mp;
    auto encode = [&](int x, int y) ->ll {
      if(x < y) std::swap(x, y);
      return (ll(x) << 32) + y;
    };
    auto csum = [&](auto &&csum, int x, int y)->ll {
      if(x < y) std::swap(x, y);
      if(x < low.size()) return low[x][y];
      if(y == 0) return 0;
      ll enc = encode(x, y);
      auto itr = mp.find(enc);
      if(itr != mp.end()) return itr->second;
      ll xy = (ll)x * y, res = xy;
      auto d = enumerate_quotients_pair(x, y);
      for(auto [a, b] : d){
        if(a == 1) continue;
        res -= (b - a) * csum(csum, x / a, y / a);
      }
      mp.emplace(enc, res);
      return res;
    };
    ll ans = 0;
    if(ra > 1 && rb > 1){
      for(auto [l, r] : enumerate_quotients_pair(ra - 1, rb - 1)){
        l = std::max(l, lc), r = std::min(r, rc);
        if(l < r) ans += (r - l) * csum(csum, (ra - 1) / l, (rb - 1) / l);
      }
    }
    if(la > 1 && rb > 1){
      for(auto [l, r] : enumerate_quotients_pair(la - 1, rb - 1)){
        l = std::max(l, lc), r = std::min(r, rc);
        if(l < r) ans -= (r - l) * csum(csum, (la - 1) / l, (rb - 1) / l);
      }
    }
    if(ra > 1 && lb > 1){
      for(auto [l, r] : enumerate_quotients_pair(ra - 1, lb - 1)){
        l = std::max(l, lc), r = std::min(r, rc);
        if(l < r) ans -= (r - l) * csum(csum, (ra - 1) / l, (lb - 1) / l);
      }
    }
    if(la > 1 && lb > 1){
      for(auto [l, r] : enumerate_quotients_pair(la - 1, lb - 1)){
        l = std::max(l, lc), r = std::min(r, rc);
        if(l < r) ans += (r - l) * csum(csum, (la - 1) / l, (lb - 1) / l);
      }
    }
    return ans;
  }
}


#include "../../data_structure/bit_sequence/dynamic_bitset_rank_select.hpp"
#include "../../data_structure/range_query/accumulate1d.hpp"

// 長さm+1のretを返す
// ret[i] := vの各要素をiで割った時の余りの最小値, ret[0]は未定義, O(Vmax * log^2(Vmax))
std::vector<int> sequence_min_of_mod(const std::vector<int> &v, int m){
  assert(!v.empty());
  int max_elem = v[0];
  for(int i : v) max_elem = std::max(max_elem, i);
  dynamic_bitset_rank_select bit(max_elem + 1);
  for(int i : v) bit.set(i, 1);
  std::vector<int> ret(m + 1, max_elem);
  for(int i = 1; i <= m; i++){
    for(int j = 0; j <= max_elem; j += i){
      int next = bit.find_next1(j);
      if(next == -1) break;
      if(next - j < i) ret[i] = std::min(ret[i], next - j);
    }
  }
  return ret;
}
// 長さm+1のretを返す
// ret[i] := vの各要素をiで割った時の余りの最大値, ret[0]は未定義, O(Vmax * log^2(Vmax))
std::vector<int> sequence_max_of_mod(const std::vector<int> &v, int m){
  assert(!v.empty());
  int max_elem = v[0];
  for(int i : v) max_elem = std::max(max_elem, i);
  dynamic_bitset_rank_select bit(max_elem + 1);
  for(int i : v) bit.set(i, 1);
  std::vector<int> ret(m + 1, -1);
  for(int i = 1; i <= m; i++){
    for(int j = 0; j <= max_elem; j += i){
      int prev = bit.find_prev1(std::min(max_elem, j + i - 1));
      if(prev == -1) break;
      if(prev >= j) ret[i] = std::max(ret[i], prev - j);
    }
  }
  return ret;
}
// 長さm+1のretを返す
// ret[i] := vの各要素をiで割った時の余りの和, ret[0]は0, O(Vmax * log(Vmax))
std::vector<long long> sequence_sum_of_mod(const std::vector<int> &v, int m){
  assert(!v.empty());
  int max_elem = v[0];
  for(int i : v) max_elem = std::max(max_elem, i);
  std::vector<int> cnt(max_elem + 1, 0);
  long long sum = 0;
  for(int i : v){
    cnt[i]++;
    sum += i;
  }
  accumulate1d<int> ac(cnt);
  std::vector<long long> ret(m + 1, 0);
  for(int i = 1; i <= m; i++){
    ret[i] = sum;
    for(int j = 0; j <= max_elem; j += i){
      int rank = ac.query(j, std::min(max_elem + 1, j + i)); // [j, j + i)の要素数
      ret[i] -= (long long)j * rank;
    }
  }
  return ret;
}
// 長さm+1のretを返す
// ret[i] := vの中でiの約数の数, ret[0]は0, O(Vmax * log(Vmax))
std::vector<int> sequence_divisor_count(const std::vector<int> &v, int m){
  assert(!v.empty());
  int max_elem = v[0];
  for(int i : v) max_elem = std::max(max_elem, i);
  std::vector<int> cnt(max_elem + 1, 0);
  for(int i : v) cnt[i]++;
  std::vector<int> ret(m + 1, 0);
  for(int i = 1; i <= std::min(m, max_elem); i++){
    for(int j = i; j <= m; j += i){
      ret[j] += cnt[i];
    }
  }
  return ret;
}
// 長さm+1のretを返す
// ret[i] := vの中でiの倍数の数, ret[0]は0, O(Vmax * log(Vmax))
std::vector<int> sequence_multiple_count(const std::vector<int> &v, int m){
  assert(!v.empty());
  int max_elem = v[0];
  for(int i : v) max_elem = std::max(max_elem, i);
  std::vector<int> cnt(max_elem + 1, 0);
  for(int i : v) cnt[i]++;
  std::vector<int> ret(m + 1, 0);
  for(int i = 1; i <= std::min(m, max_elem); i++){
    for(int j = i; j <= max_elem; j += i){
      ret[i] += cnt[j];
    }
  }
  return ret;
}

// lを固定した時に区間[l, r]のlcmの値が変わるようなrを列挙
// lcmの上限をinfとすると O(Nlog(inf))
template<typename T>
std::vector<std::vector<std::pair<T, int>>> enumerate_lcm(const std::vector<T> &v, T inf){
  int n = v.size();
  std::vector<std::vector<std::pair<T, int>>> ret(n);
  for(int i = n - 1; i >= 0; i--){
    ret[i].push_back({std::min(v[i], inf), i});
    if(i < n - 1){
      for(auto [xi, id] : ret[i + 1]){
        xi = std::min(inf, lcm_limited(ret[i].back().first, xi, inf));
        if(xi > ret[i].back().first) ret[i].push_back({xi, id});
      }
    }
  }
  return ret;
}

// lを固定した時に区間[l, r]のgcdの値が変わるようなrを列挙 O(Nlog(max(v)))
template<typename T>
std::vector<std::vector<std::pair<T, int>>> enumerate_gcd(const std::vector<T> &v){
  int n = v.size();
  std::vector<std::vector<std::pair<T, int>>> ret(n);
  for(int i = n - 1; i >= 0; i--){
    ret[i].push_back({v[i], i});
    if(i < n - 1){
      for(auto [xi, id] : ret[i + 1]){
        xi = gcd(ret[i].back().first, xi);
        if(xi < ret[i].back().first) ret[i].push_back({xi, id});
      }
    }
  }
  return ret;
}
