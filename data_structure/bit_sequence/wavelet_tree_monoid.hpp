#ifndef _WAVELET_TREE_MONOID_H_
#define _WAVELET_TREE_MONOID_H_
#include "bitvector.hpp"
#include "../segment_tree/segment_tree.hpp"
#include <vector>
#include <queue>
#include <tuple>
#include <numeric>

template<typename monoid>
struct wavelet_tree_monoid{
  using Val = typename monoid::Val;

private:
  using __bitvector = bitvector_memory;
  using __segtree = segment_tree<monoid>;
  int n, h, inf;
  std::vector<__bitvector> bv;
  std::vector<__segtree> seg;
  void build(std::vector<int> v, std::vector<Val> val){
    bv.resize((1 << h) - 1);
    seg.resize((1 << (h + 1)) - 1);
    std::queue<std::array<int, 4>> q;// (dep, v上でのl, r)
    std::vector<std::pair<int, Val>> tmp(n);
    q.push({h - 1, 0, n, 0});
    while(!q.empty()){
      auto [d, l, r, k] = q.front();
      q.pop();
      std::vector<Val> x(r - l);
      std::copy(val.begin() + l, val.begin() + r, x.begin());
      seg[k] = __segtree(x);
      if(d < 0) continue;
      std::vector<bool> bits(r - l);
      int lidx = 0, ridx = 0;
      for(int i = l; i < r; i++){
        bool b = (v[i] >> d) & 1;
        bits[i - l] = b;
        if(b) tmp[ridx++] = {v[i], val[i]};
        else{
          v[l + lidx] = v[i];
          val[l + lidx++] = val[i];
        }
      }
      for(int i = 0; i < ridx; i++) std::tie(v[l + lidx + i], val[l + lidx + i]) = tmp[i];
      bv[k] = __bitvector(bits);
      int mid = l + lidx;
      q.push({d - 1, l, mid, k * 2 + 1});
      q.push({d - 1, mid, r, k * 2 + 2});
    }
  }
  int __range_freq(int l, int r, int b, const int s, const int t, int S, int T){
    if(l == r || t <= S || T <= s) return 0;
    if(s <= S && T <= t) return r - l;
    int l0 = (!l ? 0 : bv[b].rank0(l)), l1 = l - l0;
    int r0 = bv[b].rank0(r), r1 = r - r0;
    return __range_freq(l0, r0, b * 2 + 1, s, t, S, (S + T) / 2) +
            __range_freq(l1, r1, b * 2 + 2, s, t, (S + T) / 2, T);
  }
  int __rank(int r, int x, int b, int d){
    if(!r) return 0;
    int f = (x >> d) & 1;
    r = f ? bv[b].rank1(r) : bv[b].rank0(r);
    return !d ? r : __rank(r, x, b * 2 + 1 + f, d - 1);
  }
  int __select(int k, int x, int b, int d){
    bool f = (x >> d) & 1;
    int nxt_k = f ? bv[b].select1(k) : bv[b].select0(k);
    if(nxt_k == -1) return -1;
    return !b ? nxt_k : __select(nxt_k, x, (b - 1) >> 1, d + 1);
  }
  int __access(int k, int b, int d){
    int f = bv[b].access(k);
    int nxt_k = f ? bv[b].rank1(k) : bv[b].rank0(k);
    return (f << d) + (d ? __access(nxt_k, b * 2 + 1 + f, d - 1) : 0);
  }
  int __find_next(int k, int x, int b, int d){
    int f = (x >> d) & 1;
    int nxt_k = f ? bv[b].rank1(k) : bv[b].rank0(k);
    return (f << d) + (d ? __access(nxt_k, b * 2 + 1 + f, d - 1) : 0);
  }
  int __range_kth_smallest(int l, int r, int k, int b, int d){
    int l0 = bv[b].rank0(l), l1 = l - l0;
    int r0 = bv[b].rank0(r), r1 = r - r0;
    int z = r0 - l0;
    if(z <= k) return (1 << d) + (d ? __range_kth_smallest(l1, r1, k - z, b * 2 + 2, d - 1) : 0);
    else return (d ? __range_kth_smallest(l0, r0, k, b * 2 + 1, d - 1) : 0);
  }
  std::pair<int, int> __range_kth_smallest_super(int l, int r, const int s, const int t, int S, int T, int k, int b, int d){
    if(l == r || t <= S || T <= s) return {-1, 0};
    if(s <= S && T <= t){
      if(r - l <= k) return {-1, r - l};
      if(d == -1) return {0, r - l};
    }
    int l0 = (!l ? 0 : bv[b].rank0(l)), l1 = l - l0;
    int r0 = bv[b].rank0(r), r1 = r - r0;
    auto L = __range_kth_smallest_super(l0, r0, s, t, S, (S + T) >> 1, k, b * 2 + 1, d - 1);
    if(L.first != -1) return L;
    auto R = __range_kth_smallest_super(l1, r1, s, t, (S + T) >> 1, T, k - L.second, b * 2 + 2, d - 1);
    if(R.first == -1) R.second += L.second;
    else R.first += 1 << d;
    return R;
  }
  std::pair<int, int> __range_kth_largest_super(int l, int r, const int s, const int t, int S, int T, int k, int b, int d){
    if(l == r || t <= S || T <= s) return {-1, 0};
    if(s <= S && T <= t){
      if(r - l <= k) return {-1, r - l};
      if(d == -1) return {0, r - l};
    }
    int l0 = (!l ? 0 : bv[b].rank0(l)), l1 = l - l0;
    int r0 = bv[b].rank0(r), r1 = r - r0;
    auto R = __range_kth_largest_super(l1, r1, s, t, (S + T) >> 1, T, k, b * 2 + 2, d - 1);
    if(R.first != -1){
      R.first += 1 << d;
      return R;
    }
    auto L = __range_kth_largest_super(l0, r0, s, t, S, (S + T) >> 1, k - R.second, b * 2 + 1, d - 1);
    if(L.first == -1) L.second += R.second;
    return L;
  }
  void __count_prefix(int l, int r, int x, std::vector<int> &ret){
    for(int i = h - 1, j = 0; i >= 0 && l < r; i--){
      int l0 = bv[j].rank0(l), r0 = bv[j].rank0(r);
      if((x >> i) & 1){
        ret[h - 1 - i] = r0 - l0;
        l = l - l0, r = r - r0;
        j = j * 2 + 2;
      }else{
        ret[h - 1 - i] = (r - l) - (r0 - l0);
        l = l0, r = r0;
        j = j * 2 + 1;
      }
    }
    ret[h] = r - l;
  }
  void __set(int k, int b, int d, Val z){
    seg[b].set(k, z);
    if(d == -1) return;
    int f = bv[b].access(k);
    int nxt_k = f ? bv[b].rank1(k) : bv[b].rank0(k);
    return __set(nxt_k, b * 2 + 1 + f, d - 1, z);
  }
  void __update(int k, int b, int d, Val z){
    seg[b].update(k, z);
    if(d == -1) return;
    int f = bv[b].access(k);
    int nxt_k = f ? bv[b].rank1(k) : bv[b].rank0(k);
    return __update(nxt_k, b * 2 + 1 + f, d - 1);
  }
  Val __query(int l, int r, int b, const int s, const int t, int S, int T){
    if(l == r || t <= S || T <= s) return monoid::id();
    if(s <= S && T <= t) return seg[b].query(l, r);
    int l0 = (!l ? 0 : bv[b].rank0(l)), l1 = l - l0;
    int r0 = bv[b].rank0(r), r1 = r - r0;
    return monoid::merge(__query(l0, r0, b * 2 + 1, s, t, S, (S + T) / 2),
            __query(l1, r1, b * 2 + 2, s, t, (S + T) / 2, T));
  }
  /*
  std::pair<int, int> __range_kth_largest_super(int l, int r, const int s, int S, int T, int k, int b, int d){
    if(l == r || T <= s) return {-1, 0};
    if(s <= S){
      if(r - l <= k) return {-1, r - l};
      if(d == -1) return {0, r - l};
    }
    int l0 = (!l ? 0 : bv[b].rank0(l)), l1 = l - l0;
    int r0 = bv[b].rank0(r), r1 = r - r0;
    auto R = __range_kth_largest_super(l1, r1, s, t, (S + T) >> 1, T, k, b * 2 + 2, d - 1);
    if(R.first != -1){
      R.first += 1 << d;
      return R;
    }
    auto L = __range_kth_largest_super(l0, r0, s, t, S, (S + T) >> 1, k - R.second, b * 2 + 1, d - 1);
    if(L.first == -1) L.second += R.second;
    return L;
  }
  */
public:
  wavelet_tree_monoid(){}
  wavelet_tree_monoid(const std::vector<std::pair<int, Val>> &v, int inf): n(v.size()), inf(inf){
    assert(inf >= 0);
    h = 0;
    while((1 << h) < inf) h++;
    std::vector<int> v2(n);
    std::vector<Val> v3(n);
    for(int i = 0; i < n; i++) std::tie(v2[i], v3[i]) = v[i];
    build(v2, v3);
  }
  wavelet_tree_monoid(const std::vector<int> &v, const std::vector<Val> &val, int inf): n(v.size()), inf(inf){
    assert(inf >= 0);
    assert(v.size() == val.size());
    h = 0;
    while((1 << h) < inf) h++;
    build(v, val);
  }
  int access(int k){
    assert(0 <= k && k < n);
    return __access(k, 0, h - 1);
  }
  // [0, r)のcの数
  int rank(int r, int c){
    assert(0 <= r && r <= n);
    assert(0 <= c && c < inf);
    return __rank(r, c, 0, h - 1);
  }
  // k番目のc, ない場合は-1
  int select(int k, int c){
    assert(0 <= k && k < n);
    assert(0 <= c && c < inf);
    int left_leaf = (1 << h) - 1;
    return __select(k, c, (left_leaf + c - 1) >> 1, 0);
  }
  // k以降(k含む)のc, 無い場合は-1
  int find_next(int k, int c){
    if(c < 0 || c >= inf) return -1;
    return select(rank(k, c), c);
  }
  // k以前(k含む)のc, 無い場合は-1
  int find_prev(int k, int c){
    if(c < 0 || c >= inf) return -1;
    int r = rank(k + 1, c);
    if(r == 0) return -1;
    return select(r - 1, c);
  }
  // [l, r)の[s, t)の数
  int range_freq(int l, int r, int s, int t){
    assert(0 <= l && r <= n);
    assert(0 <= s && t <= inf);
    if(l >= r || s >= t) return 0;
    return __range_freq(l, r, 0, s, t, 0, 1 << h);
  }
  // [l, r)でk番目に小さい要素(ない場合は-1)
  int range_kth_smallest(int l, int r, int k){
    assert(0 <= l && r <= n);
    if(k >= r - l) return -1;
    return __range_kth_smallest(l, r, k, 0, h - 1);
  }
  // [l, r)でk番目に大きい要素(ない場合は-1)
  int range_kth_largest(int l, int r, int k){
    assert(0 <= l && r <= n);
    if(k >= r - l) return -1;
    return __range_kth_smallest(l, r, r - l - 1 - k, 0, h - 1);
  }
  // [l, r)で値が[s, t)の要素でk番目に大きい要素とその個数(ない場合は-1)
  std::pair<int, int> range_kth_smallest_super(int l, int r, int s, int t, int k){
    assert(0 <= l && r <= n);
    assert(0 <= s && t <= inf);
    return __range_kth_smallest_super(l, r, s, t, 0, 1 << h, k, 0, h - 1);
  }
  // [l, r)で値が[s, t)の要素でk番目に大きい要素とその個数(ない場合は-1)
  std::pair<int, int> range_kth_largest_super(int l, int r, int s, int t, int k){
    assert(0 <= l && r <= n);
    assert(0 <= s && t <= inf);
    return __range_kth_largest_super(l, r, s, t, 0, 1 << h, k, 0, h - 1);
  }
  // 区間[l, r)のx以上の最小要素(ない場合は-1)
  int range_lower_bound(int l, int r, int x){
    assert(0 <= l && r <= n);
    assert(0 <= x && x < inf);
    int cnt = range_freq(l, r, 0, x);
    return range_kth_smallest(l, r, cnt).first;
  }
  // 区間[l, r)のx以下の最大要素(ない場合は-1)
  int range_lower_bound_rev(int l, int r, int x){
    assert(0 <= l && r <= n);
    assert(0 <= x && x < inf);
    int cnt = range_freq(l, r, 0, x + 1);
    if(cnt == 0) return -1;
    return range_kth_smallest(l, r, cnt - 1).first;
  }
  // [l, r)で値が[s, t)の要素でインデックス順にk番目, ない場合は-1
  int range_select(int l, int r, int s, int t, int k){
    if(range_freq(l, r, s, t) <= k) return -1;
    int L = l;
    while(r - l > 1){
      int mid_r = (l + r) / 2;
      if(range_freq(L, mid_r, s, t) > k) r = mid_r;
      else l = mid_r;
    }
    return r - 1;
  }
  // ret[i] := [l, r)でh-bitの2進数で見たときxとのlcpがちょうどiの要素の数
  std::vector<int> count_prefix(int l, int r, int x){
    std::vector<int> ret(h + 1, 0);
    __count_prefix(l, r, x, ret);
    return ret;
  }
  // val[k] <- z
  void set(int k, Val z){
    __set(k, 0, h - 1, z);
  }
  // val[k] <- merge(val[k], z)
  void update(int k, Val z){
    __update(k, 0, h - 1, z);
  }
  Val query(int l, int r, int s, int t){
    assert(0 <= l && r <= n);
    assert(0 <= s && t <= inf);
    if(l >= r || s >= t) return monoid::id();
    return __query(l, r, 0, s, t, 0, 1 << h);
  }
  /*
  template<typename F>
  int bisect_from_left(int l, int r, int s, F &f){
    assert(0 <= l && r <= n);
    assert(0 <= s && s < inf);

  }
  */
};

template<typename T, typename monoid>
struct compressed_wavelet_tree_monoid{
  using Val = typename monoid::Val;
private:
  int n;
  std::vector<T> rev;
  wavelet_tree_monoid<monoid> wt;
  int lb(T c){
    return std::lower_bound(rev.begin(), rev.end(), c) - rev.begin();
  }
public:
  compressed_wavelet_tree_monoid(const std::vector<T> &_v, const std::vector<Val> &_val): n(_v.size()){
    assert(_val.size() == n);
    std::vector<std::pair<T, int>> tmp(n);
    std::vector<int> v(n);
    for(int i = 0; i < n; i++) tmp[i].first = _v[i], tmp[i].second = i;
    std::sort(tmp.begin(), tmp.end());
    for(int i = 0; i < n; i++){
      if(i == 0 || tmp[i].first != tmp[i - 1].first) rev.push_back(tmp[i].first);
      v[tmp[i].second] = (int)rev.size() - 1;
    }
    tmp.clear();
    wt = wavelet_tree_monoid<monoid>(v, _val, rev.size());
  }
  compressed_wavelet_tree_monoid(const std::vector<std::pair<T, Val>> &_v): n(_v.size()){
    std::vector<std::pair<T, int>> tmp(n);
    std::vector<int> v(n);
    std::vector<Val> val(n);
    for(int i = 0; i < n; i++) tmp[i].first = _v[i].first, tmp[i].second = i;
    std::sort(tmp.begin(), tmp.end());
    for(int i = 0; i < n; i++){
      if(i == 0 || tmp[i].first != tmp[i - 1].first) rev.push_back(tmp[i].first);
      v[tmp[i].second] = (int)rev.size() - 1;
    }
    tmp.clear();
    for(int i = 0; i < n; i++) val[i] = _v[i].second;
    wt = wavelet_tree_monoid<monoid>(v, val, rev.size());
  }
  // V[k]
  T access(int k){
    return rev[wt.access(k)];
  }
  // r未満のcの数
  int rank(int r, T c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return 0;
    return wt.rank(r, idx);
  }
  // k番目のc, 無い場合は-1
  int select(int k, T c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return -1;
    return wt.select(k, idx);
  }
  // k以降(k含む)のc, 無い場合は-1
  int find_next(int k, T c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return -1;
    return wt.find_next(k, idx);
  }
  // k以前(k含む)のc, 無い場合は-1
  int find_prev(int k, T c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return -1;
    return wt.find_prev(k, idx);
  }
  // [l, r)の[s, t)の個数
  int range_freq(int l, int r, T s, T t){
    return wt.range_freq(l, r, lb(s), lb(t));
  }
  // [l, r)でk番目に小さい要素とインデックス(ない場合は-1)
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<T, int> range_kth_smallest(int l, int r, int k){
    auto p = wt.range_kth_smallest(l, r, k);
    return {rev[p.first], p.second};
  }
  // [l, r)でk番目に大きい要素とインデックス(ない場合は-1)
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<T, int> range_kth_largest(int l, int r, int k){
    auto p = wt.range_kth_largest(l, r, k);
    return {rev[p.first], p.second};
  }
  // [l, r)かつ[s, t)でk番目に小さい要素とその個数
  std::pair<T, int> range_kth_smallest_super(int l, int r, T s, T t, int k){
    auto x = wt.range_kth_smallest_super(l, r, lb(s), lb(t), k);
    assert(x.second != -1);
    return {rev[x.first], x.second};
  }
  // [l, r)かつ[s, t)でk番目に大きい要素とその個数
  std::pair<T, int> range_kth_largest_super(int l, int r, T s, T t, int k){
    auto x = wt.range_kth_largest_super(l, r, lb(s), lb(t), k);
    assert(x.second != -1);
    return {rev[x.first], x.second};
  }
  // [l, r)かつ[s, t)でインデックス順にk番目
  int range_select(int l, int r, T s, T t, int k){
    return wt.range_select(l, r, lb(s), lb(t), k);
  }
  // val[k] <- z
  void set(int k, Val z){
    wt.set(k, z);
  }
  // val[k] <- merge(val[k], z)
  void update(int k, Val z){
    wt.update(k, z);
  }
  // [l, r)の[ly, ry)のモノイド積
  Val query(int l, int r, T ly, T ry){
    return wt.query(l, r, lb(ly), lb(ry));
  }
};

template<typename Idx, typename monoid>
struct compressed_rectangle_wavelet_tree_monoid{
  using Val = typename monoid::Val;
private:
  int n;
  std::vector<Idx> X, Y;
  std::vector<std::pair<Idx, Idx>> XY;
  wavelet_tree_monoid<monoid> wt;
  int lb_x(Idx c){
    return std::lower_bound(X.begin(), X.end(), c) - X.begin();
  }
  int lb_y(Idx c){
    return std::lower_bound(Y.begin(), Y.end(), c) - Y.begin();
  }
public:
  compressed_rectangle_wavelet_tree_monoid(std::vector<std::tuple<Idx, Idx, Val>> _v): n(_v.size()){
    sort(_v.begin(), _v.end());
    for(int i = 0; i < n; i++){
      X.push_back(std::get<0>(_v[i]));
      Y.push_back(std::get<1>(_v[i]));
    }
    sort(Y.begin(), Y.end());
    Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
    std::vector<std::pair<int, Val>> v2(n);
    XY.resize(n);
    for(int i = 0; i < n; i++){
      XY[i] = {std::get<0>(_v[i]), std::get<1>(_v[i])};
      v2[i] = {lb_y(std::get<1>(_v[i])), std::get<2>(_v[i])};
    }
    wt = wavelet_tree_monoid<monoid>(v2, Y.size());
  }
  // (x, y)のインデックスを探す, ない場合は-1
  int find_point(Idx x, Idx y){
    int k = std::lower_bound(XY.begin(), XY.end(), std::make_pair(x, y)) - XY.begin();
    if(k == n || XY[k].first != x || XY[k].second != y) return -1;
    return k;
  }
  // [lx, rx)かつ[ly, ry)の数
  int range_freq(Idx lx, Idx rx, Idx ly, Idx ry){
    return wt.range_freq(lb_x(lx), lb_x(rx), lb_y(ly), lb_y(ry));
  }
  // val{x, y} <- z {x, y}がないとエラー
  void set(Idx x, Idx y, Val z){
    int k = find_point(x, y);
    assert(k != -1);
    wt.set(k, z);
  }
  // val{x, y} <- merge(val{x, y}, z) {x, y}がないとエラー
  void update(Idx x, Idx y, Val z){
    int k = find_point(x, y);
    assert(k != -1);
    wt.update(k, z);
  }
  Val query(Idx lx, Idx rx, Idx ly, Idx ry){
    return wt.query(lb_x(lx), lb_x(rx), lb_y(ly), lb_y(ry));
  }
};
#endif