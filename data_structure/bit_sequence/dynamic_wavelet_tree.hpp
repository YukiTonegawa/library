#ifndef _DYNAMIC_WAVELET_MATRIX_H_
#define _DYNAMIC_WAVELET_MATRIX_H_
#include "dynamic_bitvector.hpp"
#include <vector>
#include <bitset>
#include <cassert>
#include <algorithm>
#include <array>
#include <queue>
#include <numeric>
#include <cstdint>

struct dynamic_wavelet_tree{
private:
  using __bitvector = dynamic_bitvector;
  using node = __bitvector::node;
  int n, h, inf;
  std::vector<__bitvector> bv;
  void build(std::vector<int> v){
    bv.resize((1 << h) - 1);
    std::queue<std::array<int, 4>> q;// (dep, v上でのl, r)
    std::vector<int> tmp(n);
    if(h) q.push({h - 1, 0, n, 0});
    while(!q.empty()){
      auto [d, l, r, k] = q.front();
      q.pop();
      std::vector<bool> bits(r - l);
      int lidx = 0, ridx = 0;
      for(int i = l; i < r; i++){
        bool b = (v[i] >> d) & 1;
        bits[i - l] = b;
        if(b) tmp[ridx++] = v[i];
        else v[l + lidx++] = v[i];
      }
      for(int i = 0; i < ridx; i++) v[l + lidx + i] = tmp[i];
      bv[k] = __bitvector(bits);
      if(d){
        int mid = l + lidx;
        q.push({d - 1, l, mid, k * 2 + 1});
        q.push({d - 1, mid, r, k * 2 + 2});
      }
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
  int __select_erase(int k, int x, int b, int d){
    bool f = (x >> d) & 1;
    int nxt_k = f ? bv[b].select1(k) : bv[b].select0(k);
    assert(nxt_k != -1);
    bv[b].erase(nxt_k);
    return !b ? nxt_k : __select_erase(nxt_k, x, (b - 1) >> 1, d + 1);
  }
  int __access(int k, int b, int d){
    int f = bv[b].get(k);
    int nxt_k = f ? bv[b].rank1(k) : bv[b].rank0(k);
    return (f << d) + (d ? __access(nxt_k, b * 2 + 1 + f, d - 1) : 0);
  }
  int __find_next(int k, int x, int b, int d){
    int f = (x >> d) & 1;
    int nxt_k = f ? bv[b].rank1(k) : bv[b].rank0(k);
    return (f << d) + (d ? __access(nxt_k, b * 2 + 1 + f, d - 1) : 0);
  }
  void __insert(int k, int x, int b, int d){
    bool f = (x >> d) & 1;
    int nxt_k = f ? bv[b].rank1(k) : bv[b].rank0(k);
    bv[b].insert(k, f);
    if(d) __insert(nxt_k, x, b * 2 + 1 + f, d - 1);
  }
  int __erase(int k, int b, int d){
    int f = bv[b].erase(k);
    int nxt_k = f ? bv[b].rank1(k) : bv[b].rank0(k);
    return (f << d) + (d ? __erase(nxt_k, b * 2 + 1 + f, d - 1) : 0);
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
public:
  dynamic_wavelet_tree(){}
  dynamic_wavelet_tree(const std::vector<int> &v, int inf): n(v.size()), inf(inf){
    assert(inf >= 0);
    h = 0;
    while((1 << h) < inf) h++;
    build(v);
  }
  int size(){
    return n;
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
  // k番目のcの位置を返してそれを削除する
  // k < count(c)でない状態で使うと壊れる
  int select_erase(int k, int c){
    assert(0 <= k && k < n);
    assert(0 <= c && c < inf);
    int left_leaf = (1 << h) - 1;
    int ret = __select_erase(k, c, (left_leaf + c - 1) >> 1, 0);
    n = bv[0].size();
    return ret;
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
  // k番目にcを挿入
  void insert(int k, int c){
    assert(0 <= k && k <= n);
    assert(0 <= c && c < inf);
    __insert(k, c, 0, h - 1);
    n = bv[0].size();
  }
  // k番目を削除してその値を返す
  int erase(int k){
    assert(0 <= k && k < n);
    int ret = __erase(k, 0, h - 1);
    n = bv[0].size();
    return ret;
  }
  // k番目の値をcにする
  void update(int k, int c){
    erase(k);
    insert(k, c);
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
    return range_kth_smallest(l, r, cnt);
  }
  // 区間[l, r)のx以下の最大要素(ない場合は-1)
  int range_lower_bound_rev(int l, int r, int x){
    assert(0 <= l && r <= n);
    assert(0 <= x && x < inf);
    int cnt = range_freq(l, r, 0, x + 1);
    if(cnt == 0) return -1;
    return range_kth_smallest(l, r, cnt - 1);
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
};

template<typename T>
struct compressed_dynamic_wavelet_tree{
private:
  std::vector<T> rev;
  dynamic_wavelet_tree wt;
  int lb(T c){
    return std::lower_bound(rev.begin(), rev.end(), c) - rev.begin();
  }
public:
  compressed_dynamic_wavelet_tree(){}
  // 初期状態の配列
  compressed_dynamic_wavelet_tree(const std::vector<T> &_v){
    int n = _v.size();
    std::vector<int> v(n);
    std::vector<std::pair<T, int>> tmp(n);
    for(int i = 0; i < n; i++) tmp[i].first = _v[i], tmp[i].second = i;
    std::sort(tmp.begin(), tmp.end());
    for(int i = 0; i < n; i++){
      if(i == 0 || tmp[i].first != tmp[i - 1].first) rev.push_back(tmp[i].first);
      v[tmp[i].second] = (int)rev.size() - 1;
    }
    tmp.clear();
    wt = dynamic_wavelet_tree(v, rev.size());
  }
  // {初期状態の配列, 後から追加されるかもしれない値}
  compressed_dynamic_wavelet_tree(const std::vector<T> &_v, const std::vector<T> &p){
    int n = _v.size();
    for(int i = 0; i < n; i++) rev.push_back(_v[i]);
    for(int i = 0; i < p.size(); i++) rev.push_back(p[i]);
    std::sort(rev.begin(), rev.end());
    rev.erase(std::unique(rev.begin(), rev.end()), rev.end());
    std::vector<int> v(n);
    for(int i = 0; i < n; i++) v[i] = lb(_v[i]);
    wt = dynamic_wavelet_tree(v, rev.size());
  }
  // V[k]
  T access(int k){
    assert(k < wt.size());
    return rev[wt.access(k)];
  }
  // インデックスがr未満のcの数
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
  // k番目のcの位置を返してそれを削除する
  // k < count(c)でない, またはcが事前に登録されていない状態で使うと壊れる
  int select_erase(int k, T c){
    int idx = lb(c);
    assert(idx != rev.size() && rev[idx] == c);
    return wt.select_erase(k, idx);
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
  // k番目にcを挿入(cが事前に登録していない値だと壊れる)
  void insert(int k, T c){
    int idx = lb(c);
    assert(idx != rev.size() && rev[idx] == c);
    wt.insert(k, idx);
  }
  // k番目を削除してその値を返す
  T erase(int k){
    int ret = wt.erase(k);
    return rev[ret];
  }
  // k番目の値をcにする
  void update(int k, T c){
    erase(k);
    insert(k, c);
  }
  // [l, r)の[s, t)の個数
  int range_freq(int l, int r, T s, T t){
    return wt.range_freq(l, r, lb(s), lb(t));
  }
  // [l, r)でk番目に小さい要素(ない場合は-1)
  T range_kth_smallest(int l, int r, int k){
    int ret = wt.range_kth_smallest(l, r, k);
    if(ret == -1) return ret;
    return rev[ret];
  }
  // [l, r)でk番目に大きい要素(ない場合は-1)
  T range_kth_largest(int l, int r, int k){
    int ret = wt.range_kth_largest(l, r, k);
    if(ret == -1) return ret;
    return rev[ret];
  }
  // [l, r)で値が[s, t)の要素でk番目に大きい要素とその個数(ない場合は-1)
  std::pair<T, int> range_kth_smallest_super(int l, int r, T s, T t, int k){
    auto p = wt.range_kth_smallest_super(l, r, lb(s), lb(t), k);
    if(p.first == -1) return {-1, -1};
    return {rev[p.first], p.second};
  }
  // [l, r)で値が[s, t)の要素でk番目に大きい要素とその個数(ない場合は-1)
  std::pair<T, int> range_kth_largest_super(int l, int r, T s, T t, int k){
    auto p = wt.range_kth_largest_super(l, r, lb(s), lb(t), k);
    if(p.first == -1) return {-1, -1};
    return {rev[p.first], p.second};
  }
  // 区間[l, r)のx以上の最小要素(ない場合は-1)
  T range_lower_bound(int l, int r, T x){
    int xi = lb(x);
    if(xi == rev.size()) return -1;
    int ret = wt.range_lower_bound(l, r, xi);
    return ret == -1 ? -1 : rev[ret];
  }
  // 区間[l, r)のx以下の最大要素(ない場合は-1)
  T range_lower_bound_rev(int l, int r, T x){
    int xi = std::upper_bound(rev.begin(), rev.end(), x) - rev.begin();
    if(xi == 0) return -1;
    int ret = wt.range_lower_bound_rev(l, r, xi - 1);
    return ret == -1 ? -1 : rev[ret];
  }
  // [l, r)で値が[s, t)の要素でインデックス順にk番目, ない場合は-1
  int range_select(int l, int r, T s, T t, int k){
    return wt.range_select(l, r, lb(s), lb(t), k);
  }
};
#endif
