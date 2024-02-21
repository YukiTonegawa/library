#ifndef _WAVELET_MATRIX_H_
#define _WAVELET_MATRIX_H_
#include "bitvector.hpp"
#include <vector>
#include <queue>
#include <tuple>
#include <numeric>

struct wavelet_matrix{
private:
  using __bitvector = bitvector_memory;
  int n, h, inf;
  std::vector<__bitvector> bv;
  std::vector<int> bottom_idx;
  void build(std::vector<int> v){
    bv.resize(h);
    std::vector<bool> bits(n);
    std::vector<int> tmp(n), tmp_idx(n);
    std::iota(bottom_idx.begin(), bottom_idx.end(), 0);
    std::iota(tmp_idx.begin(), tmp_idx.end(), 0);
    std::queue<std::tuple<int, int, int>> q;
    if(h) q.push({h - 1, 0, n});
    while(!q.empty()){
      auto [d, l, r] = q.front();
      q.pop();
      int lidx = 0, ridx = 0;
      for(int i = l; i < r; i++){
        bool b = (v[i] >> d) & 1;
        bits[i] = b;
        if(b) tmp[ridx] = v[i], tmp_idx[ridx++] = bottom_idx[i];
        else v[l + lidx] = v[i], bottom_idx[l + lidx++] = bottom_idx[i];
      }
      for(int i = 0; i < ridx; i++){
        v[l + lidx + i] = tmp[i];
        bottom_idx[l + lidx + i] = tmp_idx[i];
      }
      if(d){
        int mid = l + lidx;
        if(mid > l) q.push({d - 1, l, mid});
        if(r > mid) q.push({d - 1, mid, r});
      }
      if(r == n) bv[d] = __bitvector(bits);
    }
  }
  // v[k]
  int __access(int k){
    assert(0 <= k && k < n);
    int L = 0, R = n, ret = 0;
    for(int i = h - 1; i >= 0; i--){
      int L0 = bv[i].rank0(L), R0 = bv[i].rank0(R), k0 = bv[i].rank0(k);
      if(bv[i].access(k)){
        k += R0 - k0;
        L += R0 - L0;
        ret += 1 << i;
      }else{
        k -= (k - L) - (k0 - L0);
        R = L + R0 - L0;
      }
    }
    return ret;
  }
  // [0, r)のcの数
  int __rank(int r, int c){
    if(c < 0 || c >= inf) return 0;
    int L = 0, R = n;
    for(int i = h - 1; i >= 0 && r > L; i--){
      int L0 = bv[i].rank0(L), R0 = bv[i].rank0(R) - L0, r0 = bv[i].rank0(r) - L0;
      if((c >> i) & 1){
        r += R0 - r0;
        L += R0;
      }else{
        r -= (r - L) - r0;
        R = L + R0;
      }
    }
    return r - L;
  }
  // k番目のc, ない場合は-1
  int __select(int k, int c){
    if(c < 0 || c >= inf) return -1;
    int L = 0, R = n;
    for(int i = h - 1; i >= 0; i--){
      int LR0 = bv[i].rank0(R) - bv[i].rank0(L);
      if((c >> i) & 1){
        L += LR0;
      }else{
        R = L + LR0;
      }
    }
    if(R - L <= k) return -1;
    return bottom_idx[L + k];
  }
  int __range_freq(int d, int l, int r, int L, int R, int s, int t, int S, int T){
    if(l == r || L == R || t <= S || T <= s) return 0;
    else if(s <= S && T <= t) return r - l;
    int L0 = bv[d].rank0(L), R0 = bv[d].rank0(R) - L0;
    int l0 = bv[d].rank0(l) - L0, r0 = bv[d].rank0(r) - L0;
    return __range_freq(d - 1, L + l0, L + r0, L, L + R0, s, t, S, (S + T) / 2) +
            __range_freq(d - 1, l + (R0 - l0), r + (R0 - r0), L + R0, R, s, t, (S + T) / 2, T);
  }
  std::pair<int, int> __range_kth_smallest(int l, int r, int k){
    if(l >= r || r - l <= k) return {-1, -1};
    int L = 0, R = n, ret = 0;
    for(int i = h - 1; i >= 0; i--){
      int L0 = bv[i].rank0(L), R0 = bv[i].rank0(R) - L0;
      int l0 = bv[i].rank0(l) - L0, r0 = bv[i].rank0(r) - L0;
      if(r0 - l0 <= k){
        k -= r0 - l0;
        ret += 1 << i;
        l += R0 - l0, r += R0 - r0, L += R0;
      }else{
        l = L + l0, r = L + r0, R = L + R0;
      }
    }
    return {ret, bottom_idx[l + k]};
  }
  std::pair<int, int> __range_kth_smallest_super(int d, int l, int r, int L, int R, int s, int t, int S, int T, int &k){
    if(l == r || L == R || T <= s || t <= S) return {-1, -1};
    else if(s <= S && T <= t){
      if(r - l > k){
        for(int i = d; i >= 0; i--){
          int L0 = bv[i].rank0(L), R0 = bv[i].rank0(R) - L0;
          int l0 = bv[i].rank0(l) - L0, r0 = bv[i].rank0(r) - L0;
          if(r0 - l0 <= k){
            k -= r0 - l0;
            l += R0 - l0, r += R0 - r0, L += R0;
            S = (S + T) / 2;
          }else{
            l = L + l0, r = L + r0, R = L + R0;
            T = (S + T) / 2;
          }
        }
        return {S, bottom_idx[l + k]};
      }
      k -= r - l;
      return {-1, -1};
    }
    int L0 = bv[d].rank0(L), R0 = bv[d].rank0(R) - L0;
    int l0 = bv[d].rank0(l) - L0, r0 = bv[d].rank0(r) - L0;
    auto p = __range_kth_smallest_super(d - 1, L + l0, L + r0, L, L + R0, s, t, S, (S + T) / 2, k);
    if(p.first != -1) return p;
    return __range_kth_smallest_super(d - 1, l + (R0 - l0), r + (R0 - r0), L + R0, R, s, t, (S + T) / 2, T, k);
  }
  std::pair<int, int> __range_kth_largest_super(int d, int l, int r, int L, int R, int s, int t, int S, int T, int &k){
    if(l == r || L == R || T <= s || t <= S) return {-1, -1};
    else if(s <= S && T <= t){
      if(r - l > k){
        for(int i = d; i >= 0; i--){
          int L0 = bv[i].rank0(L), R0 = bv[i].rank0(R) - L0;
          int l0 = bv[i].rank0(l) - L0, r0 = bv[i].rank0(r) - L0;
          if((r - l) - (r0 - l0) <= k){
            k -= (r - l) - (r0 - l0);
            l = L + l0, r = L + r0, R = L + R0;
            T = (S + T) / 2;
          }else{
            l += R0 - l0, r += R0 - r0, L += R0;
            S = (S + T) / 2;
          }
        }
        return {S, bottom_idx[r - 1 - k]};
      }
      k -= r - l;
      return {-1, -1};
    }
    int L0 = bv[d].rank0(L), R0 = bv[d].rank0(R) - L0;
    int l0 = bv[d].rank0(l) - L0, r0 = bv[d].rank0(r) - L0;
    auto p =  __range_kth_largest_super(d - 1, l + (R0 - l0), r + (R0 - r0), L + R0, R, s, t, (S + T) / 2, T, k);
    if(p.first != -1) return p;
    return  __range_kth_largest_super(d - 1, L + l0, L + r0, L, L + R0, s, t, S, (S + T) / 2, k);
  }
  void __count_prefix(int l, int r, int x, std::vector<int> &ret){
    int L = 0, R = n;
    for(int i = h - 1; i >= 0; i--){
      int L0 = bv[i].rank0(L), R0 = bv[i].rank0(R) - L0;
      int l0 = bv[i].rank0(l) - L0, r0 = bv[i].rank0(r) - L0;
      if((x >> i) & 1){
        ret[h - 1 - i] = r0 - l0;
        l += R0 - l0, r += R0 - r0, L += R0;
      }else{
        ret[h - 1 - i] = (r - l) - (r0 - l0);
        l = L + l0, r = L + r0, R = L + R0;
      }
    }
    ret[h] = r - l;
  }
public:
  wavelet_matrix(): n(0){}
  // 値が[0, inf)
  wavelet_matrix(const std::vector<int> &v, int inf): n(v.size()), inf(inf), bottom_idx(n){
    assert(inf >= 0);
    h = 0;
    while((1 << h) < inf) h++;
    build(v);
  }
  // v[k]
  int access(int k){
    return __access(k);
  }
  // [0, r)のcの数
  int rank(int r, int c){
    return __rank(r, c);
  }
  // k番目のc, ない場合は-1
  int select(int k, int c){
    return __select(k, c);
  }
  // k以降(k含む)のc, 無い場合は-1
  int find_next(int k, int c){
    if(c < 0 || c >= inf) return -1;
    return __select(__rank(k, c), c);
  }
  // k以前(k含む)のc, 無い場合は-1
  int find_prev(int k, int c){
    if(c < 0 || c >= inf) return -1;
    int r = __rank(k + 1, c);
    if(r == 0) return -1;
    return __select(r - 1, c);
  }
  // [l, r)の[s, t)の数
  int range_freq(int l, int r, int s, int t){
    assert(0 <= l && r <= n);
    assert(0 <= s && t <= inf);
    if(l >= r || s >= t) return 0;
    return __range_freq(h - 1, l, r, 0, n, s, t, 0, 1 << h);
  }
  // [l, r)でk番目に小さい要素とインデックス(ない場合は-1)
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<int, int> range_kth_smallest(int l, int r, int k){
    return __range_kth_smallest(l, r, k);
  }
  // [l, r)でk番目に大きい要素とインデックス(ない場合は-1)
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<int, int> range_kth_largest(int l, int r, int k){
    if(r - l <= k) return {-1, -1};
    return __range_kth_smallest(l, r, r - l - 1 - k);
  }
  // [l, r)で値が[s, t)の要素でk番目に小さい
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<int, int> range_kth_smallest_super(int l, int r, int s, int t, int k){
    return __range_kth_smallest_super(h - 1, l, r, 0, n, s, t, 0, 1 << h, k);
  }
  // [l, r)で値が[s, t)の要素でk番目に大きい
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<int, int> range_kth_largest_super(int l, int r, int s, int t, int k){
    return __range_kth_largest_super(h - 1, l, r, 0, n, s, t, 0, 1 << h, k);
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
    return l;
  }
  // ret[i] := [l, r)でh-bitの2進数で見たときxとのlcpがちょうどiの要素の数
  std::vector<int> count_prefix(int l, int r, int x){
    std::vector<int> ret(h + 1, 0);
    __count_prefix(l, r, x, ret);
    return ret;
  }
};

template<typename T>
struct compressed_wavelet_matrix{
private:
  std::vector<T> rev;
  wavelet_matrix wm;
  int lb(T c){
    return std::lower_bound(rev.begin(), rev.end(), c) - rev.begin();
  }
public:
  compressed_wavelet_matrix(){}
  compressed_wavelet_matrix(const std::vector<T> &_v){
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
    wm = wavelet_matrix(v, rev.size());
  }
  // V[k]
  T access(int k){
    return rev[wm.access(k)];
  }
  // インデックスがr未満のcの数
  int rank(int r, T c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return 0;
    return wm.rank(r, idx);
  }
  // k番目のc, 無い場合は-1
  int select(int k, T c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return -1;
    return wm.select(k, idx);
  }
  // k以降(k含む)のc, 無い場合は-1
  int find_next(int k, T c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return -1;
    return wm.find_next(k, idx);
  }
  // k以前(k含む)のc, 無い場合は-1
  int find_prev(int k, T c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return -1;
    return wm.find_prev(k, idx);
  }
  // [l, r)の[s, t)の個数
  int range_freq(int l, int r, T s, T t){
    return wm.range_freq(l, r, lb(s), lb(t));
  }
  // [l, r)でk番目に小さい要素と場所, ない場合は-1
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<T, int> range_kth_smallest(int l, int r, int k){
    auto p = wm.range_kth_smallest(l, r, k);
    if(p.first == -1) return {-1, -1};
    return {rev[p.first], p.second};
  }
  // [l, r)でk番目に大きい要素と場所, ない場合は-1
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<T, int> range_kth_largest(int l, int r, int k){
    auto p = wm.range_kth_largest(l, r, k);
    if(p.first == -1) return {-1, -1};
    return {rev[p.first], p.second};
  }
  // [l, r)の[s, t)の中でk番目に小さい要素と場所, ない場合は-1
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<T, int> range_kth_smallest_super(int l, int r, T s, T t, int k){
    auto p = wm.range_kth_smallest_super(l, r, lb(s), lb(t), k);
    if(p.first == -1) return {-1, -1};
    return {rev[p.first], p.second};
  }
  // [l, r)の[s, t)の中でk番目に大きい要素と場所, ない場合は-1
  // 値が同じ場合インデックスが小さいものを小さいとする
  std::pair<T, int> range_kth_largest_super(int l, int r, T s, T t, int k){
    auto p = wm.range_kth_largest_super(l, r, lb(s), lb(t), k);
    if(p.first == -1) return {-1, -1};
    return {rev[p.first], p.second};
  }
  // 区間[l, r)のx以上の最小要素(ない場合は-1)
  T range_lower_bound(int l, int r, T x){
    int xi = lb(x);
    if(xi == rev.size()) return -1;
    int ret = wm.range_lower_bound(l, r, xi);
    return ret == -1 ? -1 : rev[ret];
  }
  // 区間[l, r)のx以下の最大要素(ない場合は-1)
  T range_lower_bound_rev(int l, int r, T x){
    int xi = std::upper_bound(rev.begin(), rev.end(), x) - rev.begin();
    if(xi == 0) return -1;
    int ret = wm.range_lower_bound_rev(l, r, xi - 1);
    return ret == -1 ? -1 : rev[ret];
  }
  // [l, r)で値が[s, t)の要素でインデックス順にk番目, ない場合は-1
  int range_select(int l, int r, T s, T t, int k){
    return wm.range_select(l, r, lb(s), lb(t), k);
  }
};

template<int R>
struct wavelet_matrix_arbitrary_radix{
private:
  static constexpr int bitlen(uint64_t x){
    if(!x) return 0;
    return 64 - __builtin_clzll(x);
  }
  static constexpr int Rdiv = bitlen(R) - 1, Rmod = R - 1;
  int n, h, inf;
  std::vector<bitvector_arbitrary_radix<R>> bv;
  std::vector<int> bottom_idx;
  void build(const std::vector<int> &v){
    bv.resize(h);
    bottom_idx.resize(n);
    std::queue<std::tuple<int, int, int>> q;
    if(h) q.push({h - 1, 0, n});
    std::vector<std::pair<int, int>> cur(n);
    std::array<std::vector<std::pair<int, int>>, R> next;
    std::vector<uint8_t> bits(n);
    for(int i = 0; i < n; i++) cur[i] = {v[i], i};
    while(!q.empty()){
      auto [d, l, r] = q.front();
      q.pop();
      for(int i = l; i < r; i++){
        int dir = (cur[i].first >> (d * Rdiv)) & Rmod;
        bits[i] = dir;
        next[dir].push_back({cur[i].first, cur[i].second});
      }
      for(int i = 0, j = l; i < R; i++){
        std::copy(next[i].begin(), next[i].end(), cur.begin() + j);
        int k = j + next[i].size();
        if(d && j < k) q.push({d - 1, j, k});
        j = k;
        next[i].clear();
      }
      if(r == n) bv[d] = bitvector_arbitrary_radix<R>(bits);
    }
    for(int i = 0; i < n; i++) bottom_idx[i] = cur[i].second;
  }
  int __rank(int r, int c){
    int l = 0, _L = 0, _R = n;
    for(int d = h - 1, shift = d * Rdiv; r > l && d >= 0; d--, shift -= Rdiv){
      int dir = (c >> shift) & Rmod;
      int Ldir = bv[d].rank_lower(_L, dir);
      int Ldir2 = bv[d].rank(_L, dir);
      int Rdir = bv[d].rank_lower(_R, dir);
      _L += Rdir - Ldir;
      _R = _L + bv[d].rank(_R, dir) - Ldir2;
      l = _L + bv[d].rank(l, dir) - Ldir2;
      r = _L + bv[d].rank(r, dir) - Ldir2;
    }
    return r - l;
  }
  int __select(int k, int c){
    int _L = 0, _R = n;
    for(int d = h - 1, shift = d * Rdiv; d >= 0; d--, shift -= Rdiv){
      int dir = (c >> shift) & Rmod;
      int Ldir = bv[d].rank_lower(_L, dir), Ldir2 = bv[d].rank(_L, dir), Rdir = bv[d].rank_lower(_R, dir);
      _L += Rdir - Ldir, _R = _L + bv[d].rank(_R, dir) - Ldir2;
      if(_L + k >= _R) return -1;
    }
    return bottom_idx[_L + k];
  }
  int __range_freq(int d, int l, int r, int _L, int _R, int s, int t, int S, int T){
    if(l >= r || s >= t || t <= S || T <= s) return 0;
    else if(s <= S && T <= t) return r - l;
    if(s < S) s = S;
    if(t > T) t = T;
    int shift = d * Rdiv;
    int sb = (s - S) >> shift, tb = (t - S) >> shift; // [0, R]
    if(sb == tb){
      int Lsb = bv[d].rank_lower(_L, sb), Lsb2 = bv[d].rank(_L, sb), Rsb = bv[d].rank_lower(_R, sb);
      int Lnext = _L + Rsb - Lsb, Rnext = Lnext + bv[d].rank(_R, sb) - Lsb2;
      int lnext = Lnext + bv[d].rank(l, sb) - Lsb2, rnext = Lnext + bv[d].rank(r, sb) - Lsb2;
      return __range_freq(d - 1, lnext, rnext, Lnext, Rnext, s, t, S + (sb << shift), S + ((sb + 1) << shift));
    }
    // f(tb) - f(sb) + [tb] - [sb]
    int ans = 0;
    ans += bv[d].rank_lower(r, tb) - bv[d].rank_lower(l, tb);
    ans -= bv[d].rank_lower(r, sb) - bv[d].rank_lower(l, sb);
    // [sb << d, s)
    if((sb << shift) < s){
      int Lsb = bv[d].rank_lower(_L, sb), Lsb2 = bv[d].rank(_L, sb), Rsb = bv[d].rank_lower(_R, sb);
      int Lnext = _L + Rsb - Lsb, Rnext = Lnext + bv[d].rank(_R, sb) - Lsb2;
      int lnext = Lnext + bv[d].rank(l, sb) - Lsb2, rnext = Lnext + bv[d].rank(r, sb) - Lsb2;
      ans -= __range_freq(d - 1, lnext, rnext, Lnext, Rnext, S + (sb << shift), s, S + (sb << shift), S + ((sb + 1) << shift));
    }
    // [tb << d, t)
    if((tb << shift) < t){
      int Ltb = bv[d].rank_lower(_L, tb), Ltb2 = bv[d].rank(_L, tb), Rtb = bv[d].rank_lower(_R, tb);
      int Lnext = _L + Rtb - Ltb, Rnext = Lnext + bv[d].rank(_R, tb) - Ltb2;
      int lnext = Lnext + bv[d].rank(l, tb) - Ltb2, rnext = Lnext + bv[d].rank(r, tb) - Ltb2;
      ans += __range_freq(d - 1, lnext, rnext, Lnext, Rnext, S + (tb << shift), t, S + (tb << shift), S + ((tb + 1) << shift));
    }
    return ans;
  }
  std::pair<int, int> __range_kth_smallest(int l, int r, int k){
    int _L = 0, _R = n, ans = 0;
    for(int d = h - 1, shift = d * Rdiv; d >= 0; d--, shift -= Rdiv){
      int _l = -1, _r = R - 1;
      while(_r - _l > 1){
        int mid = (_l + _r) >> 1;
        (bv[d].rank_lower(r, mid + 1) - bv[d].rank_lower(l, mid + 1) > k ? _r : _l) = mid;
      }
      int dir = _r;
      k -= bv[d].rank_lower(r, dir) - bv[d].rank_lower(l, dir);
      ans += dir << shift;
      int Ldir = bv[d].rank_lower(_L, dir), Ldir2 = bv[d].rank(_L, dir), Rdir = bv[d].rank_lower(_R, dir);
      _L += Rdir - Ldir, _R = _L + bv[d].rank(_R, dir) - Ldir2;
      l = _L + bv[d].rank(l, dir) - Ldir2, r = _L + bv[d].rank(r, dir) - Ldir2;
    }
    return {ans, bottom_idx[l + k]};
  }
public:
  wavelet_matrix_arbitrary_radix(){}
  wavelet_matrix_arbitrary_radix(const std::vector<int> &v, int inf): n(v.size()), inf(inf){
    // 高速化のためRは2べき
    static_assert((1 << Rdiv) == R);
    h = 0;
    long long k = 1;
    while(k < inf) h++, k <<= Rdiv;
    build(v);
  }
  // r未満のcの数
  int rank(int r, int c){
    assert(0 <= r && r <= n);
    if(c < 0 || c >= inf) return 0;
    return __rank(r, c);
  }
  // k番目のcのインデックス, ない場合は-1
  int select(int k, int c){
    assert(0 <= k);
    if(c < 0 || c >= inf) return -1;
    return __select(k, c);
  }
  // [l, r)の[s, t)の数
  int range_freq(int l, int r, int s, int t){
    assert(0 <= l && r <= n);
    assert(0 <= s && t <= inf);
    if(l >= r || s >= t) return 0;
    return __range_freq(h - 1, l, r, 0, n, s, t, 0, 1 << (h * Rdiv));
  }
  // 区間[l, r)でk番目に小さい要素とそのインデックス, ない場合は-1
  // 値が同じ場合インデックスの昇順
  std::pair<int, int> range_kth_smallest(int l, int r, int k){
    if(r - l <= k) return {-1, -1};
    return __range_kth_smallest(l, r, k);
  }
  // 区間[l, r)でk番目に大きい要素とそのインデックス, ない場合は-1
  // 値が同じ場合インデックスの降順
  std::pair<int, int> range_kth_largest(int l, int r, int k){
    if(r - l <= k) return {-1, -1};
    return __range_kth_smallest(l, r, r - l - 1 - k);
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
    return l;
  }
};

template<int R, typename Val>
struct compressed_wavelet_matrix_arbitrary_radix{
private:
  int N;
  wavelet_matrix_arbitrary_radix<R> wm;
  std::vector<Val> rev;
  int lb(Val c){return std::lower_bound(rev.begin(), rev.end(), c) - rev.begin();}
public:
  compressed_wavelet_matrix_arbitrary_radix(){}
  compressed_wavelet_matrix_arbitrary_radix(const std::vector<Val> &_v): N(_v.size()){
    std::vector<int> v(N);
    std::vector<std::pair<Val, int>> tmp(N);
    for(int i = 0; i < N; i++) tmp[i] = {_v[i], i};
    std::sort(tmp.begin(), tmp.end());
    for(int i = 0; i < N; i++){
      if(i == 0 || tmp[i].first != tmp[i - 1].first) rev.push_back(tmp[i].first);
      v[tmp[i].second] = (int)rev.size() - 1;
    }
    tmp.clear();
    wm = wavelet_matrix_arbitrary_radix<R>(v, rev.size());
  }
  // val[k]
  Val access(int k){
    return rev[wm.access(k)];
  }
  // [0, r)のc
  int rank(int r, Val c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return 0;
    return wm.rank(r, idx);
  }
  // k番目のcのインデックス, ない場合は-1
  int select(int k, Val c){
    int idx = lb(c);
    if(idx == rev.size() || rev[idx] != c) return -1;
    return wm.select(k, idx);
  }
  // [l, r)の[s, t)の数
  int range_freq(int l, int r, Val s, Val t){
    return wm.range_freq(l, r, lb(s), lb(t));
  }
  // 区間[l, r)でk番目に小さい要素とそのインデックス, ない場合は-1
  // 値が同じ場合インデックスの昇順
  std::pair<Val, int> range_kth_smallest(int l, int r, int k){
    auto p = wm.range_kth_smallest(l, r, k);
    if(p.first == -1) return {-1, -1};
    return {rev[p.first], p.second};
  }
  // 区間[l, r)でk番目に大きい要素とそのインデックス, ない場合は-1
  // 値が同じ場合インデックスの降順
  std::pair<Val, int> range_kth_largest(int l, int r, int k){
    auto p = wm.range_kth_largest(l, r, k);
    if(p.first == -1) return {-1, -1};
    return {rev[p.first], p.second};
  }
  // [l, r)で値が[s, t)の要素でインデックス順にk番目, ない場合は-1
  int range_select(int l, int r, Val s, Val t, int k){
    return wm.range_select(l, r, lb(s), lb(t), k);
  }
};
#endif