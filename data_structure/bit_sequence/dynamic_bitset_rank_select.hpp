#ifndef _DYNAMIC_BITSET_RANK_SELECT_H_
#define _DYNAMIC_BITSET_RANK_SELECT_H_
#include "bit_operation.hpp"
#include "../segment_tree/binary_indexed_tree.hpp"
#include <iostream>
#include <string>

// 同じサイズのbitset同士でしか演算できない
// サイズを変更する場合resize
// n, m, im: 最大サイズ, 最大ブロック, 内部的なサイズ(これ以降に1がないことが保証されている)
// rank, selectを実行する前にブロックごとのbinary_indexed_treeを再計算
struct dynamic_bitset_rank_select{
  static constexpr int BITLEN = 64, REM = 63, SHIFT = 6;
  int n, m, im;
  std::vector<uint64_t> v;
  binary_indexed_tree<int> bit_count;

  uint64_t rightmost_mask(){
    return !(n & REM) ? ~0ULL : (1ULL << (n & REM)) - 1;
  }
  dynamic_bitset_rank_select(): n(0), m(0), im(0){}
  dynamic_bitset_rank_select(const dynamic_bitset_rank_select &r): n(r.n), m(r.m),
  im(r.im), v(r.v.begin(), r.v.begin() + r.im), bit_updated(false){}
  dynamic_bitset_rank_select(int n, bool f = 0): n(n), m((n + REM) >> SHIFT), bit_updated(false){
    if(f){
      im = m;
      v.resize(m, ~0ULL);
      v.back() &= rightmost_mask();
    }else im = 0;
  }
  // stringの左を0bit目とする
  dynamic_bitset_rank_select(const std::string &s, char one = '1'): n(s.size()), m((n + REM) >> SHIFT), im(0), bit_updated(false){
    uint64_t sum = 0;
    for(int i = 0, j = 0, t = 0; i < n; i++){
      sum += (uint64_t)(s[i] == one) << j;
      if(i == n - 1 || j == REM){
        t++;
        v.push_back(sum);
        if(sum) im = t;
        sum = j = 0;
      }else j++;
    }
  }
  int size(){
    return n;
  }
  void set(int k, bool f){
    assert(0 <= k && k < n);
    int block = k >> SHIFT;
    if(f){
      if(block >= im) im = block + 1;
      while(block >= v.size()) v.resize(std::min((int)v.size() * 2 + 1, m), 0);
      if(bit_updated && !((v[block] >> (k & REM)) & 1)) bit_count.update(block, 1);
      v[block] |= (1ULL << (k & REM));
    }else if(block < im){
      if(bit_updated && ((v[block] >> (k & REM)) & 1)) bit_count.update(block, -1);
      v[block] &= ~(1ULL << (k & REM));
    }
  }
  bool get(int k){
    return (k >> SHIFT) < im ? (v[k >> SHIFT] >> (k & REM)) & 1 : 0;
  }
  bool any(){
    if(bit_updated) return bit_count.query(m);
    for(int i = 0; i < im; i++) if(v[i]) return true;
    return false;
  }
  bool all(){
    if(bit_updated) return bit_count.query(m) == n;
    if(im != m) return false;
    for(int i = 0; i < m - 1; i++) if(v[i] != ~0ULL) return false;
    if(v[m - 1] != rightmost_mask()) return false;
    return true;
  }
  bool none(){
    return !any();
  }
  /*
  void resize(int s){
    bit_updated = false;
    m = (s + REM) >> SHIFT;
    if(m < im) v.resize(m), im = m;
    n = s;
    if(m == im) v[m - 1] &= rightmost_mask();
  }
  */
  dynamic_bitset_rank_select operator ~(){
    dynamic_bitset_rank_select res = *this;
    res.v.resize(m, ~0ULL);
    for(int i = 0; i < res.im; i++) res.v[i] = ~res.v[i];
    res.im = res.m;
    res.v[m - 1] &= res.rightmost_mask();
    return res;
  }
  dynamic_bitset_rank_select operator ^ (const dynamic_bitset_rank_select &r){
    dynamic_bitset_rank_select res(*this);
    return res ^= r;
  }
  dynamic_bitset_rank_select operator | (const dynamic_bitset_rank_select &r){
    dynamic_bitset_rank_select res(*this);
    return res |= r;
  }
  dynamic_bitset_rank_select operator & (const dynamic_bitset_rank_select &r){
    dynamic_bitset_rank_select res(*this);
    return res &= r;
  }
  dynamic_bitset_rank_select operator ^= (const dynamic_bitset_rank_select &r){
    assert(n == r.n);
    bit_updated = false;
    if(im < r.im){
      im = r.im;
      v.resize(r.im, 0);
    }
    for(int i = 0; i < r.im; i++) v[i] ^= r.v[i];
    return *this;
  }
  dynamic_bitset_rank_select operator |= (const dynamic_bitset_rank_select &r){
    assert(n == r.n);
    bit_updated = false;
    if(im < r.im){
      im = r.im;
      v.resize(r.im, 0);
    }
    for(int i = 0; i < r.im; i++) v[i] |= r.v[i];
    return *this;
  }
  dynamic_bitset_rank_select operator &= (const dynamic_bitset_rank_select &r){
    assert(n == r.n);
    bit_updated = false;
    if(im > r.im){
      im = r.im;
      v.resize(im);
    }
    for(int i = 0; i < im; i++) v[i] &= r.v[i];
    return *this;
  }
  // 全てのビットが等しく, かつサイズも同じか
  bool operator == (const dynamic_bitset_rank_select &r){
    if(n != r.n) return false;
    for(int i = 0; i < std::min(im, r.im); i++) if(v[i] != r.v[i]) return false;
    if(im > r.im){
      for(int i = r.im; i < im; i++) if(v[i]) return false;
    }else{
      for(int i = im; i < r.im; i++) if(r.v[i]) return false;
    }
    return true;
  }
  bool operator != (const dynamic_bitset_rank_select &r){
    return !(*this == r);
  }
  // l |= r << s, l = rでもok
  static void lshift_or(dynamic_bitset_rank_select &l, dynamic_bitset_rank_select &r, int s){
    assert(l.n == r.n);
    l.bit_updated = false;
    int t = s & REM, k = s >> SHIFT;
    if(!t){
      int imr = std::min(r.m, r.im + k);
      l.im = std::max(l.im, imr);
      if(imr > l.v.size()) l.v.resize(imr, 0);
      for(int i = imr - 1; i >= k; i--) l.v[i] |= r.v[i - k];
      if(imr == l.m) l.v.back() &= l.rightmost_mask();
      return;
    }
    int imr = std::min(r.m, r.im + k + 1);
    l.im = std::max(l.im, imr);
    if(imr > l.v.size()) l.v.resize(imr, 0);
    int upper_bit = 64 - t;
    for(int i = imr - 1; i >= k; i--){
      if(i == k) l.v[i] |= r.v[i - k] << t;
      else l.v[i] |= (r.v[i - k] << t) | (r.v[i - k - 1] >> upper_bit);
    }
    if(imr == l.m) l.v.back() &= l.rightmost_mask();
  }
  // l ^= r << s, l = rでもok
  static void lshift_xor(dynamic_bitset_rank_select &l, dynamic_bitset_rank_select &r, int s){
    assert(l.n == r.n);
    l.bit_updated = false;
    int t = s & REM, k = s >> SHIFT;
    if(!t){
      int imr = std::min(r.m, r.im + k);
      l.im = std::max(l.im, imr);
      if(imr > l.v.size()) l.v.resize(imr, 0);
      for(int i = imr - 1; i >= k; i--) l.v[i] ^= r.v[i - k];
      if(imr == l.m) l.v.back() &= l.rightmost_mask();
      return;
    }
    int imr = std::min(r.m, r.im + k + 1);
    l.im = std::max(l.im, imr);
    if(imr > l.v.size()) l.v.resize(imr, 0);
    int upper_bit = 64 - t;
    for(int i = imr - 1; i >= k; i--){
      if(i == k) l.v[i] ^= r.v[i - k] << t;
      else l.v[i] ^= (r.v[i - k] << t) | (r.v[i - k - 1] >> upper_bit);
    }
    if(imr == l.m) l.v.back() &= l.rightmost_mask();
  }
  // l &= r << s, l = rでもok
  static void lshift_and(dynamic_bitset_rank_select &l, dynamic_bitset_rank_select &r, int s){
    assert(l.n == r.n);
    l.bit_updated = false;
    int t = s & REM, k = s >> SHIFT;
    if(!t){
      int imr = std::min(l.im, r.im + k);
      for(int i = imr - 1; i >= 0; i--) l.v[i] &= (i >= k ? r.v[i - k] : 0);
      for(int i = l.im - 1; i >= imr; i--) l.v[i] = 0;
      l.im = imr;
      return;
    }
    int imr = std::min(l.im, r.im + k + 1);
    int upper_bit = 64 - t;
    for(int i = imr - 1; i >= 0; i--){
      if(i < k) l.v[i] = 0;
      else if(i == k) l.v[i] &= r.v[i - k] << t;
      else l.v[i] &= (r.v[i - k] << t) | (r.v[i - k - 1] >> upper_bit);
    }
    for(int i = l.im - 1; i >= imr; i--) l.v[i] = 0;
    l.im = imr;
  }
  dynamic_bitset_rank_select operator <<= (const int s){
    assert(0 <= s);
    bit_updated = false;
    int t = s & REM, k = s >> SHIFT;
    if(!t){
      im = std::min(m, im + k);
      if(im > v.size()) v.resize(im, 0);
      for(int i = im - 1; i >= 0; i--) v[i] = (i >= k ? v[i - k] : 0);
      if(im == m) v.back() &= rightmost_mask();
      return *this;
    }
    im = std::min(m, im + k + 1);
    if(im > v.size()) v.resize(im, 0);
    int upper_bit = 64 - t;
    for(int i = im - 1; i >= 0; i--){
      if(i < k) v[i] = 0;
      else if(i == k) v[i] = v[i - k] << t;
      else v[i] = (v[i - k] << t) | (v[i - k - 1] >> upper_bit);
    }
    if(im == m) v.back() &= rightmost_mask();
    return *this;
  }
  dynamic_bitset_rank_select operator >>= (const int s){
    assert(0 <= s);
    bit_updated = false;
    int t = s & REM, k = s >> SHIFT;
    im -= k;
    if(!t){
      for(int i = 0; i < im; i++) v[i] = v[i + k];
      for(int i = std::max(0, im); i < im + k; i++) v[i] = 0;
      if(im < 0) im = 0;
      return *this;
    }
    int upper_bit = 64 - t;
    for(int i = 0; i < im; i++){
      if(i + k == (int)v.size() - 1) v[i] = v[i + k] >> t;
      else v[i] = (v[i + k] >> t) | (v[i + k + 1] << upper_bit);
    }
    for(int i = std::max(0, im); i < im + k; i++) v[i] = 0;
    if(im < 0) im = 0;
    return *this;
  }
  dynamic_bitset_rank_select operator << (const int s){
    dynamic_bitset_rank_select res(*this);
    return res <<= s;
  }
  dynamic_bitset_rank_select operator >> (const int s){
    dynamic_bitset_rank_select res(*this);
    return res >>= s;
  }
  bool bit_updated;
  void bit_recalc(){
    if(bit_updated) return;
    if(bit_count.sum.size() != m + 1) bit_count = binary_indexed_tree<int>(m);
    for(int i = 0; i < m; i++){
      if(i < im) bit_count.sum[i + 1] = __builtin_popcountll(v[i]);
      else bit_count.sum[i + 1] = 0;
    }
    for(int i = 1; i <= m; i++){
      int nxt = i + (i & (-i));
      if(nxt <= m) bit_count.sum[nxt] += bit_count.sum[i];
    }
    bit_updated = true;
  }
  int pop_count(){
    if(!bit_updated) bit_recalc();
    return bit_count.query(m);
  }
  int rank1(int r){
    assert(0 <= r && r <= n);
    if(!bit_updated) bit_recalc();
    int block = r >> SHIFT;
    return bit_count.query(block) + (im <= block ? 0 : rank_64bit(v[block], r & REM));
  }
  int rank0(int r){
    assert(0 <= r && r <= n);
    return r - rank1(r);
  }
  // k個目(0-indexed)の1, 無い場合は-1
  int select1(int k){
    if(!bit_updated) bit_recalc();
    auto [b, c] = bit_count.lower_bound(k + 1);
    if(b == m) return -1;
    return (b << SHIFT) + select_64bit(v[b], k - c + __builtin_popcountll(v[b]));
  }
  // k個目(0-indexed)の0, 無い場合は-1
  int select0(int k){
    if(!bit_updated) bit_recalc();
    int x = k + 1, y = 1 << bit_count.H, h = bit_count.H;
    int s = 0, t = 0;
    while(h--){
      if(bit_count.M < y) y -= 1 << h;
      else{
        int sum_zero = s + (1 << (h + 1 + SHIFT)) - bit_count.sum[y]; // 区間長 = 2^(h + 1) * 64 = 1 << (h + 1 + SHIFT)
        if(x <= sum_zero) t = sum_zero, y -= 1 << h;
        else s = sum_zero, y += 1 << h;
      }
    }
    if(y == bit_count.M + 1) return -1;
    int b, c;
    if(x <= s + (1 << SHIFT) - bit_count.sum[y]) b = y - 1, c = s + (1 << SHIFT) - bit_count.sum[y];
    else b = y, c = t;
    int i = (b << SHIFT);
    if(b >= im) i += k - c + BITLEN;
    else i += select_64bit(~v[b], k - c + __builtin_popcountll(~v[b]));
    return i >= n ? -1 : i;
  }
  // i以降(i含む)の1, 無い場合は-1
  int find_next1(int i){
    assert(0 <= i && i < n);
    if((i >> SHIFT) >= im) return -1;
    int _i = i & REM, _j = find_next_64bit(v[i >> SHIFT], _i);
    if(_j != -1) return i - _i + _j;
    return select1(rank1(i));
  }
  // i以前(i含む)の1, 無い場合は-1
  int find_prev1(int i){
    assert(0 <= i && i < n);
    if((i >> SHIFT) >= im){
      if(!im) return -1;
      i = (im << SHIFT) - 1;
    }
    int _i = i & REM, _j = find_prev_64bit(v[i >> SHIFT], _i);
    if(_j != -1) return i - _i + _j;
    int r = rank1(i + 1);
    if(!r) return -1;
    return select1(r - 1);
  }
  // i以降(i含む)の0, 無い場合は-1
  int find_next0(int i){
    assert(0 <= i && i < n);
    if((i >> SHIFT) >= im) return i;
    int _i = i & REM, _j = find_next_64bit(~v[i >> SHIFT], _i);
    if(_j != -1){
      int res = i - _i + _j;
      return res >= n ? -1 : res;
    }
    return select0(rank0(i));
  }
  // i以前(i含む)の0, 無い場合は-1
  int find_prev0(int i){
    assert(0 <= i && i < n);
    if((i >> SHIFT) >= im) return i;
    int _i = i & REM, _j = find_prev_64bit(~v[i >> SHIFT], _i);
    if(_j != -1) return i - _i + _j;
    int r = rank0(i + 1);
    if(!r) return -1;
    return select0(r - 1);
  }
};
std::ostream &operator<<(std::ostream &dest, const dynamic_bitset_rank_select &v){
  if(!v.n) return dest;
  int cnt = 0;
  for(int i = 0; i < v.m; i++){
    if(i >= v.im){
      for(int j = 0; j < 64 && cnt < v.n; j++, cnt++) dest << 0;
    }else{
      for(int j = 0; j < 64 && cnt < v.n; j++, cnt++) dest << ((v.v[i] >> j) & 1);
    }
    if(i != v.m - 1) dest << ' ';
  }
  return dest;
}
#endif