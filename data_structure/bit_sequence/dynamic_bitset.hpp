#ifndef _DYNAMIC_BITSET_H_
#define _DYNAMIC_BITSET_H_
#include "bit_operation.hpp"
#include "../segment_tree/binary_indexed_tree.hpp"
#include <iostream>
#include <string>

// 同じサイズのbitset同士でしか演算できない
// サイズを変更する場合resize
// n, m, im: 最大サイズ, 最大ブロック, 内部的なサイズ(これ以降に1がないことが保証されている)
struct dynamic_bitset{
  static constexpr int BITLEN = 64, REM = 63, SHIFT = 6;
  int n, m, im;
  std::vector<uint64_t> v;

  uint64_t rightmost_mask(){
    return !(n & REM) ? ~0ULL : (1ULL << (n & REM)) - 1;
  }
  dynamic_bitset(): n(0), m(0), im(0){}
  dynamic_bitset(const dynamic_bitset &r): n(r.n), m(r.m), im(r.im), v(r.v.begin(), r.v.begin() + r.im){}
  dynamic_bitset(int n, bool f = 0): n(n), m((n + REM) >> SHIFT){
    if(f){
      im = m;
      v.resize(m, ~0ULL);
      v.back() &= rightmost_mask();
    }else im = 0;
  }
  // stringの左を0bit目とする
  dynamic_bitset(const std::string &s, char one = '1'): n(s.size()), m((n + REM) >> SHIFT), im(0){
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
      v[block] |= (1ULL << (k & REM));
    }else if(block < im) v[block] &= ~(1ULL << (k & REM));
  }
  bool get(int k){
    return (k >> SHIFT) < im ? (v[k >> SHIFT] >> (k & REM)) & 1 : 0;
  }
  bool any(){
    for(int i = 0; i < im; i++) if(v[i]) return true;
    return false;
  }
  int find_first(){
    for(int i = 0; i < im; i++) if(v[i]) return i * BITLEN + find_next_64bit(v[i], 0);
    return -1;
  }
  bool all(){
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
    m = (s + REM) >> SHIFT;
    if(m < im) v.resize(m), im = m;
    n = s;
    if(m == im) v[m - 1] &= rightmost_mask();
  }
  */
  dynamic_bitset operator ~(){
    dynamic_bitset res = *this;
    res.v.resize(m, ~0ULL);
    for(int i = 0; i < res.im; i++) res.v[i] = ~res.v[i];
    res.im = res.m;
    res.v[m - 1] &= res.rightmost_mask();
    return res;
  }
  dynamic_bitset operator ^ (const dynamic_bitset &r){
    dynamic_bitset res(*this);
    return res ^= r;
  }
  dynamic_bitset operator | (const dynamic_bitset &r){
    dynamic_bitset res(*this);
    return res |= r;
  }
  dynamic_bitset operator & (const dynamic_bitset &r){
    dynamic_bitset res(*this);
    return res &= r;
  }
  dynamic_bitset operator ^= (const dynamic_bitset &r){
    assert(n == r.n);
    if(im < r.im){
      im = r.im;
      v.resize(r.im, 0);
    }
    for(int i = 0; i < r.im; i++) v[i] ^= r.v[i];
    return *this;
  }
  dynamic_bitset operator |= (const dynamic_bitset &r){
    assert(n == r.n);
    if(im < r.im){
      im = r.im;
      v.resize(r.im, 0);
    }
    for(int i = 0; i < r.im; i++) v[i] |= r.v[i];
    return *this;
  }
  dynamic_bitset operator &= (const dynamic_bitset &r){
    assert(n == r.n);
    if(im > r.im){
      im = r.im;
      v.resize(im);
    }
    for(int i = 0; i < im; i++) v[i] &= r.v[i];
    return *this;
  }
  // 全てのビットが等しく, かつサイズも同じか
  bool operator == (const dynamic_bitset &r){
    if(n != r.n) return false;
    for(int i = 0; i < std::min(im, r.im); i++) if(v[i] != r.v[i]) return false;
    if(im > r.im){
      for(int i = r.im; i < im; i++) if(v[i]) return false;
    }else{
      for(int i = im; i < r.im; i++) if(r.v[i]) return false;
    }
    return true;
  }
  bool operator != (const dynamic_bitset &r){
    return !(*this == r);
  }
  // l |= r << s, l = rでもok
  static void lshift_or(dynamic_bitset &l, dynamic_bitset &r, int s){
    assert(l.n == r.n);
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
  static void lshift_xor(dynamic_bitset &l, dynamic_bitset &r, int s){
    assert(l.n == r.n);
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
  static void lshift_and(dynamic_bitset &l, dynamic_bitset &r, int s){
    assert(l.n == r.n);
    int t = s & REM, k = s >> SHIFT;
    if(!t){
      int imr = std::min(l.im, r.im + k);
      for(int i = imr - 1; i >= k; i--) l.v[i] &= r.v[i - k];
      for(int i = l.im - 1; i >= imr; i--) l.v[i] = 0;
      l.im = imr;
      return;
    }
    int imr = std::min(l.im, r.im + k + 1);
    int upper_bit = 64 - t;
    for(int i = imr - 1; i >= k; i--){
      if(i == k) l.v[i] &= r.v[i - k] << t;
      else l.v[i] &= (r.v[i - k] << t) | (r.v[i - k - 1] >> upper_bit);
    }
    for(int i = l.im - 1; i >= imr; i--) l.v[i] = 0;
    l.im = imr;
  }
  dynamic_bitset operator <<= (const int s){
    assert(0 <= s);
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
  dynamic_bitset operator >>= (const int s){
    assert(0 <= s);
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
  dynamic_bitset operator << (const int s){
    dynamic_bitset res(*this);
    return res <<= s;
  }
  dynamic_bitset operator >> (const int s){
    dynamic_bitset res(*this);
    return res >>= s;
  }
  int pop_count(){
    int res = 0;
    for(int i = 0; i < im; i++) res += __builtin_popcountll(v[i]);
    return res;
  }
};
std::ostream &operator<<(std::ostream &dest, const dynamic_bitset &v){
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