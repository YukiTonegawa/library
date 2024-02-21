#ifndef _BITSET_POWER2_H_
#define _BITSET_POWER2_H_
#include <vector>
#include <array>
#include <iostream>
#include <cassert>

std::vector<std::array<uint16_t, 1 << 16>> index_xor_table;
std::vector<std::array<uint16_t, 1 << 16>> index_or_table;
std::vector<std::array<uint16_t, 1 << 16>> index_and_table;
void index_xor_init(){
  if(!index_xor_table.empty()) return;
  index_xor_table.resize(16);
  for(int i = 0; i < (1 << 16); i++){
    for(int j = 0; j < 16; j++){
      int tmp = 0;
      for(int k = 0; k < 16; k++){
        if((i >> k) & 1){
          tmp += 1 << (k ^ j);
        }
      }
      index_xor_table[j][i] = tmp;
    }
  }
}
void index_or_init(){
  if(!index_or_table.empty()) return;
  index_or_table.resize(16);
  for(int i = 0; i < (1 << 16); i++){
    for(int j = 0; j < 16; j++){
      int tmp = 0;
      for(int k = 0; k < 16; k++){
        if((i >> k) & 1){
          tmp |= 1 << (k | j);
          //tmp |= 1 << (k ^ j);
        }
      }
      index_or_table[j][i] = tmp;
    }
  }
}
void index_and_init(){
  if(!index_and_table.empty()) return;
  index_and_table.resize(16);
  for(int i = 0; i < (1 << 16); i++){
    for(int j = 0; j < 16; j++){
      int tmp = (1 << 16) - 1;
      for(int k = 0; k < 16; k++){
        if((i >> k) & 1){
          tmp &= 1 << (k & j);
          //tmp |= 1 << (k ^ j);
        }
      }
      index_and_table[j][i] = tmp;
    }
  }
}
// resのk_bit目 = aの(k ^ x)_bit目
unsigned long long index_xor64(unsigned long long a, int x){
  static constexpr int mask16 = 0xFFFF;
  int A = a & mask16;
  int B = (a >> 16) & mask16;
  int C = (a >> 32) & mask16;
  int D = (a >> 48) & mask16;
  if(x >= 32){
    std::swap(A, C);
    std::swap(B, D);
    x -= 32;
  }
  if(x >= 16){
    std::swap(A, B);
    std::swap(C, D);
    x -= 16;
  }
  return index_xor_table[x][A] + (index_xor_table[x][B] << 16) + ((unsigned long long)index_xor_table[x][C] << 32) + ((unsigned long long)index_xor_table[x][D] << 48);
}

// resの(k | x)_bit目 = aのk_bit目
// (k | x)が存在しないものは0
unsigned long long index_or64(unsigned long long a, int x){
  static constexpr int mask16 = 0xFFFF;
  int A = a & mask16;
  int B = (a >> 16) & mask16;
  int C = (a >> 32) & mask16;
  int D = (a >> 48) & mask16;
  if(x >= 32){
    D |= B;
    B = 0;
    C |= A;
    A = 0;
    x -= 32;
  }
  if(x >= 16){
    B |= A;
    A = 0;
    D |= C;
    C = 0;
    x -= 16;
  }
  return index_or_table[x][A] + (index_or_table[x][B] << 16) + ((unsigned long long)index_or_table[x][C] << 32) + ((unsigned long long)index_or_table[x][D] << 48);
}
// resの(k & x)_bit目 = aのk_bit目
// (k & x)がないものは1
unsigned long long index_and64(unsigned long long a, int x){
  static constexpr int mask16 = 0xFFFF;
  int A = a & mask16;
  int B = (a >> 16) & mask16;
  int C = (a >> 32) & mask16;
  int D = (a >> 48) & mask16;
  if(x >= 32){
    x -= 32;
  }else{
    A &= C;
    C = mask16;
    B &= D;
    D = mask16;
  }
  if(x >= 16){
    x -= 16;
  }else{
    A &= B;
    B = mask16;
    C &= D;
    D = mask16;
  }
  return index_and_table[x][A] + (index_and_table[x][B] << 16) + ((unsigned long long)index_and_table[x][C] << 32) + ((unsigned long long)index_and_table[x][D] << 48);
}
// サイズ2^k
template<int k>
struct bitset_power2{
  static constexpr int w = 64; // ワードサイズ
  static constexpr int wmod = 63; // mod用
  static constexpr int wdiv = 6; // 除算用
  static constexpr int n = 1 << k; // bit数
  static constexpr int m = std::max(1, n >> wdiv); // 配列サイズ　
  using ull = unsigned long long;
  using Array = std::array<ull, m>;
  using Bit = bitset_power2<k>;
private:
  Array v;
public:
  bitset_power2(bool f = 0){
    static_assert(0 <= k && k < 30);
    index_xor_init();
    if(f) v.fill(~(ull)0);
    else v.fill(0);
  }
  Bit operator ^ (const Bit &r)const{Bit res(*this); return res ^= r;}
  Bit operator | (const Bit &r)const{Bit res(*this); return res |= r;}
  Bit operator & (const Bit &r)const{Bit res(*this); return res &= r;}
  Bit operator ^= (const Bit &r){for(int i = 0; i < m; i++){v[i] ^= r.v[i];}return *this;}
  Bit operator |= (const Bit &r){for(int i = 0; i < m; i++){v[i] |= r.v[i];}return *this;}
  Bit operator &= (const Bit &r){for(int i = 0; i < m; i++){v[i] &= r.v[i];}return *this;}
  bool operator == (const Bit &r)const{for(int i = 0; i < m; i++){if(v[i] != r.v[i]) return false;}return true;}
  bool operator != (const Bit &r)const{return !(*this == r);}
  bool any(){for(int i = 0; i < m; i++){if(v[i]) return true;}return false;}
  bool none(){return !any();}
  // 全て0にする
  void reset(){
    v.fill(0);
  }
  // i番目のbitをfにする
  void set(int i, bool f){
    f = f ^ ((v[i >> wdiv] >> (i & wmod)) & 1);
    v[i >> wdiv] ^= (ull)f <<  (i & wmod);
  }
  // i番目のbit
  bool get(int i)const{
    return (v[i >> wdiv] >> (i & wmod)) & 1;
  }
  // v′[i ^ x] = v[i]
  void index_xor(int x){
    assert(0 <= x && x < ((ull)1 << k));
    int xlarge = x >> wdiv, xsmall = x & wmod;
    if(xlarge){
      for(int i = 0; i < m; i++){
        int j = i ^ xlarge;
        if(i <= j){
          std::swap(v[i], v[j]);
          v[i] = index_xor64(v[i], xsmall);
          v[j] = index_xor64(v[j], xsmall);
        }
      }
    }else{
      for(int i = 0; i < m; i++) v[i] = index_xor64(v[i], xsmall);
    }
  }
  // v′[i | x] = v[i]
  // (i | x)が等しくなる場合それらのv[i]のor(他の演算にも変更できる)
  // (i | x)が存在しないものは全て0
  void index_or(int x){
    assert(0 <= x && x < ((ull)1 << k));
    int xlarge = x >> wdiv, xsmall = x & wmod;
    if(xlarge){
      for(int i = m - 1; i >= 0; i--){
        ull tmp = v[i];
        v[i] = 0;
        v[i | xlarge] |= tmp;
        // v[i | xlarge] ^= tmp; xorにしたい場合index_or64も変える　
      }
    }
    for(int i = 0; i < m; i++) v[i] = index_or64(v[i], xsmall);
  }
  // v′[i & x] = v[i]
  // (i & x)が等しくなる場合それらのv[i]のand(他の演算にも変更できる)
  // (i & x)が存在しないものは全て1
  void index_and(int x){
    assert(0 <= x && x < ((ull)1 << k));
    int xlarge = x >> wdiv, xsmall = x & wmod;
    if(xlarge){
      for(int i = 0; i < m; i++){
        ull tmp = v[i];
        v[i] = ~(ull)0;
        v[i & xlarge] &= tmp;
        // v[i & xlarge] ^= tmp; xorにしたい場合index_and64も変える　
      }
    }
    for(int i = 0; i < m; i++) v[i] = index_and64(v[i], xsmall);
  }
  Bit operator <<= (const int s){
    assert(0 <= s);
    int a = s & wmod, b = s >> wdiv;
    b = std::min(b, m);
    if(!a){
      for(int i = m - 1; i >= 0; i--) v[i] = (i >= b ? v[i - b] : 0);
      return *this;
    }
    int upper_bit = w - a;
    for(int i = m - 1; i > b; i--) v[i] = (v[i - b] << a) | (v[i - b - 1] >> upper_bit);
    if(b < m) v[b] = v[0] << a;
    for(int i = b - 1; i >= 0; i--) v[i] = 0;
    return *this;
  }
  Bit operator >>= (const int s){
    assert(0 <= s);
    int a = s & wmod, b = s >> wdiv;
    int m2 = m - b;
    if(!a){
      for(int i = 0; i < m2; i++) v[i] = v[i + b];
      for(int i = std::max(0, m2); i < m; i++) v[i] = 0;
      return *this;
    }
    int upper_bit = w - a;
    for(int i = 0; i < m2 - 1; i++) v[i] = (v[i + b] >> a) | (v[i + b + 1] << upper_bit);
    if(m2 - 1 >= 0) v[m2 - 1] = v[m - 1] >> a;
    for(int i = std::max(0, m2); i < m; i++) v[i] = 0;
    return *this;
  }
  Bit operator << (const int s){
    Bit res(*this);
    return res <<= s;
  }
  Bit operator >> (const int s){
    Bit res(*this);
    return res >>= s;
  }
  int popcount(){
    int res = 0;
    for(int i = 0; i < m; i++) res += __builtin_popcountll(v[i]);
    return res;
  }
};
template<int k>
std::ostream &operator<<(std::ostream &dest, const bitset_power2<k> &v){
  for(int i = 0; i < v.n; i++) dest << (int)v.get(i);
  return dest;
}
#endif