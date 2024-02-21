#ifndef _ENUMERATE_H_
#define _ENUMERATE_H_
#include <vector>
#include <cassert>

// 自然数の分割の列挙
// 辞書順で小さい順
// {1, 1, 1, 1} -> {1, 1, 2} -> {1, 3} -> {2, 2} -> {4}
// n = 100の190569292通りを全列挙すると2secくらい
struct enumerate_partition{
  int n;
  std::vector<int> v;
  enumerate_partition(int n): n(n), v(n, 1){
    assert(0 < n);
  }
  // 次がある場合は更新してtrue, {n}の場合は何もせずfalse
  bool next(){
    if(v.size() == 1) return false;
    int sz = v.size();
    if(v[sz - 1] - v[sz - 2] >= 2){
      int x = v[sz - 1] - 1;
      v[sz - 2]++;
      int next_sz = sz - 1 + (x / v[sz - 2]);
      v.resize(next_sz);
      for(int j = sz - 1; j < next_sz; j++) v[j] = v[sz - 2];
      v.back() = v[sz - 2] + x % v[sz - 2];
    }else{
      int x = v.back();
      v.pop_back();
      v.back() += x;
    }
    return true;
  }
  // 配列 -> ハッシュ
  // 開始の1 + (n - 1個の隙間に仕切りを入れるか(0の時入れる)) {1, 2, 3} -> 1010110
  unsigned long long encode(const std::vector<int> &val)const{
    assert(n <= 63);
    unsigned long long res = 1;
    for(int x : val){
      res <<= x - 1;
      res |= (1ULL << (x - 1)) - 1;
      res <<= 1;
    }
    return res;
  }
  // ハッシュ -> 配列
  std::vector<int> decode(unsigned long long x)const{
    assert(n <= 63 && x && x % 2 == 0);
    std::vector<int> res;
    int cnt = 0;
    for(int i = 62 - __builtin_clzll(x); i >= 0; i--){
      if((x >> i) & 1){
        cnt++;
      }else{
        res.push_back(cnt + 1);
        cnt = 0;
      }
    }
    return res;
  }
};

// 自然数の分割の列挙(順序付き)
// 辞書順で小さい順
// 2 ^ (n - 1)通り
// {1, 1, 1, 1} -> {1, 1, 2} -> {1, 2, 1} -> {1, 3} -> {2, 1, 1} -> {2, 2} -> {3, 1} -> {4}
struct enumerate_ordered_partition{
  int n, itr = 0;
  enumerate_ordered_partition(int n): n(n){
    assert(0 < n && n < 30);
  }
  // 次がある場合は更新してtrue, {n}の場合は何もせずfalse
  bool next(){
    if(itr + 1 == (1ULL << (n - 1))) return false;
    itr++;
    return true;
  }
  // 前がある場合は更新してtrue, {1, 1, 1...}の場合は何もせずfalse
  bool prev(){
    if(itr == 0) return false;
    itr--;
    return true;
  }
  // 現在の配列を取得
  std::vector<int> get(){
    return decode((itr << 1) | (1ULL << n));
  }
  // 配列 -> ハッシュ
  // 開始の1 + (n - 1個の隙間に仕切りを入れるか(0の時入れる)) {1, 2, 3} -> 1010110
  unsigned long long encode(const std::vector<int> &val)const{
    assert(n <= 63);
    unsigned long long res = 1;
    for(int x : val){
      res <<= x - 1;
      res |= (1ULL << (x - 1)) - 1;
      res <<= 1;
    }
    return res;
  }
  // ハッシュ -> 配列
  std::vector<int> decode(unsigned long long x)const{
    assert(n <= 63 && x && x % 2 == 0);
    std::vector<int> res;
    int cnt = 0;
    for(int i = 62 - __builtin_clzll(x); i >= 0; i--){
      if((x >> i) & 1){
        cnt++;
      }else{
        res.push_back(cnt + 1);
        cnt = 0;
      }
    }
    return res;
  }
};

template<typename T>
struct enumerate_permutation{
  int n;
  std::vector<T> v;
  enumerate_permutation(const std::vector<T> &v): n(v.size()), v(v){
    assert(0 < n);
  }
  // 次の状態がある場合は更新してtrue, そうでない場合は何もせずfalse
  bool next(){
    return std::next_permutation(v.begin(), v.end());
  }
  // 前の状態ある場合は更新してtrue, そうでない場合は何もせずfalse
  bool prev(){
    return std::prev_permutation(v.begin(), v.end());
  }
  // 配列 -> ハッシュ
  // 開始の1, 要素の数の1, 仕切りの0 {1, 2, 3} -> 1101101110
  // assert (任意の値 >= 0) かつ (総和 + n <= 63)
  unsigned long long encode_small_sum(const std::vector<T> &val)const{
    unsigned long long res = 1;
    T S = 0;
    for(T x : val){
      assert(x >= 0);
      res <<= x;
      res |= (1ULL << x) - 1;
      res <<= 1;
      S += x;
      assert(S + n <= 63);
    }
    return res;
  }
  // ハッシュ -> 配列
  std::vector<T> decode_small_sum(unsigned long long x)const{
    assert(x && x % 2 == 0);
    std::vector<T> res;
    int cnt = 0;
    for(int i = 62 - __builtin_clzll(x); i >= 0; i--){
      if((x >> i) & 1){
        cnt++;
      }else{
        res.push_back(cnt);
        cnt = 0;
      }
    }
    return res;
  }
};



#endif