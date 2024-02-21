#ifndef _RANDOMS_H_
#define _RANDOMS_H_
#include <vector>
#include <random>
#include <set>
#include "../traits.hpp"
// Tの[MIN, MAX]から一様ランダム, 2回目以降は同じ値が返される
template<class T, std::enable_if_t<is_intle64<T>::value>* = nullptr>
T random_once(){
  std::random_device seed_gen;
  std::mt19937_64 engine(seed_gen());
  static T res = (T)engine();
  return res;
}
// Tの[MIN, MAX]から一様ランダム
template<class T, std::enable_if_t<is_intle64<T>::value>* = nullptr>
T random_number(){
  std::random_device seed_gen;
  std::mt19937_64 engine(seed_gen());
  return (T)engine();
}
// [l, r]から一様ランダム
template<class T, std::enable_if_t<is_signed_intle64<T>::value>* = nullptr>
T random_number(T l, T r){
  std::random_device seed_gen;
  std::mt19937_64 engine(seed_gen());
  assert(l <= r);
  __int128_t len = (__int128_t)r - l + 1;
  return (T)((__int128_t)engine() % len + l);
}
// [l, r]から一様ランダム
template<class T, std::enable_if_t<is_unsigned_intle64<T>::value>* = nullptr>
T random_number(T l, T r){
  std::random_device seed_gen;
  std::mt19937_64 engine(seed_gen());
  assert(l <= r);
  if(l == std::numeric_limits<T>::min() && r == std::numeric_limits<T>::max()) return (T)engine();
  T len = (r - l + 1);
  return (T)engine() % len + l;
}
// 128bit, 小数は未実装

// n要素のvector(各要素は[l, r]から独立して一様ランダムに選ばれる)
template<class T, std::enable_if_t<std::is_integral<T>::value>* = nullptr>
std::vector<T> random_vector(int n, T l, T r){
  std::vector<T> res(n);
  for(int i = 0; i < n; i++) res[i] = random_number<T>(l, r);
  return res;
}

// [l, r]の各要素が1つずつあるくじからn回引く
template<class T, std::enable_if_t<std::is_integral<T>::value>* = nullptr>
std::vector<T> random_choice(int n, T l, T r){
  static constexpr T inf = std::numeric_limits<T>::max();
  static constexpr T minf = std::numeric_limits<T>::min();
  assert(n >= 0);
  bool sparse = false;
  if(l == minf && r == inf){
    sparse = true;
  }else{
    T len = (r - l + 1);
    assert(n <= len);
    sparse = (2 * n <= len);
  }
  std::vector<T> res;
  if(sparse){
    std::set<T> used;
    while(res.size() < n){
      T x = random_number<T>(l, r);
      if(used.find(x) == used.end()){
        res.push_back(x);
        used.insert(x);
      }
    }
  }else{
    std::vector<std::pair<unsigned long long, T>> v;
    for(T i = l; i <= r; i++){
      v.push_back({random_number<unsigned long long>(), i});
    }
    std::sort(v.begin(), v.end());
    for(int i = 0; i < n; i++) res.push_back(v[i].second);
  }
  return res;
}

// [l, r]のn要素の集合(ソート済)
template<class T, std::enable_if_t<std::is_integral<T>::value>* = nullptr>
std::vector<T> random_set(int n, T l, T r){
  auto res = random_choice<T>(n, l, r);
  std::sort(res.begin(), res.end());
  return res;
}

// [1, n]の順列
std::vector<int> random_permutation(int n){
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::vector<int> res(n);
  std::iota(res.begin(), res.end(), 1);
  std::shuffle(res.begin(), res.end(), engine);
  return res;
}
template<class T>
void random_shuffle(std::vector<T> &v){
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::shuffle(v.begin(), v.end(), engine);
}

// n文字のstring(各文字が[a, z]から一様ランダム)
std::string random_string(int n, char a = 'a', char z = 'z'){
  assert(a <= z);
  auto tmp = random_vector<int>(n, (int)a, (int)z);
  std::string res = "";
  for(int i = 0; i < n; i++) res += (char)tmp[i];
  return res;
}
// l <= a <= b <= rを満たすa, bを返す(閉区間)
template<class T, std::enable_if_t<is_signed_intle64<T>::value>* = nullptr>
std::pair<T, T> random_segment1(T l, T r){
  assert(l <= r);
  T res1 = random_number<T>(l, r);
  T res2 = random_number<T>(l, r);
  if(res1 > res2) std::swap(res1, res2);
  return {res1, res2};
}
// l <= a <= b <= rを満たすa, bを返す(閉区間)
// 距離を固定してからそれを満たすl, rを返す(random_segment1と比べて大きい区間が出来やすい)
template<class T, std::enable_if_t<is_signed_intle64<T>::value>* = nullptr>
std::pair<T, T> random_segment2(T l, T r){
  assert(l <= r);
  unsigned long long d = random_number<unsigned long long>(0, (__int128_t)r - l);
  T res_l = random_number<T>(l, (__int128_t)r - d);
  return {res_l, (__int128_t)res_l + d};
}

// [l, r]をk分割
template<class T, std::enable_if_t<is_signed_intle64<T>::value>* = nullptr>
std::vector<std::pair<T, T>> random_partition(int k, T l, T r){
  assert(k >= 1);
  if(k == 1) return {{l, r}};
  auto p = random_set<T>(k - 1, l, r - 1);
  p.push_back(r);
  T prel = l;
  std::vector<std::pair<T, T>> res;
  for(auto x : p){
    res.push_back({prel, x});
    prel = x + 1;
  }
  return res;
}
/*
unsigned long long random_prime(unsigned long long l, unsigned long long r){
  return generate_random_prime(l, r);
}
*/
namespace random_tree{
  // 直径の最小値を決める
  std::vector<std::pair<int, int>> random_tree_line(int n){
    int d = random_number<int>(1, n);
    std::vector<std::pair<int, int>> res;
    for(int i = 1; i + 1 <= d; i++) res.push_back({i, i + 1});
    for(int i = d + 1; i <= n; i++) res.push_back({i, random_number<int>(1, i - 1)});
    assert(res.size() == n - 1);
    auto p = random_permutation(n);
    for(int i = 0; i < n - 1; i++){
      res[i].first = p[res[i].first - 1];
      res[i].second = p[res[i].second - 1];
    }
    return res;
  }
  // ウニの中心の最小次数を決める
  std::vector<std::pair<int, int>> random_tree_star(int n){
    int d = random_number<int>(1, n);
    std::vector<std::pair<int, int>> res;
    for(int i = 2; i <= d; i++) res.push_back({1, i});
    for(int i = d + 1; i <= n; i++) res.push_back({i, random_number<int>(1, i - 1)});
    assert(res.size() == n - 1);
    auto p = random_permutation(n);
    for(int i = 0; i < n - 1; i++){
      res[i].first = p[res[i].first - 1];
      res[i].second = p[res[i].second - 1];
    }
    return res;
  }
  // 次数の最大値を決める
  std::vector<std::pair<int, int>> random_tree_binary(int n, int max_deg = 3){
    assert(max_deg >= 2);
    std::vector<int> deg(n + 1, 0);
    std::vector<std::pair<int, int>> res;
    for(int i = 2; i <= n; i++){
      int x = random_number<int>(1, i - 1);
      while(deg[x] >= max_deg) x = random_number<int>(1, i - 1);
      deg[x]++;
      deg[i]++;
      res.push_back({x, i});
    }
    assert(res.size() == n - 1);
    auto p = random_permutation(n);
    for(int i = 0; i < n - 1; i++){
      res[i].first = p[res[i].first - 1];
      res[i].second = p[res[i].second - 1];
    }
    return res;
  }
};
#endif