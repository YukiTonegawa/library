#ifndef _TEMPLATE_H_
#define _TEMPLATE_H_
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <stack>
#include <queue>
#include <deque>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <bitset>
#include <cmath>
#include <functional>
#include <cassert>
#include <climits>
#include <iomanip>
#include <numeric>
#include <memory>
#include <random>
#include <thread>
#include <chrono>
#define allof(obj) (obj).begin(), (obj).end()
#define range(i, l, r) for(int i=l;i<r;i++)
#define unique_elem(obj) obj.erase(std::unique(allof(obj)), obj.end())
#define bit_subset(i, S) for(int i=S, zero_cnt=0;(zero_cnt+=i==S)<2;i=(i-1)&S)
#define bit_kpop(i, n, k) for(int i=(1<<k)-1,x_bit,y_bit;i<(1<<n);x_bit=(i&-i),y_bit=i+x_bit,i=(!i?(1<<n):((i&~y_bit)/x_bit>>1)|y_bit))
#define bit_kth(i, k) ((i >> k)&1)
#define bit_highest(i) (i?63-__builtin_clzll(i):-1)
#define bit_lowest(i) (i?__builtin_ctzll(i):-1)
#define sleepms(t) std::this_thread::sleep_for(std::chrono::milliseconds(t))
using ll = long long;
using ld = long double;
using ul = uint64_t;
using pi = std::pair<int, int>;
using pl = std::pair<ll, ll>;
using namespace std;

template<typename F, typename S>
std::ostream &operator<<(std::ostream &dest, const std::pair<F, S> &p){
  dest << p.first << ' ' << p.second;
  return dest;
}
template<typename T>
std::ostream &operator<<(std::ostream &dest, const std::vector<std::vector<T>> &v){
  int sz = v.size();
  if(sz==0) return dest;
  for(int i=0;i<sz;i++){
    int m = v[i].size();
    for(int j=0;j<m;j++) dest << v[i][j] << (i!=sz-1&&j==m-1?'\n':' ');
  }
  return dest;
}
template<typename T>
std::ostream &operator<<(std::ostream &dest, const std::vector<T> &v){
  int sz = v.size();
  if(sz==0) return dest;
  for(int i=0;i<sz-1;i++) dest << v[i] << ' ';
  dest << v[sz-1];
  return dest;
}
template<typename T, size_t sz>
std::ostream &operator<<(std::ostream &dest, const std::array<T, sz> &v){
  if(sz==0) return dest;
  for(int i=0;i<sz-1;i++) dest << v[i] << ' ';
  dest << v[sz-1];
  return dest;
}
template<typename T>
std::ostream &operator<<(std::ostream &dest, const std::set<T> &v){
  for(auto itr=v.begin();itr!=v.end();){
    dest << *itr;
    itr++;
    if(itr!=v.end()) dest << ' ';
  }
  return dest;
}
template<typename T, typename E>
std::ostream &operator<<(std::ostream &dest, const std::map<T, E> &v){
  for(auto itr=v.begin();itr!=v.end();){
    dest << '(' << itr->first << ", " << itr->second << ')';
    itr++;
    if(itr!=v.end()) dest << '\n';
  }
  return dest;
}
std::ostream &operator<<(std::ostream &dest, __int128_t value) {
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[128];
    char *d = std::end(buffer);
    do {
      --d;
      *d = "0123456789"[tmp % 10];
      tmp /= 10;
    } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}
template<typename T>
vector<T> make_vec(size_t sz, T val){return std::vector<T>(sz, val);}
template<typename T, typename... Tail>
auto make_vec(size_t sz, Tail ...tail){
  return std::vector<decltype(make_vec<T>(tail...))>(sz, make_vec<T>(tail...));
}
template<typename T>
vector<T> read_vec(size_t sz){
  std::vector<T> v(sz);
  for(int i=0;i<(int)sz;i++) std::cin >> v[i];
  return v;
}
template<typename T, typename... Tail>
auto read_vec(size_t sz, Tail ...tail){
  auto v = std::vector<decltype(read_vec<T>(tail...))>(sz);
  for(int i=0;i<(int)sz;i++) v[i] = read_vec<T>(tail...);
  return v;
}
long long max(long long a, int b){return std::max(a, (long long)b);}
long long max(int a, long long b){return std::max((long long)a, b);}
long long min(long long a, int b){return std::min(a, (long long)b);}
long long min(int a, long long b){return std::min((long long)a, b);}
long long modulo(long long a, long long m){a %= m; return a < 0 ? a + m : a;}

void io_init(){
  std::cin.tie(nullptr);
  std::ios::sync_with_stdio(false);
}
#endif