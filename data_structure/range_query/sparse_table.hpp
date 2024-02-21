#ifndef _SPARSE_TABLE_H_
#define _SPARSE_TABLE_H_
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <limits>

template<typename T, T (*merge)(T, T), T (*id)()>
struct sparse_table{
  int n;
  std::vector<std::vector<T>> table;
  sparse_table(): n(0){}
  sparse_table(const std::vector<T> &v): n(v.size()){
    int m = 0;
    while((1 << m) <= n) m++;
    table.resize(m, std::vector<T>(n));
    table[0] = v;
    for(int i = 1; i < m; i++){
      for(int j = 0; j + (1 << (i - 1)) < n; j++){
        table[i][j] = merge(table[i - 1][j], table[i - 1][j + (1 << (i - 1))]);
      }
    }
  }
  // 最左bit
  int msb(int x){
    return 31 - __builtin_clz(x);
  }
  T query(int l, int r){
    l = std::max(l, 0);
    r = std::min(r, n);
    if(l >= r) return id();
    int len = r - l;
    int b = msb(len);
    return merge(table[b][l], table[b][r - (1 << b)]);
  }
};

template<typename T, T (*merge)(T, T), T (*id)()>
struct sparse_table_short{
  int n;
  std::vector<std::vector<T>> table;
  sparse_table_short(): n(0){}
  sparse_table_short(const std::vector<T> &v, int lim): n(v.size()){
    int m = 0;
    while((1 << m) <= n) m++;
    m = std::min(m, lim);
    table.resize(m, std::vector<T>(n));
    table[0] = v;
    for(int i = 1; i < m; i++){
      for(int j = 0; j + (1 << (i - 1)) < n; j++){
        table[i][j] = merge(table[i - 1][j], table[i - 1][j + (1 << (i - 1))]);
      }
    }
  }
  // 最左bit
  int msb(int x){
    return (x == 1 ? 0 : 31 - __builtin_clz(x - 1));
  }
  T query(int l, int r){
    l = std::max(l, 0);
    r = std::min(r, n);
    if(l >= r) return id();
    int len = r - l;
    int b = msb(len);
    return merge(table[b][l], table[b][r - (1 << b)]);
  }
};

template<typename T, T (*merge)(T, T), T (*id)()>
struct sparse_table_memory{
private:
  static constexpr int block_size = 8;
  static constexpr int block_size_log = 3;
  sparse_table<T, merge, id> st;
  sparse_table_short<T, merge, id> st_short;
public:
  int n;
  sparse_table_memory(): n(0){}
  sparse_table_memory(const std::vector<T> &v): n(v.size()){
    if(v.empty()) return;
    int m = (n + block_size - 1) / block_size;
    std::vector<T> v2(m, id());
    for(int i = 0; i < m; i++){
      v2[i] = v[i * block_size];
      for(int j = (i * block_size) + 1, k = 1; j < n && k < block_size; j++, k++){
        v2[i] = merge(v2[i], v[j]);
      }
    }
    st = sparse_table<T, merge, id>(v2);
    st_short = sparse_table_short<T, merge, id>(v, block_size_log);
  }
  T query(int l, int r){
    l = std::max(l, 0);
    r = std::min(r, n);
    if(l >= r) return id();
    int lb = (l + block_size - 1) / block_size, rb = r / block_size;
    int lceil = lb * block_size, rfloor = rb * block_size;
    T ret = st.query(lb, rb);
    ret = merge(ret, st_short.query(l, std::min(r, lceil)));
    ret = merge(ret, st_short.query(std::max(l, rfloor), r));
    return ret;
  }
  void clear(){
    st.table.clear();
    st_short.clear();
  }
};

template<typename T>
struct rmq{
private:
  static constexpr T min_func(T a, T b){
    return std::min(a, b);
  }
  static constexpr T id(){
    return std::numeric_limits<T>::max();
  }
  static constexpr int block_size = 8;
  static constexpr int block_size_log = 3;
  sparse_table<T, min_func, id> st;
  sparse_table_short<T, min_func, id> st_short;
public:
  int n;
  rmq(): n(0){}
  rmq(const std::vector<T> &v): n(v.size()){
    if(v.empty()) return;
    int m = (n + block_size - 1) / block_size;
    std::vector<T> v2(m);
    for(int i = 0; i < m; i++){
      v2[i] = v[i * block_size];
      for(int j = (i * block_size) + 1, k = 1; j < n && k < block_size; j++, k++){
        v2[i] = min_func(v2[i], v[j]);
      }
    }
    st = sparse_table<T, min_func, id>(v2);
    st_short = sparse_table_short<T, min_func, id>(v, block_size_log);
  }
  T query(int l, int r){
    l = std::max(l, 0);
    r = std::min(r, n);
    if(l >= r) return id();
    int lb = (l + block_size - 1) / block_size, rb = r / block_size;
    int lceil = lb * block_size, rfloor = rb * block_size;
    T ret = st.query(lb, rb);
    ret = min_func(ret, st_short.query(l, std::min(r, lceil)));
    ret = min_func(ret, st_short.query(std::max(l, rfloor), r));
    return ret;
  }
  void clear(){
    st.table.clear();
    st_short.clear();
  }
};
#endif
