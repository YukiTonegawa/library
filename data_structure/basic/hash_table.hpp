#ifndef _HASH_TABLE_H_
#define _HASH_TABLE_H_
#include <vector>
#include <cassert>
#include <limits>
#include "../../traits.hpp"
#include "../../misc/random_number.hpp"
#include "../../misc/ceillog2.hpp"

// uint32 : uKey ∈ [0, 2^31 - 1)
// uint64 : uKey ∈ [0, 2^63 - 1)
template<typename uKey>
struct uhash_set{
  static_assert(is_unsigned_intle64<uKey>::value, "uKey must be unsigned integer (32 or 64bit)");
private:
  using Key = typename std::make_signed<uKey>::type;
  static constexpr int initial_size = 8; // 指定しなかった場合の初期サイズ
  static constexpr float alpha = 0.7; // 占有率がこれを超えると再構築
  static constexpr Key null_val = std::numeric_limits<Key>::min();
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static const unsigned long long r;
  int table_size, table_size_mod;
  int occupied_cnt, elem_cnt; // 使われている要素の数, そのうちまだ消されていない要素の数
  std::vector<Key> v;
  size_t hash(unsigned long long x){
    x += r;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return (x ^ (x >> 31)) & table_size_mod;
  }
  void expand(){
    std::vector<Key> tmp(elem_cnt, null_val);
    for(int i = 0, j = 0; i < table_size; i++) if(v[i] >= 0) tmp[j++] = v[i];
    occupied_cnt = elem_cnt = 0;
    table_size <<= 1;
    table_size_mod = table_size - 1;
    v = std::vector<Key>(table_size, null_val);
    for(int i = 0; i < tmp.size(); i++) insert_inner(tmp[i]);
  }
  void insert_inner(Key x){
    if(occupied_cnt > table_size * alpha) expand();
    int i = hash(x);
    while(true){
      if(v[i] == null_val){
        v[i] = x;
        occupied_cnt++;
        elem_cnt++;
        return;
      }else if(v[i] == x || v[i] == -x - 1){
        elem_cnt += (v[i] == -x - 1);
        v[i] = x;
        return;
      }
      i = (i + 1) & table_size_mod;
    }
  }
  void erase_inner(Key x){
    int i = hash(x);
    while(true){
      if(v[i] == null_val) return;
      else if(v[i] == x || v[i] == -x - 1){
        elem_cnt -= (v[i] == x);
        v[i] = -x - 1;
        return;
      }
      i = (i + 1) & table_size_mod;
    }
  }
  bool find_inner(Key x){
    int i = hash(x);
    while(true){
      if(v[i] == null_val) return false;
      else if(v[i] == x || v[i] == -x - 1) return v[i] == x;
      i = (i + 1) & table_size_mod;
    }
  }
public:
  uhash_set(): table_size(initial_size), table_size_mod(table_size - 1), occupied_cnt(0), elem_cnt(0), v(table_size, null_val){}
  uhash_set(int _sz): table_size(1 << ceillog2(_sz)), table_size_mod(table_size - 1), occupied_cnt(0), elem_cnt(0), v(table_size, null_val){}
  // 追加, すでにある場合は無視
  void insert(uKey x){
    assert(x < inf);
    insert_inner(x);
  }
  // 削除, 無い場合は無視
  void erase(uKey x){
    assert(x < inf);
    erase_inner(x);
  }
  // 検索
  bool find(uKey x){
    assert(x < inf);
    return find_inner(x);
  }
  int size(){
    return elem_cnt;
  }
  bool empty(){
    return elem_cnt == 0;
  }
  std::vector<uKey> enumerate(){
    std::vector<uKey> res;
    for(int i = 0; i < table_size; i++) if(v[i] >= 0) res.push_back(v[i]);
    return res;
  }
  void clear(){
    table_size = initial_size;
    table_size_mod = table_size - 1;
    occupied_cnt = elem_cnt = 0;
    v = std::vector<Key>(table_size, null_val);
  }
};
template<typename uKey>
const unsigned long long uhash_set<uKey>::r = random_number();


// uint32 : uKey ∈ [0, 2^31 - 1)
// uint64 : uKey ∈ [0, 2^63 - 1)
template<typename uKey, typename Val>
struct uhash_map{
  static_assert(is_unsigned_intle64<uKey>::value, "uKey must be unsigned integer (32 or 64bit)");
private:
  using Key = typename std::make_signed<uKey>::type;
  static constexpr int initial_size = 8; // 指定しなかった場合の初期サイズ
  static constexpr float alpha = 0.7; // 占有率がこれを超えると再構築
  static constexpr Key null_val = std::numeric_limits<Key>::min();
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static const unsigned long long r;
  int table_size, table_size_mod;
  int occupied_cnt, elem_cnt; // 使われている要素の数, そのうちまだ消されていない要素の数
  std::vector<std::pair<Key, Val>> v;
  size_t hash(unsigned long long x){
    x += r;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return (x ^ (x >> 31)) & table_size_mod;
  }
  void expand(){
    std::vector<std::pair<Key, Val>> tmp(elem_cnt, {null_val, Val()});
    for(int i = 0, j = 0; i < table_size; i++) if(v[i].first >= 0) tmp[j++] = v[i];
    occupied_cnt = elem_cnt = 0;
    table_size <<= 1;
    table_size_mod = table_size - 1;
    v = std::vector<std::pair<Key, Val>>(table_size, {null_val, Val()});
    for(int i = 0; i < tmp.size(); i++) emplace_inner(tmp[i].first, tmp[i].second);
  }
  void emplace_inner(Key x, Val y){
    if(occupied_cnt > table_size * alpha) expand();
    int i = hash(x);
    while(true){
      if(v[i].first == null_val){
        v[i] = {x, y};
        occupied_cnt++;
        elem_cnt++;
        return;
      }else if(v[i].first == x || v[i].first == -x - 1){
        if(v[i].first == -x - 1){
          elem_cnt++;
          v[i] = {x, y};
        }
        return;
      }
      i = (i + 1) & table_size_mod;
    }
  }
  void emplace_replace_inner(Key x, Val y){
    if(occupied_cnt > table_size * alpha) expand();
    int i = hash(x);
    while(true){
      if(v[i].first == null_val){
        v[i] = {x, y};
        occupied_cnt++;
        elem_cnt++;
        return;
      }else if(v[i].first == x || v[i].first == -x - 1){
        elem_cnt += (v[i].first == -x - 1);
        v[i] = {x, y};
        return;
      }
      i = (i + 1) & table_size_mod;
    }
  }
  void erase_inner(Key x){
    int i = hash(x);
    while(true){
      if(v[i].first == null_val) return;
      else if(v[i].first == x || v[i].first == -x - 1){
        elem_cnt -= (v[i].first == x);
        v[i].first = -x - 1;
        return;
      }
      i = (i + 1) & table_size_mod;
    }
  }
  bool find_inner(Key x){
    int i = hash(x);
    while(true){
      if(v[i].first == null_val) return false;
      else if(v[i].first == x || v[i].first == -x - 1) return v[i].first == x;
      i = (i + 1) & table_size_mod;
    }
  }
  std::pair<bool, Val> at_inner(Key x){
    int i = hash(x);
    while(true){
      if(v[i].first == null_val) return std::make_pair(false, v[i].second);
      else if(v[i].first == x || v[i].first == -x - 1) return (v[i].first == x ? std::make_pair(true, v[i].second) : std::make_pair(false, v[i].second));
      i = (i + 1) & table_size_mod;
    }
  }
public:
  uhash_map(): table_size(initial_size), table_size_mod(table_size - 1), occupied_cnt(0), elem_cnt(0), v(table_size, {null_val, Val()}){}
  uhash_map(int _sz): table_size(1 << ceillog2(_sz)), table_size_mod(table_size - 1), occupied_cnt(0), elem_cnt(0), v(table_size, {null_val, Val()}){}
  // 追加, すでにある場合は無視
  void emplace(uKey x, Val y){
    assert(x < inf);
    emplace_inner(x, y);
  }
  // 追加, すでにある場合は置き換える
  void emplace_replace(uKey x, Val y){
    assert(x < inf);
    emplace_replace_inner(x, y);
  }
  // 削除, 無い場合は無視
  void erase(uKey x){
    assert(x < inf);
    erase_inner(x);
  }
  // 検索
  bool find(uKey x){
    assert(x < inf);
    return find_inner(x);
  }
  // 存在するか, 存在する場合その値
  std::pair<bool, Val> at(uKey x){
    assert(x < inf);
    return at_inner(x);
  }
  int size(){
    return elem_cnt;
  }
  bool empty(){
    return elem_cnt == 0;
  }
  std::vector<std::pair<uKey, Val>> enumerate(){
    std::vector<std::pair<uKey, Val>> res;
    for(int i = 0; i < table_size; i++) if(v[i].first >= 0) res.push_back(v[i]);
    return res;
  }
  void clear(){
    table_size = initial_size;
    table_size_mod = table_size - 1;
    occupied_cnt = elem_cnt = 0;
    v = std::vector<std::pair<Key, Val>>(table_size, {null_val, Val()});
  }
};
template<typename uKey, typename Val>
const unsigned long long uhash_map<uKey, Val>::r = random_number();
#endif