#ifndef _ERASABLE_HEAP_H_
#define _ERASABLE_HEAP_H_
#include "binary_heap.hpp"
#include "interval_heap.hpp"
#include <unordered_map>

template<typename Key, typename cmp>
struct erasable_binary_heap{
private:
  int n;
  std::unordered_map<Key, int> mp;
  binary_heap<Key, cmp> h;
public:
  erasable_binary_heap(): n(0){}
  erasable_binary_heap(const std::vector<Key> &_v): n(_v.size()), h(_v){}
  int size(){
    return n;
  }
  bool empty(){
    return !n;
  }
  void push(const Key &x){
    auto itr = mp.find(x);
    if(itr == mp.end()) mp.emplace(x, 1);
    else itr->second++;
    n++;
    h.push(x);
  }
  int count(const Key &x){
    auto itr = mp.find(x);
    return (itr == mp.end() ? 0 : itr->second);
  }
  void erase(const Key &x){
    auto itr = mp.find(x);
    if(itr != mp.end() && itr->second){
      itr->second--;
      n--;
    }
  }
  Key pop_min(){
    Key ans = h.pop_min();
    while(mp[ans] == 0) ans = h.pop_min();
    n--;
    mp[ans]--;
    return ans;
  }
  Key min(){
    while(mp[h.min()] == 0) h.pop_min();
    return h.min();
  }
};

template<typename Key>
struct erasable_interval_heap{
private:
  int n;
  std::unordered_map<Key, int> mp;
  interval_heap<Key> h;
public:
  erasable_interval_heap(): n(0){}
  erasable_interval_heap(const std::vector<Key> &_v): n(_v.size()), h(_v){}
  int size(){
    return n;
  }
  bool empty(){
    return !n;
  }
  void push(const Key &x){
    auto itr = mp.find(x);
    if(itr == mp.end()) mp.emplace(x, 1);
    else itr->second++;
    n++;
    h.push(x);
  }
  int count(const Key &x){
    auto itr = mp.find(x);
    return (itr == mp.end() ? 0 : itr->second);
  }
  void erase(const Key &x){
    auto itr = mp.find(x);
    if(itr != mp.end() && itr->second){
      itr->second--;
      n--;
    }
  }
  Key pop_min(){
    Key ans = h.pop_min();
    while(mp[ans] == 0) ans = h.pop_min();
    n--;
    mp[ans]--;
    return ans;
  }
  Key pop_max(){
    Key ans = h.pop_max();
    while(mp[ans] == 0) ans = h.pop_max();
    n--;
    mp[ans]--;
    return ans;
  }
  Key min(){
    while(mp[h.min()] == 0) h.pop_min();
    return h.min();
  }
  Key max(){
    while(mp[h.max()] == 0) h.pop_max();
    return h.max();
  }
};
#endif