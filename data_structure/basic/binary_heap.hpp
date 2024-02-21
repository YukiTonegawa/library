#ifndef _BINARY_HEAP_H_
#define _BINARY_HEAP_H_

#include <vector>
#include <algorithm>
template<typename Key, typename cmp = std::less<Key>>
struct binary_heap{
private:
  int n;
  std::vector<Key> v;
  void rotate_down(int k){
    while(true){
      int l = k * 2 + 1, r = l + 1;
      if(l >= n) break;
      if(r >= n){
        if(cmp()(v[l], v[k])) std::swap(v[k], v[l]), k = l;
        else break;
      }else{
        if(cmp()(v[l], v[k]) || cmp()(v[r], v[k])){
          if(cmp()(v[l], v[r])) std::swap(v[k], v[l]), k = l;
          else std::swap(v[k], v[r]), k = r;
        }else break;
      }
    }
  }
  void rotate_up(int k){
    while(k){
      k = (k - 1) / 2;
      int l = k * 2 + 1, r = l + 1;
      if(r >= n){
        if(cmp()(v[l], v[k])) std::swap(v[k], v[l]);
        else break;
      }else{
        if(cmp()(v[l], v[k]) || cmp()(v[r], v[k])){
          if(cmp()(v[l], v[r])) std::swap(v[k], v[l]);
          else std::swap(v[k], v[r]);
        }else break;
      }
    }
  }
  void heapify(const std::vector<Key> &_v){
    n = _v.size();
    v = _v;
    for(int i = (n - 1) / 2; i >= 0; i--) rotate_down(i);
  }
public:
  binary_heap(): n(0){}
  binary_heap(const std::vector<Key> &_v){heapify(_v);}
  void push(Key x){
    static constexpr int vector_init_size = 4;
    int vs = v.size();
    if(vs == n) v.resize(!vs ? vector_init_size : vs << 1, x);
    else v[n] = x;
    rotate_up(n++);
  }
  int size(){
    return n;
  }
  int empty(){
    return !n;
  }
  Key min(){
    assert(n);
    return v[0];
  }
  Key pop_min(){
    assert(n);
    Key ret = v[0];
    if(n == 1){
      n = 0;
    }else{
      v[0] = v[--n];
      rotate_down(0);
    }
    return ret;
  }
};
template<typename Key, typename Val, typename cmp>
struct binary_heap_map{
private:
  using p = std::pair<Key, Val>;
  int n;
  std::vector<p> v;
  void rotate_down(int k){
    while(true){
      int l = k * 2 + 1, r = l + 1;
      if(l >= n) break;
      if(r >= n){
        if(cmp()(v[l].first, v[k].first)) std::swap(v[k], v[l]), k = l;
        else break;
      }else{
        if(cmp()(v[l].first, v[k].first) || cmp()(v[r].first, v[k].first)){
          if(cmp()(v[l].first, v[r].first)) std::swap(v[k], v[l]), k = l;
          else std::swap(v[k], v[r]), k = r;
        }else break;
      }
    }
  }
  void rotate_up(int k){
    while(k){
      k = (k - 1) / 2;
      int l = k * 2 + 1, r = l + 1;
      if(r >= n){
        if(cmp()(v[l].first, v[k].first)) std::swap(v[k], v[l]);
        else break;
      }else{
        if(cmp()(v[l].first, v[k].first) || cmp()(v[r].first, v[k].first)){
          if(cmp()(v[l].first, v[r].first)) std::swap(v[k], v[l]);
          else std::swap(v[k], v[r]);
        }else break;
      }
    }
  }
  void heapify(const std::vector<p> &_v){
    n = _v.size();
    v = _v;
    for(int i = (n - 1) / 2; i >= 0; i--) rotate_down(i);
  }
public:
  binary_heap_map(): n(0){}
  binary_heap_map(const std::vector<p> &_v){heapify(_v);}
  void push(Key x_k, Val x_v){
    push(std::make_pair(x_k, x_v));
  }
  void push(p x){
    static constexpr int vector_init_size = 4;
    int vs = v.size();
    if(vs == n) v.resize(!vs ? vector_init_size : vs << 1, x);
    else v[n] = x;
    rotate_up(n++);
  }
  int size(){
    return n;
  }
  int empty(){
    return !n;
  }
  p min(){
    assert(n);
    return v[0];
  }
  p pop_min(){
    assert(n);
    p ret = v[0];
    if(n == 1){
      n = 0;
    }else{
      v[0] = v[--n];
      rotate_down(0);
    }
    return ret;
  }
};
#endif
