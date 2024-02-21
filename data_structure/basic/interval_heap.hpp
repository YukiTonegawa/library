#ifndef _INTERVAL_HEAP_H_
#define _INTERVAL_HEAP_H_

#include <vector>

template<typename Key>
struct interval_heap{
private:
  inline int parent(int k){
    return ((k >> 1) - 1) & ~1;
  }
  int n;
  std::vector<Key> v;
public:
  interval_heap(): n(0){}
  interval_heap(const std::vector<Key> &_v): n(_v.size()), v(_v){
    make_heap();
  }
  int size(){
    return n;
  }
  bool empty(){
    return !n;
  }
  void push(const Key &x){
    static constexpr int vector_init_size = 4;
    int vs = v.size();
    if(vs == n) v.resize(!vs ? vector_init_size : vs << 1, x);
    else v[n] = x;
    up(n++);
  }
  Key pop_min(){
    Key res = min();
    if(n < 3){
      n--;
	  }else{
	    std::swap(v[1], v[--n]);
	    int k = down(1);
	    up(k);
    }
    return res;
  }
  Key pop_max(){
    Key res = max();
    if(n < 2){
      n--;
    }else{
	    std::swap(v[0], v[--n]);
	    int k = down(0);
	    up(k);
    }
    return res;
  }
  Key min(){
    return n < 2 ? v[0] : v[1];
  }
  Key max(){
    return v[0];
  }
  void make_heap(){
    for(int i = n; i--; ){
      if (i & 1 && v[i - 1] < v[i]) std::swap(v[i - 1], v[i]);
      int k = down(i);
      up(k, i);
    }
  }
  int down(int k){
    if(k & 1){
	    while(2 * k + 1 < n){
        int c = 2 * k + 3;
        if(n <= c || v[c - 2] < v[c]) c -= 2;
        if(c < n && v[c] < v[k]){
          std::swap(v[k], v[c]);
          k = c;
        }else break;
	    }
    }else{
      while(2 * k + 2 < n) {
        int c = 2 * k + 4;
        if (n <= c || v[c] < v[c - 2]) c -= 2;
        if (c < n && v[k] < v[c]){
          std::swap(v[k], v[c]);
          k = c;
        }else break;
      }
    }
    return k;
  }
  int up(int k, int root = 1){
    if((k | 1) < n && v[k & ~1] < v[k | 1]) {
      std::swap(v[k & ~1], v[k | 1]);
      k ^= 1;
    }
    int p;
    while(root < k && v[p = parent(k)] < v[k]){
      std::swap(v[p], v[k]);
      k = p;
    }
    while(root < k && v[k] < v[p = parent(k) | 1]){
      std::swap(v[p], v[k]);
      k = p;
    }
    return k;
  }
};

template<typename Key, typename Val>
struct interval_heap_map{
private:
  inline int parent(int k){
    return ((k >> 1) - 1) & ~1;
  }
  using p = std::pair<Key, Val>;
  int n;
  std::vector<p> v;
public:
  interval_heap_map(): n(0){}
  interval_heap_map(const std::vector<p> &_v): n(_v.size()), v(_v){
    make_heap();
  }
  int size(){
    return n;
  }
  bool empty(){
    return !n;
  }
  void push(const p &x){
    static constexpr int vector_init_size = 4;
    int vs = v.size();
    if(vs == n) v.resize(!vs ? vector_init_size : vs << 1, x);
    else v[n] = x;
    up(n++);
  }
  p pop_min(){
    p res = min();
    if(n < 3){
      n--;
	  }else{
	    std::swap(v[1], v[--n]);
	    int k = down(1);
	    up(k);
    }
    return res;
  }
  p pop_max(){
    p res = max();
    if(n < 2){
      n--;
    }else{
	    std::swap(v[0], v[--n]);
	    int k = down(0);
	    up(k);
    }
    return res;
  }
  p min(){
    return n < 2 ? v[0] : v[1];
  }
  p max(){
    return v[0];
  }
  void make_heap(){
    for(int i = n; i--; ){
      if (i & 1 && v[i - 1] < v[i]) std::swap(v[i - 1], v[i]);
      int k = down(i);
      up(k, i);
    }
  }
  int down(int k){
    if(k & 1){
	    while(2 * k + 1 < n){
        int c = 2 * k + 3;
        if(n <= c || v[c - 2] < v[c]) c -= 2;
        if(c < n && v[c] < v[k]){
          std::swap(v[k], v[c]);
          k = c;
        }else break;
	    }
    }else{
      while(2 * k + 2 < n) {
        int c = 2 * k + 4;
        if (n <= c || v[c] < v[c - 2]) c -= 2;
        if (c < n && v[k] < v[c]){
          std::swap(v[k], v[c]);
          k = c;
        }else break;
      }
    }
    return k;
  }
  int up(int k, int root = 1){
    if((k | 1) < n && v[k & ~1] < v[k | 1]) {
      std::swap(v[k & ~1], v[k | 1]);
      k ^= 1;
    }
    int p;
    while(root < k && v[p = parent(k)] < v[k]){
      std::swap(v[p], v[k]);
      k = p;
    }
    while(root < k && v[k] < v[p = parent(k) | 1]){
      std::swap(v[p], v[k]);
      k = p;
    }
    return k;
  }
};
#endif
