#ifndef _BINARY_TRIE_H_
#define _BINARY_TRIE_H_
#include <vector>
#include "../../traits.hpp"
template<typename Key, int KeyLog = 30, std::enable_if_t<is_intle64<Key>::value>* = nullptr>
struct binary_trie{
  static constexpr int bitlen = sizeof(Key) * 8;
private:
  struct node{
    Key k;
    int l, r, count;
    int ch[2];
    node(Key _k = 0, int _l = 0, int _r = 0):k(_k),l(_l),r(_r),count(0){
      ch[0] = ch[1] = -1;
    }
  };
  int sz = 0;
  Key xor_lazy;
  std::vector<node> nodes{node()};
  Key mask(int r){
    return ((1LL << (r + 1)) - 1);
  }
  int make_node(Key x, int l, int r){
    nodes.push_back(node(x, l, r));
    return ++sz;
  }
  int copy_node(int k){
    nodes.push_back(nodes[k]);
    return ++sz;
  }
  int next_bit(int v, Key k){
    Key x = (nodes[v].k ^ k) & mask(nodes[v].r);
    if(bitlen > 32) return (x ? bitlen - 1 - __builtin_clzll(x) : -1);
    return (x ? bitlen - 1 - __builtin_clz(x) : -1);
  }
  void __insert(Key x, int k = 1){
    x ^= xor_lazy;
    int v = 0, bit = KeyLog - 1;
    while(bit >= 0){
      nodes[v].count += k;
      int nxt = (x >> bit) & 1;
      if(nodes[v].ch[nxt]==-1){
        nodes[v].ch[nxt] = make_node(x, 0, bit);
        nodes[nodes[v].ch[nxt]].count = k;
        return;
      }
      v = nodes[v].ch[nxt];
      int diff_bit = std::max(nodes[v].l - 1, next_bit(v, x));
      if(diff_bit == nodes[v].l - 1){
        bit = diff_bit;
      }else{
        nxt = (x >> diff_bit) & 1;
        int tmp = copy_node(v);
        nodes[v].l = diff_bit + 1;
        nodes[tmp].r = diff_bit;
        nodes[v].count += k;
        nodes[v].ch[!nxt] = tmp;
        nodes[v].ch[nxt] = make_node(x, 0, diff_bit);
        nodes[nodes[v].ch[nxt]].count = k;
        return;
      }
    }
    nodes[v].count += k;
  }
  int __find(Key x, int v, int bit){
    v = nodes[v].ch[(x >> bit) & 1];
    if(!v) return 0;
    int diff_bit = std::max(nodes[v].l - 1, next_bit(v, x));
    if(diff_bit != nodes[v].l - 1) return 0;
    if(nodes[v].l == 0) return nodes[v].count;
    return __find(x, v, diff_bit);
  }
  int __erase(Key x, int k, int v, int bit){
    v = nodes[v].ch[(x >> bit) & 1];
    if(size(v)==0) return 0;
    int diff_bit = std::max(nodes[v].l - 1, next_bit(v, x));
    if(diff_bit > nodes[v].l - 1) return 0;
    if(nodes[v].l == 0){
      k = std::min(k, nodes[v].count);
      nodes[v].count -= k;
      return k;
    }
    k = __erase(x, k, v, diff_bit);
    nodes[v].count -= k;
    return k;
  }
  Key __kth_element(int k){
    assert(k < size());
    int bit = KeyLog - 1, v = 0;
    while(bit >= 0){
      int left = (xor_lazy >> bit) & 1;
      if(k < size(nodes[v].ch[left])){
        v = nodes[v].ch[left];
        bit = nodes[v].l - 1;
      }else{
        k -= size(nodes[v].ch[left]);
        v = nodes[v].ch[!left];
        bit = nodes[v].l - 1;
      }
    }
    return nodes[v].k ^ xor_lazy;
  }
  int __lower_count(Key k){
    int bit = KeyLog - 1, ret = 0, v = 0;
    Key kx = k ^ xor_lazy;
    while(bit >= 0){
      if(v == -1) return ret;
      if((k >> bit) & 1){
        ret += size(nodes[v].ch[(xor_lazy >> bit) & 1]);
        v = nodes[v].ch[!((xor_lazy >> bit) & 1)];
      }else{
        v = nodes[v].ch[(xor_lazy >> bit) & 1];
      }
      if(v != -1){
        int diff_bit = std::max(nodes[v].l - 1, next_bit(v, kx));
        if(diff_bit != nodes[v].l - 1){
          ret += ((k >> diff_bit) & 1 ? size(v) : 0);
          return ret;
        }
        bit = nodes[v].l - 1;
      }
    }
    return ret;
  }
public:
  binary_trie(){}
  int size(int v = 0){
    return (v == -1 ? 0 : nodes[v].count);
  }
  void insert(Key x, int k = 1){
    __insert(x, k);
  }
  void xor_all(Key x){
    xor_lazy ^= x;
  }
  int find(Key x){
    return __find(x ^ xor_lazy, 0, KeyLog-1);
  }
  int erase(Key x, int k = std::numeric_limits<int>::max()){
    k = __erase(x ^ xor_lazy, k, 0, KeyLog - 1);
    nodes[0].count -= k;
    return k;
  }
  Key kth_element(int k){
    return __kth_element(k);
  }
  int lower_count(Key k){
    return __lower_count(k);
  }
};
#endif