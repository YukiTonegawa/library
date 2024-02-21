#ifndef _PERSISTENT_ARRAY_H_
#define _PERSISTENT_ARRAY_H_
#include <vector>
#include <array>
#include <cassert>


template<typename Val, typename Idx = int>
struct persistent_array_iter;

template<typename Val, typename Idx = int>
struct persistent_array{
  constexpr static int Split = 8, Mask = Split - 1, Shift = 31 -__builtin_clz(Split);
private:
  struct node{
    Val val;
    std::array<int, Split> ch;
    node(){ch.fill(-1);}
    node(Val _val) : val(_val){ch.fill(-1);}
  };
  Val z;
  std::vector<int> version{-1};
  std::vector<node> nodes;
  int make_node(Val val){
    nodes.push_back(node(val));
    return (int)nodes.size() - 1;
  }
  int copy_node(int k){
    if(k == -1) nodes.push_back(node(z));
    else nodes.push_back(nodes[k]);
    return (int)nodes.size() - 1;
  }
  int inner_set(int par, Idx idx, Val x){
    int v = copy_node(par);
    if(!idx) nodes[v].val = x;
    else nodes[v].ch[idx & Mask] = inner_set(nodes[v].ch[idx & Mask], idx >> Shift, x);
    return v;
  }
  void init(const std::vector<Val> &v, Idx idx, int k, int s, int las){
    if(k == 0 || las) nodes[k].val = v[idx];
    Idx add = (Idx)1 << s;
    for(int i = 0; i < Split; i++){
      Idx nxt = idx + i * add;
      if(nxt < v.size() && add < v.size()){
        nodes[k].ch[i] = make_node(z);
        init(v, nxt, nodes[k].ch[i], s + Shift, i);
      }
    }
  }
  void init(const std::vector<Val> &v){
    assert(nodes.empty() && v.size());
    make_node(z);
    version[0] = 0;
    init(v, 0, 0, 0, 1);
  }
  int set(int ver, Idx k, Val x){
    version.push_back(inner_set(version[ver], k, x));
    return (int)version.size() - 1;
  }
  Val get(int ver, Idx k){
    int v = version[ver];
    if(v == -1) return z;
    while(k){
      v = nodes[v].ch[k & Mask];
      k >>= Shift;
      if(v == -1) return z;
    }
    return nodes[v].val;
  }
  persistent_array(){}
  persistent_array(Val z): z(z){}
  friend persistent_array_iter<Val, Idx>;
};

template<typename Val, typename Idx>
struct persistent_array_iter{

  int id;
  persistent_array<Val, Idx> *p;
  persistent_array_iter(int id, persistent_array<Val, Idx> *p): id(id), p(p){}
public:
  persistent_array_iter(){}
  persistent_array_iter(Val z): id(0){
    p = new persistent_array<Val, Idx>(z);
  }
  persistent_array_iter(Val z, const std::vector<Val> &v): id(0){
    p = new persistent_array<Val, Idx>(z);
    p->init(v);
  }
  persistent_array_iter<Val, Idx> set(Idx idx, Val x){
    return persistent_array_iter<Val, Idx>(p->set(id, idx, x), p);
  }
  Val get(Idx idx){
    return p->get(id, idx);
  }
};
/*
#include <iostream>
int main(){
  std::cin.tie(nullptr);
  std::ios::sync_with_stdio(false);
  int q;
  std::cin >> q;
  persistent_array_iter<long long, long long> pai(0);

  for(int i = 0; i < q; i++){
    int t;
    std::cin >> t;
    if(t == 0){
      long long k, x;
      std::cin >> k >> x;
      pai = pai.set(k, x);
    }else{
      long long k;
      std::cin >> k;
      std::cout << pai.get(k) << '\n';
    }
  }
  return 0;
}
*/
#endif