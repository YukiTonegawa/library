#ifndef _SEGMENT_CONNECTIVITY_H_
#define _SEGMENT_CONNECTIVITY_H_
#include "../../data_structure/union_find/union_find.hpp"

#include <vector>
#include <iostream>
#include <cassert>
#include <numeric>
#include <algorithm>
#include <array>

template<typename Key>
struct treap_set{
  struct node{
    Key key, min, max;
    const unsigned int p;
    node *l, *r;
    node(Key key, unsigned int p): key(key), min(key), max(key), p(p), l(nullptr), r(nullptr){}
  };
  using pn = std::pair<node*, node*>;
  static unsigned int xor128(){
    static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 86675123;
    unsigned int t = (x ^ (x << 11));
    x = y, y = z, z = w;
    return (w = (w ^ (w >> 19)) ^ (t ^ ( t >> 8)));
  }
  inline void update(node *v){
    v->min = v->l ? v->l->min : v->key;
    v->max = v->r ? v->r->max : v->key;
  }
  node *make_node(Key k){
    return new node(k, xor128());
  }
  node *insert(node *v, node *n){
    if(!v) return n;
    if(v->p < n->p){
      pn s = split(v, n->key);
      n->l = s.first;
      n->r = s.second;
      update(n);
      return n;
    }else if(n->key < v->key){
      v->l = insert(v->l, n);
    }else{
      v->r = insert(v->r, n);
    }
    update(v);
    return v;
  }
  node *erase(node *v, Key k){
    if(!v) return nullptr;
    if(v->key == k){
      return merge(v->l, v->r);
    }else if(k < v->key){
      v->l = erase(v->l, k);
      update(v);
      return v;
    }else{
      v->r = erase(v->r, k);
      update(v);
      return v;
    }
  }
  node *merge(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    if(a->p > b->p){
      a->r = merge(a->r, b);
      update(a);
      return a;
    }else{
      b->l = merge(a, b->l);
      update(b);
      return b;
    }
  }
  // vをa, bに分割
  // aの任意の要素のkey < k かつ k <= bの任意の要素のkey
  pn split(node *v, Key k){
    if(!v) return {nullptr, nullptr};
    if(k <= v->key){
      pn s = split(v->l, k);
      v->l = s.second;
      update(v);
      return {s.first, v};
    }else{
      pn s = split(v->r, k);
      v->r = s.first;
      update(v);
      return {v, s.second};
    }
  }
  // [-inf, l), [l, r), [r, inf]
  std::array<node*, 3> split_range(node *v, Key l, Key r){
    if(!v) return {nullptr, nullptr, nullptr};
    if(v->key < l){
      auto s = split_range(v->r, l, r);
      v->r = s[0];
      update(v);
      return {v, s[1], s[2]};
    }else if(r <= v->key){
      auto s = split_range(v->l, l, r);
      v->l = s[2];
      update(v);
      return {s[0], s[1], v};
    }else{
      auto sl = split(v->l, l);
      auto sr = split(v->r, r);
      v->l = sl.second;
      v->r = sr.first;
      update(v);
      return {sl.first, v, sr.second};
    }
  }
};

struct range_set{
private:
  using treap_node = treap_set<int>::node;
  treap_set<int> t;
  struct node{
    int cmp;
    const unsigned int p;
    node *l, *r;
    treap_node *tn;
    node(int cmp, unsigned int p, treap_node *tn): cmp(cmp), p(p), l(nullptr), r(nullptr), tn(tn){}
    inline int ly(){return tn->min;}
    inline int ry(){return tn->max;}
  };
  using pn = std::pair<node*, node*>;
  static unsigned int xor128(){
    static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 86675123;
    unsigned int t = (x ^ (x << 11));
    x = y, y = z, z = w;
    return (w = (w ^ (w >> 19)) ^ (t ^ ( t >> 8)));
  }
  node *root = nullptr;
  node *make_node(int y, int cmp){
    return new node(cmp, xor128(), t.make_node(y));
  }
  node *make_node(treap_node *tn, int cmp){
    return new node(cmp, xor128(), tn);
  }
  // yを含むノードを探す
  node *find(node *v, int y){
    if(!v) return nullptr;
    if(v->ly() <= y && y <= v->ry()) return v;
    else if(y < v->ly()) return find(v->l, y);
    else return find(v->r, y);
  }
  pn split(node *v, int y){
    if(!v) return {nullptr, nullptr};
    if(y <= v->ly()){
      pn s = split(v->l, y);
      v->l = s.second;
      return {s.first, v};
    }else{
      pn s = split(v->r, y);
      v->r = s.first;
      return {v, s.second};
    }
  }
  node *merge(node *a, node *b){
    if(!a || !b) return !a ? b : a;
    if(a->p > b->p){
      a->r = merge(a->r, b);
      return a;
    }else{
      b->l = merge(a, b->l);
      return b;
    }
  }
  node *insert(node *v, node *n){
    if(!v) return n;
    if(v->p < n->p){
      pn s = split(v, n->ly());
      n->l = s.first;
      n->r = s.second;
      return n;
    }else if(n->ly() < v->ly()){
      v->l = insert(v->l, n);
    }else{
      v->r = insert(v->r, n);
    }
    return v;
  }
  node *erase(node *v, int y){
    if(!v) return nullptr;
    if(v->ly() <= y && y <= v->ry()){
      v->tn = t.erase(v->tn, y);
      if(!v->tn) return merge(v->l, v->r);
      else return v;
    }else if(y < v->ly()){
      v->l = erase(v->l, y);
      return v;
    }else{
      v->r = erase(v->r, y);
      return v;
    }
  }
  // [-inf, ly), [ly, ry), [ry, inf]で分割し, 真ん中の区間を全てマージ
  // 質問区間を完全に含むようなノードが存在する場合, 分割せずにtrueを返す
  std::tuple<node*, node*, node*, treap_node*, bool> enumerate_and_set(node *v, int ly, int ry, std::vector<int> &cmp){
    if(!v) return {nullptr, nullptr, nullptr, nullptr, false};
    if(v->ly() < ly && ry <= v->ry()){ // vが質問区間を含む
      auto [l, mid, r] = t.split_range(v->tn, ly, ry);
      if(mid) cmp.push_back(v->cmp);
      v->tn = l;
      node *rsp = make_node(r, v->cmp);
      return {v, rsp, nullptr, mid, true};
    }else if(v->ly() < ly){
      auto [l, mid, r, cross, f] = enumerate_and_set(v->r, ly, ry, cmp);
      if(f) return {v, mid, nullptr, cross, f};
      v->r = l;
      if(ly <= v->ry()){// 左側で交差
        cmp.push_back(v->cmp);
        auto [vtn, cross_left] = t.split(v->tn, ly);
        assert(vtn);
        v->tn = vtn;
        cross = t.merge(cross_left, cross);
      }
      return {v, mid, r, cross, false};
    }else if(ry <= v->ry()){
      auto [l, mid, r, cross, f] = enumerate_and_set(v->l, ly, ry, cmp);
      if(f) return {v, mid, nullptr, cross, f};
      v->l = r;
      if(v->ly() < ry){ // 右側で交差
        cmp.push_back(v->cmp);
        auto [cross_right, vtn] = t.split(v->tn, ry);
        assert(vtn);
        v->tn = vtn;
        cross = t.merge(cross, cross_right);
      }
      return {l, mid, v, cross, false};
    }else{ // vが[ly, ry)に完全に含まれる
      cmp.push_back(v->cmp);
      auto [l1, mid1, r1, cross_left, f1] = enumerate_and_set(v->l, ly, ry, cmp);
      auto [l2, mid2, r2, cross_right, f2] = enumerate_and_set(v->r, ly, ry, cmp);
      // assert(!r1 && !f1 && !l2 && !f2);
      // delete mid1
      // delete mid2
      return {l1, v, r2, t.merge(cross_left, t.merge(v->tn, cross_right)), false};
    }
  }
public:
  void insert(int y, int cmp){
    node *v = find(root, y);
    if(v){ // yをまたぐ区間がある場合分割
      auto [ls, rs] = t.split(v->tn, y);
      v->tn = ls;
      root = insert(root, make_node(rs, v->cmp));
    }
    root = insert(root, make_node(y, cmp));
  }
  void erase(int y){
    root = erase(root, y);
  }
  std::vector<int> enumerate_and_set(int ly, int ry, int cmp){
    std::vector<int> ret;
    auto [l, mid, r, cross, f] = enumerate_and_set(root, ly, ry, ret);
    if(f){
      if(cross) root = insert(root, make_node(cross, cmp));
      root = insert(root, mid);
      return ret;
    }
    if(mid){
      assert(cross);
      mid->l = mid->r = nullptr;
      mid->tn = cross;
    }else if(cross){
      mid = make_node(cross, cmp);
    }
    root = merge(l, merge(mid, r));
    return ret;
  }
};

struct point{
  int x, y;
};
struct line{
  point a, b;
  line(){}
};


// x軸またはy軸に平行な線分がm本ある
// find(a) : 番号aの線分の連結成分番号
// same(a, b): 線分aとbが同じ連結成分か
struct segment_connectivity{
  union_find uf;
  segment_connectivity(std::vector<line> L): uf(L.size()){
    struct query{
      char type; //0:縦線挿入, 1:横線, 2:縦線削除,
      int x, ly, ry, id;
    };
    std::vector<query> Q;
    int id = 0;
    for(line &l : L){
      if(l.a.x == l.b.x){
        Q.push_back(query{1, l.a.x, l.a.y, l.b.y, id});
      }else{
        Q.push_back(query{0, l.a.x, l.a.y, -1, id});
        Q.push_back(query{2, l.b.x, l.a.y, -1, id});
      }
      id++;
    }
    std::sort(Q.begin(), Q.end(), [&](query &a, query &b){
      if(a.x == b.x) return a.type < b.type;
      return a.x < b.x;
    });
    range_set rs;
    for(query q : Q){
      if(q.type == 0){
        rs.insert(q.ly, q.id);
      }else if(q.type == 2){
        rs.erase(q.ly);
      }else{
        auto S = rs.enumerate_and_set(q.ly, q.ry + 1, q.id);
        for(auto c : S) uf.unite(q.id, c);
      }
    }
  }
  int find(int a){
    return uf.find(a);
  }
  bool same(int a, int b){
    return uf.same(a, b);
  }
};
#endif
