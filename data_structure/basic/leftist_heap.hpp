#ifndef _LEFTIST_HEAP_H_
#define _LEFTIST_HEAP_H_
#include <algorithm>
#include <vector>
#include <cassert>
template<typename Val>
struct leftist_heap{
private:
  struct node{
    node *l, *r;
    int s, sz;
    Val x;
    node(Val x): l(nullptr), r(nullptr), s(0), sz(1), x(x){}
  };
  node *root;
  static int size(node *v){
    return v ? v->sz : 0;
  }
  static bool is_null(node *v){
    return !v || !v->sz;
  }
  static node *meld_inner(node *a, node *b){
    if(is_null(a)) return b;
    if(is_null(b)) return a;
    if(a->x > b->x) std::swap(a, b);
    a->r = meld_inner(a->r, b);
    a->sz = 1 + size(a->l) + size(a->r);
    if(is_null(a->l) || a->l->s < a->r->s) std::swap(a->l, a->r);
    a->s = (is_null(a->r) ? 0 : a->r->s) + 1;
    return a;
  }
public:
  leftist_heap(): root(nullptr){}
  leftist_heap(Val x): root(new node(x)){}
  int size(){
    return (is_null(root) ? 0 : root->sz);
  }
  bool empty(){
    return size() == 0;
  }
  // マージして全要素を移動(永続でないためhの要素は全部消える)
  void meld(leftist_heap<Val> &h){
    root = meld_inner(root, h.root);
    h.root = nullptr;
  }
  Val min(){
    assert(!is_null(root));
    return root->x;
  }
  Val pop_min(){
    assert(!is_null(root));
    Val res = root->x;
    root = meld_inner(root->l, root->r);
    return res;
  }
  void push(Val x){
    root = meld_inner(root, new node(x));
  }
};


template<typename Val>
struct persistent_leftist_heap_iter;

template<typename Val>
struct persistent_leftist_heap{
private:
  struct node{
    node *l, *r;
    int s, sz;
    Val x;
    node(Val x): l(nullptr), r(nullptr), s(0), sz(1), x(x){}
  };
  static int size(node *v){
    return v ? v->sz : 0;
  }
  static bool is_null(node *v){
    return !v || !v->sz;
  }
  static node *copy_node(node *v){
    if(!v) return nullptr;
    return new node(*v);
  }
  static node *meld_inner(node *a, node *b){
    if(is_null(a)) return copy_node(b);
    if(is_null(b)) return copy_node(a);
    if(a->x > b->x) std::swap(a, b);
    a = copy_node(a);
    a->r = meld_inner(a->r, b);
    a->sz = 1 + size(a->l) + size(a->r);
    if(is_null(a->l) || a->l->s < a->r->s) std::swap(a->l, a->r);
    a->s = (is_null(a->r) ? 0 : a->r->s) + 1;
    return a;
  }
  static node *meld_build(node *a, node *b){
    if(is_null(a)) return b;
    if(is_null(b)) return a;
    if(a->x > b->x) std::swap(a, b);
    a->r = meld_inner(a->r, b);
    a->sz = 1 + size(a->l) + size(a->r);
    if(is_null(a->l) || a->l->s < a->r->s) std::swap(a->l, a->r);
    a->s = (is_null(a->r) ? 0 : a->r->s) + 1;
    return a;
  }
public:
  static node *build(const std::vector<Val> &v){
    node *res = nullptr;
    for(auto x : v) res = meld_build(res, new node(x));
    return res;
  }
  static node *make_heap(){
    return nullptr;
  }
  static node *make_heap(Val x){
    return new node(x);
  }
  static node *meld(node *a, node *b){
    return meld_inner(a, b);
  }
  static Val min(node *v){
    assert(size(v));
    return v->x;
  }
  static node *pop_min(node *v){
    assert(size(v));
    return meld(v->l, v->r);
  }
  static node *push(node *v, Val x){
    return meld(v, new node(x));
  }
  friend persistent_leftist_heap_iter<Val>;
};

template<typename Val>
struct persistent_leftist_heap_iter{
  using iter = persistent_leftist_heap_iter<Val>;
  using pheap = persistent_leftist_heap<Val>;
  using node = typename pheap::node;
  node *v;
  persistent_leftist_heap_iter(node *u): v(u){}
public:
  persistent_leftist_heap_iter(): v(nullptr){}
  persistent_leftist_heap_iter(Val x): v(pheap::make_heap(x)){}
  persistent_leftist_heap_iter(const std::vector<Val> &_v): v(pheap::build(_v)){}
  int size(){
    return pheap::size(v);
  }
  bool empty(){
    return size() == 0;
  }
  iter meld(iter &b){
    return pheap::meld(v, b.v);
  }
  Val min(){
    return pheap::min(v);
  }
  iter pop_min(){
    return pheap::pop_min(v);
  }
  iter push(Val x){
    return pheap::push(v, x);
  }
};

#endif