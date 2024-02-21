#ifndef _SCAPEGOAT_TREE_H_
#define _SCAPEGOAT_TREE_H_
#include <vector>
#include "../algebraic_structure/monoid.hpp"

template<typename Key, typename monoid>
struct scapegoat_tree{
  using Val = typename monoid::Val;
  static constexpr auto id = monoid::id;
  static constexpr auto merge_val = monoid::merge;
private:
  static constexpr float alpha = 0.6;
  struct node{
    bool is_alive;
    int size;
    node *l, *r;
    Key key;
    Val val, sum;
    node(Key key, Val val): is_alive(true), size(1), l(nullptr), r(nullptr), key(key), val(val), sum(val){}
  };
  node *make_node(Key key, Val val = id()){return new node(key);}
  int max_node_count;
  node *root;
  inline int size(node *v){return !v ? 0 : v->size;}
  inline bool check_alpha_weight_balanced(node *v){
    int threshold = alpha * v->size;
    return size(v->l) <= threshold && size(v->r) <= threshold;
  }
  void update(node *v){
    v->size = v->is_alive;
    v->sum = v->is_alive ? v->val : id();
    if(v->l){
      v->size += v->l->size;
      v->sum = merge_val(v->l->sum, v->sum);
    }
    if(v->r){
      v->size += v->r->size;
      v->sum = merge_val(v->sum, v->r->sum);
    }
  }
  void subtree_enumerate(node *v, std::vector<node*> &nodes, int &cnt){
    if(v->l && v->l->size) subtree_enumerate(v->l, nodes, cnt);
    if(v->is_alive) nodes[cnt++] = v;
    if(v->r && v->r->size) subtree_enumerate(v->r, nodes, cnt);
  }
  node *build(int l, int r, const std::vector<node*> &nodes){
    int mid = (l + r) / 2;
    node *v = nodes[mid];
    if(l < mid) v->l = build(l, mid, nodes);
    else v->l = nullptr;
    if(mid + 1 < r) v->r = build(mid + 1, r, nodes);
    else v->r = nullptr;
    update(v);
    return v;
  }
  node* rebalance(node *v){
    if(!v || !v->size) return nullptr;
    std::vector<node*> ret(v->size);
    int cnt = 0;
    subtree_enumerate(v, ret, cnt);
    return build(0, v->size, ret);
  }
  node *insert(node *v, Key key, Val val){
    if(!v) return make_node(key, val);
    if(v->key == key) return v;
    else if(v->key < key) v->r = insert(v->r, key, val);
    else v->l = insert(v->l, key, val);
    if(!check_alpha_weight_balanced(v)) return rebalance(v);
    update(v);
    return v;
  }
  void erase(node *v, Key key){
    if(!v) return;
    if(v->key == key) v->is_alive = false;
    else if(v->key < key) erase(v->r, key);
    else erase(v->l, key);
    update(v);
  }
public:
  scapegoat_tree(): max_node_count(0), root(nullptr){}
  int size(){
    return size(root);
  }
  void insert(Key key, Val val){
    root = insert(root, key, val);
    max_node_count = std::max(max_node_count, size(root));
  }
  void erase(Key key){
    erase(root, key);
    if(size(root) <= alpha * max_node_count) root = rebalance(root);
  }
};
#endif
