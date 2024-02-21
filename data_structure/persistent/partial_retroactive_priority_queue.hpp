#ifndef _PARTIAL_RETROACTIVE_PRIORITY_QUEUE_H_
#define _PARTIAL_RETROACTIVE_PRIORITY_QUEUE_H_
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
template<typename Key, typename KeySum>
struct partial_retroactive_priority_queue{
public:
  static constexpr Key inf = std::numeric_limits<Key>::max();
  static constexpr Key minf = std::numeric_limits<Key>::min();
private:
  struct node{
    node *l, *r, *p;
    int h, cnt_erased, cnt_current, cnt_pop, diffmin;
    __int8_t mode; // 0: 消されたノード, 1: 消されていないノード, 2: popノード　
    Key key;
    KeySum sum_current;
    node *max_erased, *min_current;
    // insert用
    node(Key _key): l(nullptr), r(nullptr), p(nullptr), h(1), cnt_erased(0), cnt_current(1), cnt_pop(0), diffmin(0), mode(1), key(_key), sum_current(_key), min_current(this){}
    // pop用
    node(): l(nullptr), r(nullptr), p(nullptr), h(1), cnt_erased(0), cnt_current(0), cnt_pop(1), diffmin(-1), mode(2), key(inf), sum_current(0){}
    int balanace_factor(){return (l ? l->h : 0) - (r ? r->h : 0);}
  };
  node *dummy_max, *dummy_min, *tmp_node, *root;

  // pop用
  node *make_node(){
    node *res = new node();
    res->max_erased = dummy_max;
    res->min_current = dummy_min;
    return res;
  }
  // insert用
  node *make_node(Key _key){
    node *res = new node(_key);
    res->max_erased = dummy_max;
    return res;
  }
  int size_subtree(node *v){
    return !v ? 0 : v->cnt_erased + v->cnt_current + v->cnt_pop;
  }
  void update(node *v){
    v->h = std::max(v->l ? v->l->h : 0,  v->r ? v->r->h : 0) + 1;
    if(v->mode == 0){
      // 消された
      v->cnt_erased = 1;
      v->cnt_current = 0;
      v->cnt_pop = 0;
      v->diffmin = 1;
      v->sum_current = 0;
      v->max_erased = v;
      v->min_current = dummy_min;
    }else if(v->mode == 1){
      // 存在
      v->cnt_erased = 0;
      v->cnt_current = 1;
      v->cnt_pop = 0;
      v->diffmin = 0;
      v->sum_current = v->key;
      v->max_erased = dummy_max;
      v->min_current = v;
    }else{
      // popクエリ
      v->cnt_erased = 0;
      v->cnt_current = 0;
      v->cnt_pop = 1;
      v->diffmin = -1;
      v->sum_current = 0;
      v->max_erased = dummy_max;
      v->min_current = dummy_min;
    }
    if(v->l){
      v->cnt_erased += v->l->cnt_erased;
      v->cnt_current += v->l->cnt_current;
      v->cnt_pop += v->l->cnt_pop;
      v->diffmin = std::min(v->l->diffmin, v->l->cnt_erased - v->l->cnt_pop + v->diffmin);
      v->sum_current += v->l->sum_current;
      if(v->max_erased->key < v->l->max_erased->key) v->max_erased = v->l->max_erased;
      if(v->min_current->key > v->l->min_current->key) v->min_current = v->l->min_current;
    }
    if(v->r){
      v->diffmin = std::min(v->diffmin, v->cnt_erased - v->cnt_pop + v->r->diffmin);
      v->cnt_erased += v->r->cnt_erased;
      v->cnt_current += v->r->cnt_current;
      v->cnt_pop += v->r->cnt_pop;
      v->sum_current += v->r->sum_current;
      if(v->max_erased->key < v->r->max_erased->key) v->max_erased = v->r->max_erased;
      if(v->min_current->key > v->r->min_current->key) v->min_current = v->r->min_current;
    }
  }
  node *rotate_right(node *v){
    node *l = v->l;
    if((v->l = l->r)) v->l->p = v;
    if((l->r = v)) l->r->p = l;
    update(v);
    update(l);
    return l;
  }
  node *rotate_left(node *v){
    node *r = v->r;
    if((v->r = r->l)) v->r->p = v;
    if((r->l = v)) r->l->p = r;
    update(v);
    update(r);
    return r;
  }
  node *balance(node *v){
    int bf = v->balanace_factor();
    assert(-2 <= bf && bf <= 2);
    if(bf == 2){
      if(v->l->balanace_factor() == -1){
        if((v->l = rotate_left(v->l))) v->l->p = v;
        update(v);
      }
      return rotate_right(v);
    }else if(bf == -2){
      if(v->r->balanace_factor() == 1){
        if((v->r = rotate_right(v->r))) v->r->p = v;
        update(v);
      }
      return rotate_left(v);
    }
    return v;
  }
  node *cut_leftmost(node *v){
    if(v->l){
      if((v->l = cut_leftmost(v->l))) v->l->p = v;
      update(v);
      return balance(v);
    }
    tmp_node = v;
    return v->r;
  }
  // k番目にuを追加
  node *insert_inner(node *v, int k, node *u){
    if(!v) return u;
    int szl = size_subtree(v->l);
    if(k <= szl){
      if((v->l = insert_inner(v->l, k, u))) v->l->p = v;
    }else if(k > szl){
      if((v->r = insert_inner(v->r, k - szl - 1, u))) v->r->p = v;
    }
    update(v);
    return balance(v);
  }
  // k番目のノードを消してそのノードをtmp_nodeに格納する
  node *erase_inner(node *v, int k){
    assert(v);
    int szl = size_subtree(v->l);
    if(k < szl){
      if((v->l = erase_inner(v->l, k))) v->l->p = v;
    }else if(k > szl){
      if((v->r = erase_inner(v->r, k - szl - 1))) v->r->p = v;
    }else{
      if(v->r){
        if((v->r = cut_leftmost(v->r))) v->r->p = v;
        if((tmp_node->l = v->l)) v->l->p = tmp_node;
        if((tmp_node->r = v->r)) v->r->p = tmp_node;
        std::swap(v, tmp_node);
        update(v);
        return balance(v);
      }
      tmp_node = v;
      return v->l;
    }
    update(v);
    return balance(v);
  }
  // sum[0, r]が最小となるようなrのうち最左
  int leftmost_diffmin(node *v){
    static constexpr int inf_diff = 1 << 30;
    int lsum = (v->l ? v->l->cnt_erased - v->l->cnt_pop : 0);
    int ldiffmin = (v->l ? v->l->diffmin : inf_diff);
    int rdiffmin = (v->r ? lsum - (v->mode - 1) + v->r->diffmin : inf_diff);
    if(ldiffmin <= rdiffmin && ldiffmin <= lsum - (v->mode - 1)){
      return leftmost_diffmin(v->l);
    }else if(lsum - (v->mode - 1) <= rdiffmin){
      return size_subtree(v->l);
    }else{
      return size_subtree(v->l) + 1 + leftmost_diffmin(v->r);
    }
  }
  // sum[0, r]が最小となるようなrのうち最右
  int rightmost_diffmin(node *v){
    static constexpr int inf_diff = 1 << 30;
    int lsum = (v->l ? v->l->cnt_erased - v->l->cnt_pop : 0);
    int ldiffmin = (v->l ? v->l->diffmin : inf_diff);
    int rdiffmin = (v->r ? lsum - (v->mode - 1) + v->r->diffmin : inf_diff);
    if(rdiffmin <= ldiffmin && rdiffmin <= lsum - (v->mode - 1)){
      return size_subtree(v->l) + 1 + rightmost_diffmin(v->r);
    }else if(lsum - (v->mode - 1) <= ldiffmin){
      return size_subtree(v->l);
    }else{
      return rightmost_diffmin(v->l);
    }
  }
  // l以降の消された要素のうちの最大値を探す
  node* find_max_erased_suffix(node *v, int l){
    if(!v) return dummy_max;
    int szl = size_subtree(v->l);
    if(szl + 1 <= l){
      return find_max_erased_suffix(v->r, l - szl - 1);
    }else{
      node *lmax = find_max_erased_suffix(v->l, l);
      node *rmax = v->r ? v->r->max_erased : dummy_max;
      if(v->mode == 0 && v->key >= std::max(lmax->key, rmax->key)) return v;
      return (lmax->key > rmax->key ? lmax : rmax);
    }
  }
  // [0, r)の最小ノードを消す
  node *find_min_current_prefix(node *v, int r){
    if(r <= 0 || !v) return dummy_min;
    int szl = size_subtree(v->l);
    if(r <= szl){
      return find_min_current_prefix(v->l, r);
    }else{
      node *lmin = v->l ? v->l->min_current : dummy_min;
      node *rmin = find_min_current_prefix(v->r, r - szl - 1);
      if(v->mode == 1 && v->key <= std::min(lmin->key, rmin->key)) return v;
      return (lmin->key < rmin->key ? lmin : rmin);
    }
  }
public:
  partial_retroactive_priority_queue(): root(nullptr){
    dummy_max = make_node(minf);
    dummy_min = make_node(inf);
    tmp_node = make_node();
  }
  int len_sequence(){
    return size_subtree(root);
  }
  void insert_op_insert(int k, Key x){
    assert(0 <= k && k <= len_sequence());
    node *u = make_node(x);
    u->mode = 0;
    update(u);
    if((root = insert_inner(root, k, u))) root->p = nullptr;
    int l = (root->diffmin == 0 ? rightmost_diffmin(root) : 0);
    u = find_max_erased_suffix(root, l);
    u->mode = 1;
    update(u);
    while(u->p){
      u = u->p;
      update(u);
    }
  }
  void insert_op_popmin(int k){
    assert(0 <= k && k <= len_sequence());
    node *v = make_node();
    if((root = insert_inner(root, k, v))) root->p = nullptr;
    int r = leftmost_diffmin(root);
    v = find_min_current_prefix(root, r);
    assert(v->mode == 1);
    v->mode = 0;
    update(v);
    while(v->p){
      v = v->p;
      update(v);
    }
  }
  // k番目の操作を削除
  void erase_op(int k){
    assert(0 <= k && k < len_sequence());
    //std::cout << root->cnt_current << " " << root->cnt_erased << " " << root->cnt_pop << '\n';
    if((root = erase_inner(root, k))) root->p = nullptr;
    if(tmp_node->mode == 1){
      // mode == 1の場合
      // pushされた要素が消されていなかった -> eraseして終了
      return;
    }else if(tmp_node->mode == 0){
      // この時点で全体のdiffminが-1になっている
      // これを0にするために[0, 最左のdiffmin)の残っている最小ノードを消す
      int r = leftmost_diffmin(root);
      node *v = find_min_current_prefix(root, r);
      assert(v->mode == 1);
      v->mode = 0;
      update(v);
      while(v->p){
        v = v->p;
        update(v);
      }
    }else{
      int l = (root->diffmin == 0 ? rightmost_diffmin(root) : 0);
      node *v = find_max_erased_suffix(root, l);
      v->mode = 1;
      update(v);
      while(v->p){
        v = v->p;
        update(v);
      }
    }
  }
  // 現在のmin 空の場合はinf
  Key min(){
    return root ? root->min_current->key : inf;
  }
  // 現在のkeyのsum
  KeySum sum(){
    return root ? root->sum_current : 0;
  }
};
template<typename Key, typename KeySum>
using prpq = partial_retroactive_priority_queue<Key, KeySum>;

#endif