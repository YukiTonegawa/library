/*
#include "../template.hpp"
#include "../string/string.hpp"

struct lazy_dynamic_bitvector{
private:
  using ull = unsigned long long;
  using rh = montgomery_rolling_hash;
  struct Lazy{
    bool __or, __and, __xor;
    Lazy(): __or(0), __and(1), __xor(0){}
    Lazy(bool a, bool b, bool c): __or(a), __and(b), __xor(c){}
    bool is_id(){return __or == 0 && __and == 1 && __xor == 0;}
    void reset(){__or = 0, __and = 1, __xor = 0;}
  };
  static constexpr bool or_and_xor(bool i, const Lazy &lz){
    return ((i | lz.__or) & lz.__and) ^ lz.__xor;
  }
  static constexpr void propagate_lazy(Lazy &l, const Lazy &r){
    if(!r.__or && r.__and){
      l.__xor ^= r.__xor;
      return;
    }
    l.__or  = ((l.__or & r.__or) ^ l.__xor) | r.__or;
    l.__and = ((l.__and ^ l.__xor) | r.__or) & r.__and;
    l.__xor = r.__xor;
  }
  struct node{
    node *l, *r;
    int sz, pop_sum;
    bool val, rev;
    Lazy lazy;
    ull hsum, hsum_rev;
    node(bool _val): l(nullptr), r(nullptr), sz(1),
    pop_sum(_val), val(_val), rev(false), lazy(), hsum(_val), hsum_rev(_val){}
  };
  node *root;
  int size(node *v){
    return !v ? 0 : v->sz;
  }
  void update(node *v){
    v->pop_sum = v->val + (v->l ? v->l->pop_sum : 0) + (v->r ? v->r->pop_sum : 0);
    v->hsum = v->hsum_rev = v->val;
    v->sz = 1;
    if(v->l){
      v->hsum = rh::add(v->l->hsum, rh::mul(v->hsum, rh::rpow[v->l->sz]));
      v->hsum_rev = rh::add(rh::mul(v->l->hsum_rev, rh::rpow[v->sz]), v->hsum_rev);
      v->sz += v->l->sz;
    }
    if(v->r){
      v->hsum = rh::add(v->hsum, rh::mul(v->r->hsum, rh::rpow[v->sz]));
      v->hsum_rev = rh::add(rh::mul(v->hsum_rev, rh::rpow[v->r->sz]), v->r->hsum_rev);
      v->sz += v->r->sz;
    }
  }
  void push_down(node *v){
    if(!v) return;
    if(!v->lazy.is_id()){
      if(v->l) propagate(v->l, v->lazy);
      if(v->r) propagate(v->r, v->lazy);
      v->lazy.reset();
    }
    if(v->rev){
      if(v->l) reverse(v->l);
      if(v->r) reverse(v->r);
      v->rev = false;
    }
  }
  void propagate(node *v, Lazy x){
    propagate_lazy(v->lazy, x);
    bool l = or_and_xor(0, x), r = or_and_xor(1, x);
    v->val = v->val ? r : l;
    ull all_pop = rh::mul(rh::sub(rh::rpow[v->sz], 1), rh::run_inv[1]);
    if(l && r){
      v->pop_sum = v->sz;
      v->hsum = v->hsum_rev = all_pop;
    }else if(l){
      v->pop_sum = v->sz - v->pop_sum;
      v->hsum = all_pop - v->hsum;
      v->hsum_rev = all_pop - v->hsum_rev;
    }else if(r){
      return;
    }else{
      v->pop_sum = 0;
      v->hsum = v->hsum_rev = 0;
    }
  }
  void reverse(node *v){
    std::swap(v->l, v->r);
    std::swap(v->hsum, v->hsum_rev);
    v->rev ^= 1;
  }
  // vの左の子をvの位置に持ってくる
  node *rotate_right(node *v){
    node *l = v->l;
    v->l = l->r;
    l->r = v;
    update(v);
    update(l);
    return l;
  }
  // vの右の子をvの位置に持ってくる
  node *rotate_left(node *v){
    node *r = v->r;
    v->r = r->l;
    r->l = v;
    update(v);
    update(r);
    return r;
  }
  // zig-zig, zig-zag
  node *splay_top_down(node *v, int k, node* &u){
    push_down(v);
    int szl = v->l ? v->l->sz : 0;
    if(k == szl) return u = v;
    if(k < szl){
      v->l = splay_top_down(v->l, k, u);
      update(v);
      if(v->l == u) return v;
      if(v->l->l == u) v = rotate_right(v);
      else v->l = rotate_left(v->l);
      return rotate_right(v);
    }else{
      v->r = splay_top_down(v->r, k - szl - 1, u);
      update(v);
      if(v->r == u) return v;
      if(v->r->r == u) v = rotate_left(v);
      else v->r = rotate_right(v->r);
      return rotate_left(v);
    }
    return v;
  }
  // zig
  node *splay_top_down(node *v, int k){
    node *u = nullptr;
    v = splay_top_down(v, k, u);
    if(v->l == u) return rotate_right(v);
    else if(v->r == u) return rotate_left(v);
    return v;
  }
  // 左がkノードになるように分割
  std::pair<node*, node*> split(node *v, int k){
    int n = size(v);
    if(k >= n) return {v, nullptr};
    v = splay_top_down(v, k);
    node *l = v->l;
    v->l = nullptr;
    update(v);
    return {l, v};
  }
  node *merge(node *l, node *r){
    if(!l || !r) return !l ? r : l;
    r = splay_top_down(r, 0);
    r->l = l;
    update(r);
    return r;
  }
  std::tuple<node*, node*, node*> split3(node *v, int l, int r){
    if(l == 0){
      auto [b, c] = split(v, r);
      return {nullptr, b, c};
    }
    v = splay_top_down(v, l - 1); //    (l - 1個)  /  v  / (残り)
    auto [b, c] = split(v->r, r - l); // cがnullptrまたはcの左が空
    v->r = nullptr; // vの右が空
    update(v);
    return {v, b, c};
  }
  // split3によって分割された組でないと壊れる
  node *merge3(node *a, node *b, node *c){
    node *v = merge(b, c); // O(1)
    if(!a) return v;
    a->r = v; // O(1)
    update(a);
    return a;
  }
  node *set_inner(node *v, int k, bool x){
    v = splay_top_down(v, k);
    v->val = x;
    update(v);
    return v;
  }
  node *get_inner(node *v, int k, bool &x){
    v = splay_top_down(v, k);
    x = v->val;
    return v;
  }
  node *update_inner(node *v, int l, int r, Lazy x){
    if(r == l) return v;
    auto [a, b, c] = split3(v, l, r);
    propagate(b, x);
    return merge3(a, b, c);
  }
  node *query_inner(node *v, int l, int r, int &res){
    if(r == l) return v;
    auto [a, b, c] = split3(v, l, r);
    res = b->pop_sum;
    return merge3(a, b, c);
  }
  node *reverse_inner(node *v, int l, int r){
    if(r == l) return v;
    auto [a, b, c] = split3(v, l, r);
    if(b) reverse(b);
    return merge3(a, b, c);
  }
  node *insert_inner(node *v, int k, node *u){
    if(k == size(v)){
      u->l = v;
      update(u);
      return u;
    }
    if(k == 0){
      u->r = v;
      update(u);
      return u;
    }
    v = splay_top_down(v, k);
    u->l = v->l;
    v->l = u;
    update(u);
    update(v);
    return v;
  }
  node *erase_inner(node *v, int k){
    v = splay_top_down(v, k);
    return merge(v->l, v->r);
  }
  node *build(const std::vector<bool> &v, int l, int r){
    int m = (l + r) >> 1;
    node *u = new node(v[m]);
    if(m > l) u->l = build(v, l, m);
    if(r > m + 1) u->r = build(v, m + 1, r);
    update(u);
    return u;
  }
  template<typename F>
  int bisect_from_left(node *v, F &f, int &ok, bool unpop){
    if(!v) return -1;
    push_down(v);
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    int m = ok + (unpop ? v->sz - v->pop_sum : v->pop_sum);
    if(!f(m)){
      ok = m;
      return -1;
    }
    int x = bisect_from_left(v->l, f, ok, unpop);
    if(x != -1) return x;
    ok += (unpop ? !v->val : v->val);
    if(f(ok)) return szl;
    int res = bisect_from_left(v->r, f, ok, unpop);
    return res == -1 ? res : res + szm;
  }
  template<typename F>
  int bisect_from_right(node *v, F &f, int &ok, bool unpop){
    if(!v) return -1;
    push_down(v);
    int szl = v->l ? v->l->sz : 0, szm = szl + 1;
    int m = ok + (unpop ? v->sz - v->pop_sum : v->pop_sum);
    if(!f(m)){
      ok = m;
      return -1;
    }
    int x = bisect_from_right(v->r, f, ok, unpop);
    if(x != -1) return x + szm;
    ok += (unpop ? !v->val : v->val);
    if(f(ok)) return szl;
    return bisect_from_right(v->l, f, ok, unpop);
  }
  // f(sum[l, r])がtrueになる最左のr. ない場合は-1
  template<typename F>
  int bisect_from_left(int l, F f, bool unpop){
    if(l >= size()) return -1;
    int ok =0;
    if(!l){
      int i = bisect_from_left(root, f, ok, unpop);
      if(i != -1) root = splay_top_down(root, i);
      return i;
    }
    root = splay_top_down(root, l - 1);
    int i = bisect_from_left(root->r, f, ok, unpop);
    if(i != -1){
      i += l;
      root = splay_top_down(root, i);
    }
    return i;
  }
  // f(sum[l, r])がtrueになる最右のl. ない場合は-1
  template<typename F>
  int bisect_from_right(int r, F f, bool unpop){
    if(r < 0) return -1;
    int ok = 0;
    if(r + 1 == size()){
      int i = bisect_from_right(root, f, ok, unpop);
      if(i != -1) root = splay_top_down(root, i);
      return i;
    }
    root = splay_top_down(root, r + 1);
    int i = bisect_from_right(root->l, f, ok, unpop);
    if(i != -1) root = splay_top_down(root, i);
    return i;
  }
  void to_string(node *v, int left_size, std::string &s){
    s[left_size + size(v->l)] = v->val ? '1' : '0';
    if(v->l) to_string(v->l, left_size, s);
    if(v->r) to_string(v->r, left_size + size(v->l) + 1, s);
  }
  node *range_hash_inner(node *v, int l, int r, ull &res){
    if(r == l) return v;
    auto [a, b, c] = split3(v, l, r);
    res = b->hsum;
    return merge3(a, b, c);
  }
  lazy_dynamic_bitvector(node *r): root(r){}
public:
  lazy_dynamic_bitvector(): root(nullptr){}
  lazy_dynamic_bitvector(const std::string &s, char one = '1'): root(nullptr){
    if(s.size() > 0){
      int n = s.size();
      std::vector<bool> tmp(n);
      for(int i = 0; i < n; i++) tmp[i] = (s[i] == one);
      root = build(tmp, 0, n);
    }
  }
  lazy_dynamic_bitvector(const std::vector<bool> &v): root(nullptr){
    if(!v.empty()) root = build(v, 0, v.size());
  }
  int size(){return size(root);}
  void set(int k, bool x){
    assert(0 <= k && k < size());
    root = set_inner(root, k, x);
  }
  bool get(int k){
    assert(0 <= k && k < size());
    bool res;
    root = get_inner(root, k, res);
    return res;
  }
  // [l, r)をxにする
  void range_set(int l, int r, bool x){
    assert(0 <= l && r <= size());
    if(x) root = update_inner(root, l, r, Lazy(1, 1, 0));
    else root = update_inner(root, l, r, Lazy(0, 0, 0));
  }
  // [l, r)のビットの0/1を反転
  void flip(int l, int r){
    assert(0 <= l && r <= size());
    root = update_inner(root, l, r, Lazy(0, 1, 1));
  }
  // [l, r)の位置を反転
  void reverse(int l, int r){
    assert(0 <= l && r <= size());
    root = reverse_inner(root, l, r);
  }
  void insert(int k, bool x){
    assert(0 <= k && k <= size());
    root = insert_inner(root, k, new node(x));
  }
  void erase(int k){
    assert(0 <= k && k < size());
    root = erase_inner(root, k);
  }
  void erase_range(int l, int r){
    assert(0 <= l && r <= size());
    auto [a, b, c] = split3(root, l, r);
    root = merge(a, c);
  }
  void cyclic_lshift(int l, int r, int k){
    assert(0 <= l && r <= size());
    int len = (r - l);
    if(!len || (k %= len) == 0) return;
    auto [a, b, c] = split3(root, l, r);
    auto [b1, b2] = split(b, k);
    root = merge3(a, merge(b2, b1), c);
  }
  void cyclic_rshift(int l, int r, int k){
    assert(0 <= l && r <= size());
    int len = (r - l);
    if(!len || (k %= len) == 0) return;
    cyclic_lshift(l, r, len - k);
  }
  int pop_count(int l, int r){
    assert(0 <= l && r <= size());
    int res = 0;
    root = query_inner(root, l, r, res);
    return res;
  }
  int rank1(int r){
    assert(r <= size());
    if(r == 0) return 0;
    auto [a, b] = split(root, r);
    int res = a->pop_sum;
    root = merge(a, b);
    return res;
  }
  int rank0(int r){
    return r - rank1(r);
  }
  int select1(int k){
    return bisect_from_left(0, [&](int pop){return pop > k;}, false);
  }
  int select0(int k){
    return bisect_from_left(0, [&](int pop){return pop > k;}, true);
  }
  int find_next1(int k){
    return bisect_from_left(k, [&](int pop){return pop > 0;}, false);
  }
  int find_prev1(int k){
    return bisect_from_right(k, [&](int pop){return pop > 0;}, false);
  }
  int find_next0(int k){
    return bisect_from_left(k, [&](int pop){return pop > 0;}, true);
  }
  int find_prev0(int k){
    return bisect_from_right(k, [&](int pop){return pop > 0;}, true);
  }
  // 01列に変換
  std::string to_string(){
    std::string s = "";
    to_string(root, 0, s);
    return s;
  }

  ull range_hash(int l, int r){
    assert(0 <= l && r <= size());
    ull res = 0;
    root = range_hash_inner(root, l, r, res);
    return res;
  }
  static bool compare_lexicographical(){

  }
  // -1: x < y
  //  0: x == y
  //  1: x > y
  static int compare_number(lazy_dynamic_bitvector &x, lazy_dynamic_bitvector &y){
    int xs = x.size(), ys = y.size();
    assert(xs && ys);
    if(xs > ys) return -1 * compare_number(y, x);
    // x.size() <= y.size()
    int diff_size = ys - xs;
    if(y.range_hash(0, diff_size)) return -1;
    // x[0, xs), y[diff_size, ys) (同じ長さ)で比較
    
  }
};
#include "../data_structure/bit_sequence/dynamic_bitset.hpp"
int main(){
  montgomery_rolling_hash::init(400);
  io_init();
  int n;
  std::cin >> n;
  string s;
  std::cin >> s;
  dynamic_bitset_rank_select b(s, 'p');
  dynamic_string ds(stovector(s));
  dynamic_string ans(stovector(s));
  int ans_i = 0, ans_j = 0;
  //std::cout << ans.to_vector() << '\n';

  range(i, 0, n){
    if(s[i] != 'p') continue;
    // 初めてpが現れる位置
    int j = b.select1(0);
    assert(j != -1 && j <= i);

    // [j, i]でreverse
    ds.reverse(j, i + 1);
    //std::cout << ds.to_vector() << '\n';

    if(ds < ans){
      ans.reverse(ans_j, ans_i + 1);
      ans.reverse(j, i + 1);
      ans_i = i, ans_j = j;
    }
  }
  auto v = ans.to_vector();
  //std::cout << v << '\n';
}
*/