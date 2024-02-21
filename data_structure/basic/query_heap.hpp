#ifndef _QUERY_HEAP_H_
#define _QUERY_HEAP_H_
#include "../../algebraic_structure/abelian_group.hpp"
#include "interval_heap.hpp"
#include <queue>
#include <vector>

// キーが一定(定数)以下または下位k項(定数)の総和/サイズなどをO(1)で求められる
// split_policyには{サイズ, lowの最大のキー, その値, 一時的な下位の値の総和}が与えられる
// push, popした後にsplit_policyがtrueになるまで要素が移される
// 同じキーの要素が複数挿入された場合, 最初に入れた方を小さいとみなす

// 可能なヒープの分割の例
// {キーがx未満, x以上}
// {サイズがk未満, それ以降}
// {要素の総和がx未満, x以上(和が単調な場合)}

// 取得可能な情報の例
// 下位/上位のサイズ, 最大/最小のキー, 総和
template<typename Key, typename abelian_group>
struct query_heap{
  using Val = typename abelian_group::Val;
  static constexpr auto id    = abelian_group::id;
  static constexpr auto inv   = abelian_group::inv;
  static constexpr auto merge = abelian_group::merge;
protected:
  using p = std::pair<Key, Val>;
  using func = std::function<bool(const int, const Key&, const Val&, const Val&)>;
  func split_policy;
  interval_heap_map<Key, Val> low, hi;
  Val low_sum, hi_sum;
public:
  query_heap(){}
  query_heap(func f): split_policy(f), low_sum(id()), hi_sum(id()){}
  query_heap(const std::vector<p> &_v, func f): split_policy(f), low_sum(id()), hi_sum(id()){
    int n = _v.size();
    for(int i = 0; i < n; i++) push(_v[i]);
  }
  void push(p x){
    if(!hi.empty() && hi.min().first <= x.first){
      hi.push(x);
      hi_sum = merge(hi_sum, x.second);
      return;
    }
    low.push(x);
    low_sum = merge(low_sum, x.second);
    while(!low.empty() && !split_policy(low.size(), low.max().first, low.max().second, low_sum)){
      p y = low.pop_max();
      low_sum = merge(low_sum, inv(y.second));
      hi.push(y);
      hi_sum = merge(hi_sum, y.second);
    }
  }
  int size(){
    return low.size() + hi.size();
  }
  int low_size(){
    return low.size();
  }
  int hi_size(){
    return hi.size();
  }
  int empty(){
    return low.empty() && hi.empty();
  }
  p min(){
    if(low.empty()) return hi.min();
    return low.min();
  }
  p max(){
    if(hi.empty()) return low.max();
    return hi.max();
  }
  p low_max(){
    return low.max();
  }
  p hi_min(){
    return hi.min();
  }
  p pop_min(){
    if(low.empty()){
      p y = hi.pop_min();
      hi_sum = hi.empty() ? id() : merge(hi_sum, inv(y.second));
      return y;
    }
    p x = low.pop_min();
    low_sum = (low.empty() ? id() : merge(low_sum, inv(x.second)));
    while(!hi.empty()){
      p y = hi.min();
      if(split_policy(low.size() + 1, y.first, y.second, merge(low_sum, y.second))){
        hi.pop_min();
        hi_sum = hi.empty() ? id() : merge(hi_sum, inv(y.second));
        low.push(y);
        low_sum = merge(low_sum, y.second);
      }else break;
    }
    return x;
  }
  p pop_max(){
    if(!hi.empty()){
      p y = hi.pop_max();
      hi_sum = hi.empty() ? id() : merge(hi_sum, inv(y.second));
      return y;
    }
    p x = low.pop_max();
    low_sum = (low.empty() ? id() : merge(low_sum, inv(x.second)));
    return x;
  }
  Val low_query(){
    return low_sum;
  }
  Val hi_query(){
    return hi_sum;
  }
};

// 削除不能版, 1.3倍くらい早い
template<typename Key, typename abelian_group>
struct query_heap_unpoppable{
  using Val = typename abelian_group::Val;
  static constexpr auto id    = abelian_group::id;
  static constexpr auto inv   = abelian_group::inv;
  static constexpr auto merge = abelian_group::merge;
protected:
  struct p{
    Key key;
    Val val;
    p(const std::pair<Key, Val> &x): key(x.first), val(x.second){}
    bool operator < (const p& r)const{
      return key < r.key;
    }
    bool operator > (const p& r)const{
      return key > r.key;
    }
  };
  using func = std::function<bool(const int, const Key&, const Val&, const Val&)>;
  func split_policy;
  std::priority_queue<p, std::vector<p>, std::less<p>> low;
  std::priority_queue<p, std::vector<p>, std::greater<p>> hi;
  Val low_sum, hi_sum;
public:
  query_heap_unpoppable(){}
  query_heap_unpoppable(func f): split_policy(f), low_sum(id()), hi_sum(id()){}
  query_heap_unpoppable(const std::vector<p> &_v, func f): split_policy(f), low_sum(id()), hi_sum(id()){
    int n = _v.size();
    for(int i = 0; i < n; i++) push(_v[i]);
  }
  void push(std::pair<Key, Val> x){
    if(!hi.empty() && hi.top().key <= x.first){
      hi.push(x);
      hi_sum = merge(hi_sum, x.second);
      return;
    }
    low.push(x);
    low_sum = merge(low_sum, x.second);
    while(!low.empty() && !split_policy(low.size(), low.top().key, low.top().val, low_sum)){
      p y = low.top();
      low.pop();
      low_sum = merge(low_sum, inv(y.val));
      hi.push(y);
      hi_sum = merge(hi_sum, y.val);
    }
  }
  int size(){
    return low.size() + hi.size();
  }
  int low_size(){
    return low.size();
  }
  int hi_size(){
    return hi.size();
  }
  int empty(){
    return low.empty() && hi.empty();
  }
  p low_max(){
    return low.top();
  }
  p hi_min(){
    return hi.top();
  }
  Val low_query(){
    return low_sum;
  }
  Val hi_query(){
    return hi_sum;
  }
};



// 下位k要素をマージ, 同じキーの要素が複数ある場合先に入れた方を優先
template<typename Key, typename abelian_group>
struct query_heap_k_elements : query_heap<Key, abelian_group>{
  using Val = typename abelian_group::Val;
  static constexpr auto id = abelian_group::id;
  query_heap_k_elements(int k){
    this->low_sum = id();
    this->hi_sum = id();
    this->split_policy = [k](const int a, const Key& b, const Val& c, const Val& d){
      return a <= k;
    };
  }
  query_heap_k_elements(const std::vector<std::pair<Key, Val>> &_v, int k){
    this->low_sum = id();
    this->hi_sum = id();
    this->split_policy = [k](const int a, const Key& b, const Val& c, const Val& d){
      return a <= k;
    };
    int n = _v.size();
    for(int i = 0; i < n; i++) push(_v[i]);
  }
};

// 下位k要素をマージ, 同じキーの要素が複数ある場合先に入れた方を優先
template<typename Key, typename abelian_group>
struct query_heap_k_elements_unpoppable : query_heap_unpoppable<Key, abelian_group>{
  using Val = typename abelian_group::Val;
  static constexpr auto id = abelian_group::id;
  query_heap_k_elements_unpoppable(int k){
    this->low_sum = id();
    this->hi_sum = id();
    this->split_policy = [k](const int a, const Key& b, const Val& c, const Val& d){
      return a <= k;
    };
  }
  query_heap_k_elements_unpoppable(const std::vector<std::pair<Key, Val>> &_v, int k){
    this->low_sum = id();
    this->hi_sum = id();
    this->split_policy = [k](const int a, const Key& b, const Val& c, const Val& d){
      return a <= k;
    };
    int n = _v.size();
    for(int i = 0; i < n; i++) push(_v[i]);
  }
};

// キーを昇順に並べた時に値の和がx未満になるように分割, 同じキーの要素が複数ある場合先に入れた方を優先
template<typename Key, typename abelian_group>
struct query_heap_lower_bound_sum : query_heap<Key, abelian_group>{
  using Val = typename abelian_group::Val;
  static constexpr auto id = abelian_group::id;
  query_heap_lower_bound_sum(Val x){
    this->low_sum = id();
    this->hi_sum = id();
    this->split_policy = [x](const int a, const Key& b, const Val& c, const Val& d){
      return d < x;
    };
  }
  query_heap_lower_bound_sum(const std::vector<std::pair<Key, Val>> &_v, Val x){
    this->low_sum = id();
    this->hi_sum = id();
    this->split_policy = [x](const int a, const Key& b, const Val& c, const Val& d){
      return d < x;
    };
    int n = _v.size();
    for(int i = 0; i < n; i++) push(_v[i]);
  }
};
// キーを昇順に並べた時に値の和がx未満になるように分割, 同じキーの要素が複数ある場合先に入れた方を優先
template<typename Key, typename abelian_group>
struct query_heap_lower_bound_sum_unpoppable : query_heap_unpoppable<Key, abelian_group>{
  using Val = typename abelian_group::Val;
  static constexpr auto id = abelian_group::id;
  query_heap_lower_bound_sum_unpoppable(Val x){
    this->low_sum = id();
    this->hi_sum = id();
    this->split_policy = [x](const int a, const Key& b, const Val& c, const Val& d){
      return d < x;
    };
  }
  query_heap_lower_bound_sum_unpoppable(const std::vector<std::pair<Key, Val>> &_v, Val x){
    this->low_sum = id();
    this->hi_sum = id();
    this->split_policy = [x](const int a, const Key& b, const Val& c, const Val& d){
      return d < x;
    };
    int n = _v.size();
    for(int i = 0; i < n; i++) push(_v[i]);
  }
};
#endif