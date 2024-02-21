#ifndef _PERSISTENT_SET_H_
#define _PERSISTENT_SET_H_
#include <vector>
#include <cassert>
#include <numeric>
#include "../../minior/persistent_set_inner.hpp"

template<typename Key>
struct persistent_set_iter{
  static constexpr Key inf = persistent_set<Key>::inf;
private:
  using pst = persistent_set<Key>;
  using node = typename persistent_set<Key>::node;
  using iter = persistent_set_iter<Key>;
  node *root;
  persistent_set_iter(node *v): root(v){}
public:
  persistent_set_iter(): root(nullptr){}
  // ソート済みでないと壊れる
  persistent_set_iter(const std::vector<Key> &v): root(pst::build_sorted(v)){}
  int size(){
    return pst::size(root);
  }
  bool empty(){
    return pst::empty(root);
  }
  iter insert(Key x){
    return iter(pst::insert(root, x));
  }
  iter erase(Key x){
    return iter(pst::erase(root, x));
  }
  Key min(){
    return pst::min(root);
  }
  Key max(){
    return pst::max(root);
  }
  bool find(Key x){
    return pst::find(root, x);
  }
  // x未満の要素数
  int low_count(Key x){
    return pst::low_count(root, x);
  }
  // x以上の最小, ない場合はinf
  Key lower_bound(Key x){
    return pst::lower_bound(root, x);
  }
  // x以下の最大, ない場合はinf
  Key lower_bound_rev(Key x){
    return pst::lower_bound_rev(root, x);
  }
  // k番目に小さい, ない場合はinf
  Key kth_smallest(int k){
    return pst::kth_smallest(root, k);
  }
};

template<typename Key>
struct persistent_multiset_iter{
  static constexpr Key inf = persistent_multiset<Key>::inf;
private:
  using pst = persistent_multiset<Key>;
  using node = typename persistent_multiset<Key>::node;
  using iter = persistent_multiset_iter<Key>;
  node *root;
  persistent_multiset_iter(node *v): root(v){}
public:
  persistent_multiset_iter(): root(nullptr){}
  // ソート済みでないと壊れる
  persistent_multiset_iter(const std::vector<Key> &v): root(pst::build_sorted(v)){}
  // 要素の総和
  int size(){
    return pst::size(root);
  }
  // 要素の種類数
  int size_unique(){
    return pst::size_unique(root);
  }
  bool empty(){
    return pst::empty(root);
  }
  iter insert(Key x, int cnt = 1){
    return iter(pst::insert(root, x, cnt));
  }
  iter erase(Key x, int cnt = 1){
    return iter(pst::erase(root, x, cnt));
  }
  Key min(){
    return pst::min(root);
  }
  Key max(){
    return pst::max(root);
  }
  int find(Key x){
    return pst::find(root, x);
  }
  // x未満の要素数
  int low_count_sum(Key x){
    return pst::low_count_sum(root, x);
  }
  // x未満の要素の種類数
  int low_count_unique(Key x){
    return pst::low_count_unique(root, x);
  }
  // x以上の最小, ない場合はinf
  Key lower_bound(Key x){
    return pst::lower_bound(root, x);
  }
  // x以下の最大, ない場合はinf
  Key lower_bound_rev(Key x){
    return pst::lower_bound_rev(root, x);
  }
  // (同じ要素を重複して数えて)k番目に小さい, ない場合はinf
  Key kth_smallest_sum(int k){
    return pst::kth_smallest_sum(root, k);
  }
  // (同じ要素を重複して数えずに)k番目に小さい, ない場合はinf
  Key kth_smallest_unique(int k){
    return pst::kth_smallest_unique(root, k);
  }
};


template<typename Key, typename Val>
struct persistent_map_iter{
  static constexpr Key inf = persistent_map<Key, Val>::inf;
  static constexpr Val inf_val = persistent_map<Key, Val>::inf_val;
private:
  using pmp = persistent_map<Key, Val>;
  using node = typename persistent_map<Key, Val>::node;
  using iter = persistent_map_iter<Key, Val>;
  node *root;
  persistent_map_iter(node *v): root(v){}
public:
  persistent_map_iter(): root(nullptr){}
  // ソート済みでないと壊れる
  // キーがユニークでない場合前にあるものを優先
  persistent_map_iter(const std::vector<std::pair<Key, Val>> &v): root(pmp::build_sorted_unsafe(v)){}
  int size(){
    return pmp::size(root);
  }
  bool empty(){
    return pmp::empty(root);
  }
  // ない場合はinf_val
  Val at(Key x){
    return pmp::at(root, x);
  }
  // key=xがある場合何もしない
  iter emplace(Key x, Val val){
    return iter(pmp::emplace(root, x, val));
  }
  // key=xがある場合置き換える
  iter emplace_replace(Key x, Val val){
    return iter(pmp::emplace_replace(root, x, val));
  }
  iter erase(Key x){
    return iter(pmp::erase(root, x));
  }
  std::pair<Key, Val> min(){
    return pmp::min(root);
  }
  std::pair<Key, Val> max(){
    return pmp::max(root);
  }
  bool find(Key x){
    return pmp::find(root, x);
  }
  // x未満の要素数
  int low_count(Key x){
    return pmp::low_count(root, x);
  }
  // x以上の最小, ない場合はinf
  std::pair<Key, Val> lower_bound(Key x){
    return pmp::lower_bound(root, x);
  }
  // x以下の最大, ない場合はinf
  std::pair<Key, Val> lower_bound_rev(Key x){
    return pmp::lower_bound_rev(root, x);
  }
  // k番目に小さい, ない場合はinf
  std::pair<Key, Val> kth_smallest(int k){
    return pmp::kth_smallest(root, k);
  }
};
#endif