#ifndef _HUFFMAN_TREE_H_
#define _HUFFMAN_TREE_H_
#include <vector>
#include <string>
#include <algorithm>
#include "../data_structure/basic/binary_heap.hpp"

// n - 1個の追加ノードを含むグラフを構築
// 根は2 * n - 2, [0, n)は葉
// ans[i] := {左ノード, 右ノード} 葉の場合{-1, -1}
// ∑{0 <= i < n} v[i] * depth[i]が最小になる

// 全ての要素が1以上でハフマン木を構築したとき, 高さkのノードの値はフィボナッチ数列のk項目以上
// -> 高さがO(log(N * 値の最大値))
template<typename Val>
std::vector<std::pair<int, int>> huffman_tree(const std::vector<Val> &v){
  int n = v.size();
  assert(n);
  std::vector<std::pair<Val, int>> tmp(n);
  for(int i = 0; i < n; i++) tmp[i] = {v[i], i};
  binary_heap_map<Val, int> h(tmp);
  std::vector<std::pair<int, int>> ans(2 * n - 1, {-1, -1});
  while(h.size() >= 2){
    auto [x, xi] = h.pop_min();
    auto [y, yi] = h.pop_min();
    ans[n] = {yi, xi};
    h.push({x + y, n++});
  }
  return ans;
}
// ハフマン符号, ∑長さ×重み　が最小
template<typename Val>
std::vector<std::string> huffman_code(const std::vector<Val> &v){
  int n = v.size();
  assert(n);
  if(n == 1) return {"0"};
  auto t = huffman_tree<Val>(v);
  std::vector<std::string> ans(2 * n - 1, "");
  for(int i = 2 * n - 2; i >= n; i--){
    ans[t[i].first] = ans[i] + '0';
    ans[t[i].second] = ans[i] + '1';
  }
  ans.resize(n);
  return ans;
}

/*
template<typename Val>
struct dynamic_huffman_tree{
private:
  Val swd;
  int n, r;
  std::vector<Val> v;
  struct node_info{
    int p, l, r;
    node_info(): p(-1), l(-1), r(-1){}
    node_info(int p, int l, int r): p(p), l(l), r(r){}
  };
  std::vector<node_info> u; // {親, 左, 右}
  //std::unordered_map<Val, std::unordered_set<int>> val_idx; // val_idx[x] := 値がxの要素のインデックス
  //std::unordered_map<Val, int> max_pair; // max_pair[x] := 値がxのペアとして最大のもののインデックス
  //std::unordered_map<Val, int> min_pair; // min_pair[x] := 値がxのペアとして最小のもののインデックス

  bool __is_root(int i){
    return i == r;
  }
  void __modify(int i){
    if(v[u[i].l] < v[u[i].r]) std::swap(u[i].l, u[i].r);
  }
  // i, jをswap
  void __swap_elem(int i, int j){
    int pi = u[i].p, pj = u[j].p;
    if(pi == pj){
      __modify(pi);
      return;
    }
    Val x, y;
    if(u[pi].l == i){
      u[pi].l = j;
      x = u[pi].r;
    }else{
      u[pi].r = j;
      x = u[pi].l;
    }
    if(u[pj].l == j){
      u[pj].l = i;
      y = u[pj].r;
    }else{
      u[pj].r = i;
      y = u[pj].l;
    }
    std::swap(u[i].p, u[j].p);
    __modify(pi);
    __modify(pj);

    // {i, x}, {j, y}を消して{i, y}, {j, x}を追加

  }
public:
  dynamic_huffman_tree(const std::vector<Val> &v){

  }
  int root(){
    return r;
  }
  // 末尾に1を追加
  void push_back(){

  }
  // v[i]++
  void add(int i){
    Val x = v[i]++;
    swd++;
    if(__is_root(i)) return;

    auto itr = max_pair.find(x);
    int k = *itr;
    if(v[k] == x){

    }
    
    //__swap_elem(i, j);
    add(plr[i].p);
  }
  // v[i]--
  void del(int i){
    Val x = v[i]--;
    assert(v[i]); // 0になるとハフマン符号が壊れる
    swd--;
    if(__is_root(i)) return;
    //auto []
    int j = -1;
    __swap_elem(i, j);
    del(plr[i].p);
  }
  // v[i]のハフマン符号, O(log(S))
  std::string encode(int i){
    std::string ret = "";
    while(true){
      int p = u[i].p;
      if(p == -1) break;
      ret += (i == u[p].l ? '0' : '1');
      i = p;
    }
    std::reverse(ret);
    return ret;
  }
  // ∑ depth[i] * v[i]
  Val sum_weighted_depth(){
    return swd;
  }
};
*/
#endif