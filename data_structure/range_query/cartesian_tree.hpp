#ifndef _CARTESIAN_TREE_H_
#define _CARTESIAN_TREE_H_
#include <vector>
#include <cassert>
#include <numeric>

// cmp = std::less          -> 最も小さい中で最左
// cmp = std::less_equal    -> 最も小さい中で最右
// cmp = std::greater       -> 最も大きい中で最左
// cmp = std::greater_equal -> 最も大きい中で最右
template<typename T, typename cmp>
struct cartesian_tree{
private:
  int n;
  std::vector<T> val;
  std::vector<int> l, r, p, same;
  int find_same(int i){
    i = same[i];
    if(p[i] == -1 || val[p[i]] != val[i]) return i;
    return same[i] = find_same(p[i]);
  }
  void build(const std::vector<T> &v){
    val = v;
    n = v.size();
    l.resize(n);
    r.resize(n, n);
    p.resize(n, -1);
    same.resize(n);
    std::iota(same.begin(), same.end(), 0);
    std::vector<int> st;
    for(int i = 0; i < n; i++){
      int last = -1;
      while(!st.empty() && cmp()(v[i], v[st.back()])){
        last = st.back();
        r[st.back()] = i;
        st.pop_back();
      }
      if(last != -1) p[last] = i;
      if(!st.empty()) p[i] = st.back();
      l[i] = st.empty() ? 0 : st.back() + 1;
      st.push_back(i);
    }
  }
public:
  cartesian_tree(const std::vector<T> &v){build(v);}
  // iの親のインデックス
  int parent(int i){
    return p[i];
  }
  // 部分木サイズ
  int size_subtree(int i){
    return r[i] - l[i];
  }
  // v[i]の部分木を表す半開区間[l, r)
  std::pair<int, int> segment_subtree(int i){
    return {l[i], r[i]};
  }
  // v[i]が最小であるような最大の半開区間[l, r)
  // 値がユニークなら range_subtree = range_minimum
  std::pair<int, int> segment_minumum(int i){
    return segment_subtree(find_same(i));
  }
};
#endif