#ifndef _STRING_HPP_
#define _STRING_HPP_

#include <numeric>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <array>
#include <cassert>
#include <queue>
#include "string_basic.hpp"
#include "../data_structure/range_query/sparse_table.hpp"

template<typename Val, Val (*id)()>
struct trie{
  struct node{
    Val val;
    node *p;
    std::array<node*, 26> c;
    node(Val _val = id()): val(_val), p(nullptr){c.fill(nullptr);}
  };
  node *root;
  trie(): root(new node()){}

  node *get_root(){return root;}
  // sを追加して終端ノードを返す
  node *insert(const std::string &s){
    node *v = root, *u;
    for(int i = 0; i < s.size(); i++){
      if(v->c[s[i] - 'a']){
        v = v->c[s[i] - 'a'];
      }else{
        u = (v->c[s[i] - 'a'] = new node());
        u->p = v;
        v = u;
      }
    }
    return v;
  }
  // sを検索, ない場合はnullptr
  node *find(const std::string &s){
    node *v = root;
    for(int i = 0; i < s.size(); i++){
      if(v->c[s[i] - 'a']) v = v->c[s[i] - 'a'];
      else return nullptr;
    }
    return v;
  }
  // 親がある場合書き換えてtrueを返す
  bool parent(node* &v){
    return v = v->p;
  }
  // 子がある場合書き換えてtrueを返す
  bool next(node* &v, char c){
    return v = v->c[c - 'a'];
  }
};
template<typename Val, Val (*id)()>
struct aho_corasick{
  struct node{
    bool is_end; // いずれかのパターンの終端
    Val val;
    node *p, *s;
    std::array<node*, 26> c;
    node(Val _val = id()): is_end(false), val(_val), p(nullptr), s(nullptr){c.fill(nullptr);}
  };
  node *root;
  // 文字cに対応するノードvの最長サフィックスリンクを貼る
  void find_longest_suffix(node *v, char c){
    // vの最長サフィックスリンクは [1]根 / [2]c1文字 / [3]pの最長サフィックスリンクを何回か辿ったノードの子　のいずれか
    // [2]のパターンは[3]に含まれる.
    // 最長のものが欲しいので, pの最長サフィックスリンクを何回か辿っていく過程でcに対応する子ノードを持つ最初のノードが答え
    // 存在しなければ根が答え

    // 文字列Sを追加する過程で最長サフィックスリンクを辿る回数はO(|S|)
    // i + 1文字目を追加するとき, 可能な接頭辞の長さは高々1増える.
    // 新たなノードの最長サフィックスリンクを決める過程で親のリンクを1 + k回辿ると, 可能な接頭辞の長さはk以上減る
    // よってリンクを辿る回数は2|S|回

    assert(v);
    node *p = v->p;
    // vまたはpが根なら根
    if(!p || p == root){
      v->s = root;
      return;
    }
    int k = c - 'a';
    while(p != root){//
      p = p->s;
      if(p->c[k]){
        v->s = p->c[k];
        return;
      }
    }
    // 見つからなかった(根にする)
    v->s = root;
  }
  aho_corasick(): root(new node()){}
  node *get_root(){return root;}
  // sを追加して終端ノードを返す
  // 動的に使うとサフィックスリンクが壊れる
  node *insert(const std::string &s){
    node *v = root, *u;
    for(int i = 0; i < s.size(); i++){
      if(v->c[s[i] - 'a']){
        v = v->c[s[i] - 'a'];
      }else{
        u = (v->c[s[i] - 'a'] = new node());
        u->p = v;
        v = u;
      }
    }
    v->is_end = true;
    return v;
  }
  // sを検索, ない場合はnullptr
  node *find(const std::string &s){
    node *v = root;
    for(int i = 0; i < s.size(); i++){
      if(v->c[s[i] - 'a']) v = v->c[s[i] - 'a'];
      else return nullptr;
    }
    return v;
  }
  // 親がある場合書き換えてtrueを返す
  bool parent(node* &v){
    return v = v->p;
  }
  // 子がある場合書き換えてtrueを返す
  bool next(node* &v, char c){
    return v = v->c[c - 'a'];
  }
  void build_suffix(){
    std::queue<node*> q;
    q.push(root);
    root->s = root;
    while(!q.empty()){
      node *v = q.front();
      q.pop();
      for(int i = 0; i < v->c.size(); i++){
        if(v->c[i]){
          find_longest_suffix(v->c[i], i + 'a');
          q.push(v->c[i]);
        }
      }
    }
  }
  // 現在のノード + cの最長サフィックスに書き換える
  // ない場合は根にしてfalse
  bool next_suffix(node* &v, char c){
    while(v != root && v->c[c - 'a'] == nullptr) v = v->s;
    if(v->c[c - 'a']) return v = v->c[c - 'a'];
    return false;
  }
  // 文字列SがパターンTを連続部分文字列として含む場合, あるiが存在して S_0, S_1...S_iの接尾辞の集合はTを含む
  // 現在のノードに対応する接尾辞の集合が含むパターンを列挙する
  // 実装上, 接尾辞に相当するノードを全て返して終端ノードかの判定は自分でやる
  std::vector<node*> find_suffix(node *v){
    std::vector<node*> res;
    while(v != root){
      res.push_back(v);
      v = v->s;
    }
    return res;
  }
  // 現在の接尾辞がいずれかのパターンとマッチしたか
  bool match_any(node *v){
    while(v != root){
      if(v->is_end) return true;
      v = v->s;
    }
    return false;
  }
};
//S[i]を中心とする最大回文半径
//偶数文字の回文を検出したい場合$を挟む必要がある
template<typename T>
std::vector<int> manacher(const std::vector<T> &s){
  int n = s.size(), c = 0;
  std::vector<int> rad(n);
  for(int r = 0; r < n; r++){
    int l = c - (r - c);
    if(r + rad[l] < c + rad[c]) rad[r] = rad[l];
    else{
      int rr = c + rad[c];
      int rl = r - (rr - r);
      while(rl >= 0 && rr < n && s[rl] == s[rr]) rr++, rl--;
      rad[r] = rr - r;
      c = r;
    }
  }
  return rad;
}
template<typename T>
std::vector<T> insert_dollar(const std::vector<T> &s){
  int n = s.size();
  std::vector<T> res(n * 2 + 1, '$');
  for(int i = 0; i < n; i++) res[i * 2 + 1] = s[i];
  return res;
}

std::vector<int> sa_naive(const std::vector<int>& s){
  int n = int(s.size());
  std::vector<int> sa(n);
  std::iota(sa.begin(), sa.end(), 0);
  std::sort(sa.begin(), sa.end(), [&](int l, int r){
    if (l == r) return false;
    while (l < n && r < n){
      if (s[l] != s[r]) return s[l] < s[r];
      l++;
      r++;
    }
    return l == n;
  });
  return sa;
}

std::vector<int> sa_doubling(const std::vector<int>& s) {
  int n = int(s.size());
  std::vector<int> sa(n), rnk = s, tmp(n);
  std::iota(sa.begin(), sa.end(), 0);
  for(int k = 1; k < n; k *= 2){
    auto cmp = [&](int x, int y){
      if (rnk[x] != rnk[y]) return rnk[x] < rnk[y];
      int rx = x + k < n ? rnk[x + k] : -1;
      int ry = y + k < n ? rnk[y + k] : -1;
      return rx < ry;
    };
    std::sort(sa.begin(), sa.end(), cmp);
    tmp[sa[0]] = 0;
    for (int i = 1; i < n; i++){
      tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
    }
    std::swap(tmp, rnk);
  }
  return sa;
}

// SA-IS, linear-time suffix array construction
// Reference:
// G. Nong, S. Zhang, and W. H. Chan,
// Two Efficient Algorithms for Linear Time Suffix Array Construction
template<int THRESHOLD_NAIVE = 10, int THRESHOLD_DOUBLING = 40>
std::vector<int> sa_is(const std::vector<int>& s, int upper){
  int n = int(s.size());
  if (n == 0) return {};
  if (n == 1) return {0};
  if (n == 2) {
    if (s[0] < s[1]) {
      return {0, 1};
    } else {
      return {1, 0};
    }
  }
  if (n < THRESHOLD_NAIVE) {
    return sa_naive(s);
  }
  if (n < THRESHOLD_DOUBLING) {
    return sa_doubling(s);
  }

  std::vector<int> sa(n);
  std::vector<bool> ls(n);
  for (int i = n - 2; i >= 0; i--) {
    ls[i] = (s[i] == s[i + 1]) ? ls[i + 1] : (s[i] < s[i + 1]);
  }
  std::vector<int> sum_l(upper + 1), sum_s(upper + 1);
  for (int i = 0; i < n; i++) {
    if (!ls[i]) {
      sum_s[s[i]]++;
    } else {
      sum_l[s[i] + 1]++;
    }
  }
  for (int i = 0; i <= upper; i++) {
    sum_s[i] += sum_l[i];
    if (i < upper) sum_l[i + 1] += sum_s[i];
  }

  auto induce = [&](const std::vector<int>& lms) {
    std::fill(sa.begin(), sa.end(), -1);
    std::vector<int> buf(upper + 1);
    std::copy(sum_s.begin(), sum_s.end(), buf.begin());
    for (auto d : lms) {
      if (d == n) continue;
      sa[buf[s[d]]++] = d;
    }
    std::copy(sum_l.begin(), sum_l.end(), buf.begin());
    sa[buf[s[n - 1]]++] = n - 1;
    for (int i = 0; i < n; i++) {
      int v = sa[i];
      if (v >= 1 && !ls[v - 1]) {
        sa[buf[s[v - 1]]++] = v - 1;
      }
    }
    std::copy(sum_l.begin(), sum_l.end(), buf.begin());
    for (int i = n - 1; i >= 0; i--) {
      int v = sa[i];
      if (v >= 1 && ls[v - 1]) {
        sa[--buf[s[v - 1] + 1]] = v - 1;
      }
    }
  };

  std::vector<int> lms_map(n + 1, -1);
  int m = 0;
  for (int i = 1; i < n; i++) {
    if (!ls[i - 1] && ls[i]) {
      lms_map[i] = m++;
    }
  }
  std::vector<int> lms;
  lms.reserve(m);
  for (int i = 1; i < n; i++) {
    if (!ls[i - 1] && ls[i]) {
      lms.push_back(i);
    }
  }

  induce(lms);

  if (m) {
    std::vector<int> sorted_lms;
    sorted_lms.reserve(m);
    for (int v : sa) {
      if (lms_map[v] != -1) sorted_lms.push_back(v);
    }
    std::vector<int> rec_s(m);
    int rec_upper = 0;
    rec_s[lms_map[sorted_lms[0]]] = 0;
    for (int i = 1; i < m; i++) {
      int l = sorted_lms[i - 1], r = sorted_lms[i];
      int end_l = (lms_map[l] + 1 < m) ? lms[lms_map[l] + 1] : n;
      int end_r = (lms_map[r] + 1 < m) ? lms[lms_map[r] + 1] : n;
      bool same = true;
      if (end_l - l != end_r - r) {
        same = false;
      } else {
        while (l < end_l) {
          if (s[l] != s[r]) {
            break;
          }
          l++;
          r++;
        }
        if (l == n || s[l] != s[r]) same = false;
      }
      if (!same) rec_upper++;
      rec_s[lms_map[sorted_lms[i]]] = rec_upper;
    }

    auto rec_sa = sa_is<THRESHOLD_NAIVE, THRESHOLD_DOUBLING>(rec_s, rec_upper);

    for (int i = 0; i < m; i++) sorted_lms[i] = lms[rec_sa[i]];
    induce(sorted_lms);
  }
  return sa;
}
std::vector<int> suffix_array(const std::vector<int>& s, int upper) {
  assert(0 <= upper);
  for (int d : s) assert(0 <= d && d <= upper);
  auto sa = sa_is(s, upper);
  return sa;
}

template <class T>
std::vector<int> suffix_array(const std::vector<T>& s) {
  int n = int(s.size());
  std::vector<int> idx(n);
  iota(idx.begin(), idx.end(), 0);
  sort(idx.begin(), idx.end(), [&](int l, int r) { return s[l] < s[r]; });
  std::vector<int> s2(n);
  int now = 0;
  for (int i = 0; i < n; i++) {
    if (i && s[idx[i - 1]] != s[idx[i]]) now++;
    s2[idx[i]] = now;
  }
  return sa_is(s2, now);
}

std::vector<int> suffix_array(const std::string& s) {
  int n = int(s.size());
  std::vector<int> s2(n);
  for (int i = 0; i < n; i++) s2[i] = s[i];
  return sa_is(s2, 255);
}

// lcp[i] := substr(sa[i] ... n)とsubstr(sa[i + 1]...n)の接頭辞の共通する文字の数
// lcp[n - 1]は0
template <class T>
std::vector<int> lcp_array(const std::vector<T>& s, const std::vector<int>& sa) {
  int n = int(s.size());
  assert(n >= 1);
  std::vector<int> rnk(n);
  for (int i = 0; i < n; i++) rnk[sa[i]] = i;
  std::vector<int> lcp(n - 1);
  int h = 0;
  for (int i = 0; i < n; i++) {
    if (h > 0) h--;
    if (rnk[i] == 0) continue;
    int j = sa[rnk[i] - 1];
    for (; j + h < n && i + h < n; h++) if (s[j + h] != s[i + h]) break;
    lcp[rnk[i] - 1] = h;
  }
  lcp.push_back(0);
  return lcp;
}
std::vector<int> lcp_array(const std::string& s, const std::vector<int>& sa) {
  int n = int(s.size());
  std::vector<int> s2(n);
  for (int i = 0; i < n; i++) s2[i] = s[i];
  return lcp_array(s2, sa);
}
// substr(i...n)が辞書順で何番目か
std::vector<int> rank_array(const std::vector<int>& sa){
  int n = (int)sa.size();
  std::vector<int> rank(n);
  for(int i = 0; i < n; i++) rank[sa[i]] = i;
  return rank;
}
void sa_debug(const std::string& s){
  std::vector<int> sa = suffix_array(s);
  std::vector<int> lcp = lcp_array(s, sa);
  int n = (int)s.size();
  for(int i = 0; i < n; i++) std::cout << s.substr(sa[i]) << " " << lcp[i] << '\n';
}

struct lcp_arbitrary_pair{
  static constexpr int _min(int a, int b){return std::min(a, b);}
  std::vector<int> sa;
  std::vector<int> rank;
  rmq<int> st;

  template<typename T>
  lcp_arbitrary_pair(const std::vector<T> &s): sa(suffix_array(s)), rank(rank_array(sa)), st(lcp_array(s, sa)){}
  lcp_arbitrary_pair(const std::string &s): sa(suffix_array(s)), rank(rank_array(sa)), st(lcp_array(s, sa)){}
  lcp_arbitrary_pair(const std::vector<int> _sa, const std::vector<int> _lcp, const std::vector<int> _rank): sa(_sa), rank(_rank), st(_lcp){}
  // v[i...n)とv[j...n)のlcp
  int lcp(int i, int j){
    if(i == j) return rank.size() - i;
    i = rank[i], j = rank[j];
    if(i > j) std::swap(i, j);
    return st.query(i, j);
  }
  // v[i...n) <= v[j...n)か
  // i == jのときのみ等号が成り立つ
  bool compare(int i, int j){
    if(i == j) return true;
    int len = lcp(i, j);
    int n = rank.size();
    if(i + len == n) return true;
    if(j + len == n) return false;
    return rank[i + len] <= rank[j + len];
  }
};

// z[i] := lcp(substr(0...n), substr(i...n))
template<typename T>
std::vector<int> z_algorithm(const std::vector<T>& s){
  int n = int(s.size());
  if (n == 0) return {};
  std::vector<int> z(n);
  z[0] = 0;
  for (int i = 1, j = 0; i < n; i++) {
    int& k = z[i];
    k = (j + z[j] <= i) ? 0 : std::min(j + z[j] - i, z[i - j]);
    while (i + k < n && s[k] == s[i + k]) k++;
    if (j + z[j] < i + z[i]) j = i;
  }
  z[0] = n;
  return z;
}

std::vector<int> z_algorithm(const std::string& s) {
  int n = int(s.size());
  std::vector<int> s2(n);
  for (int i = 0; i < n; i++) {
    s2[i] = s[i];
  }
  return z_algorithm(s2);
}

// res[i] := lcp(s[i以降], pat)
std::vector<int> find_lcp(const std::string &s, const std::string &pat){
  std::string s2 = pat + '$' + s;
  std::vector<int> res = z_algorithm(s2);
  res.erase(res.begin(), res.begin() + (int)pat.size() + 1);
  return res;
}

// res[i] := s_0, s_1....s_k と s_i-k ...s_i-1, s_iが一致するような最大のk
// ただし, k != i
// res[0] = -1
std::vector<int> mp_algorithm(const std::string &s){
  int n = s.size();
  std::vector<int> res(n + 1);
  res[0] = -1;
  int j = -1;
  for(int i = 0; i < n; i++){
    while(j >= 0 && s[i] != s[j]) j = res[j];
    res[i + 1] = ++j;
  }
  return res;
}
template<typename T>
std::vector<int> mp_algorithm(const std::vector<T> &s){
  int n = s.size();
  std::vector<int> res(n + 1);
  res[0] = -1;
  int j = -1;
  for(int i = 0; i < n; i++){
    while(j >= 0 && s[i] != s[j]) j = res[j];
    res[i + 1] = ++j;
  }
  return res;
}

// パターンの検索テーブルを作る O(|pat|)
// Sの中のS_i, S_i+1...がパターンと一致するiを列挙 O(|S|)
template<typename T>
struct kmp_algorithm{
private:
  std::vector<T> pat;
  std::vector<int> t;
public:
  kmp_algorithm(const std::vector<T> &pat){
    init(pat);
  }
  void init(const std::vector<T> &_pat){
    pat = _pat;
    int n = pat.size();
    assert(n);
    t.resize(n + 1);
    t[0] = -1;
    int j = -1;
    for(int i = 0; i < n;){
      while(j >= 0 && pat[i] != pat[j]) j = t[j];
      i++, j++;
      if(i < n && pat[i] == pat[j]) t[i] = t[j];
      else t[i] = j;
    }
  }
  // 一致する場所を返す
  std::vector<int> find(const std::vector<T> &s){
    int n = s.size();
    std::vector<int> res;
    for(int i = 0, j = 0; i < n;){
      while(j >= 0 && s[i] != pat[j]) j = t[j];
      i++, j++;
      if(j == (int)pat.size()){
        res.push_back(i - j);
        j = t[j];
      }
    }
    return res;
  }
};
#endif
