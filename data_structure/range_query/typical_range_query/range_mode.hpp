#ifndef _RANGE_MODE_H_
#define _RANGE_MODE_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include "../../../traits.hpp"
#include "../../basic/hash_table.hpp"

template<typename T, is_unsigned_intle64<T>* = nullptr>
struct range_mode_query{
  int n, sq;
  uhash_map<T, int> cmp;
  std::vector<T> rev; // rev[i] := 座圧後にiである元の数
  std::vector<int> val;
  std::vector<std::vector<int>> pos; // pos[x] := 値xの位置の集合
  std::vector<int> rank; // rank[i] := 位置iの数xが[0, i)にいくつあるか
  std::vector<std::vector<std::pair<int, int>>> large; // large[i][j] := [i*sq, j*sq)の{mode, freq}
  range_mode_query(const std::vector<T> &v): n(v.size()), sq(sqrt(n)), cmp(v.size() * 2), val(n), pos(n), rank(n), large(n / sq + 3, std::vector<std::pair<int, int>>(n / sq + 3)){
    for(int i = 0, j = 0; i < n; i++){
      auto [f, x] = cmp.at(v[i]);
      if(f) val[i] = x;
      else{
        cmp.emplace(v[i], j);
        rev.push_back(v[i]);
        val[i] = j++;
      }
      rank[i] = pos[val[i]].size();
      pos[val[i]].push_back(i);
    }
    for(int i = 0; i * sq < n; i++){
      int next = std::min(n, i * sq + sq);
      std::vector<int> cnt(n, 0);
      int mode = 0;
      for(int j = i * sq, k = i + 1; j < n; j++){
        cnt[val[j]]++;
        if(cnt[val[j]] > cnt[mode]) mode = val[j];
        if(j + 1 == next){
          large[i][k++] = {mode, cnt[mode]};
          next = std::min(n, next + sq);
        }
      }
    }
  }
  // {最頻値のうちの一つ, その頻度}
  // 最頻値のうちの最大/最小値を求めることもできる
  std::pair<T, int> query(int l, int r){
    int lblock = (l + sq - 1) / sq, rblock = r / sq;
    if(lblock >= rblock){
      int mode = 0, modecnt = 0;
      for(int i = l; i < r; i++){
        int x = val[i], k = rank[i];
        if(k + modecnt >= pos[x].size() || pos[x][k + modecnt] >= r) continue;
        int j = k + modecnt;
        while(j < pos[x].size() && pos[x][j] < r) j++;
        if(modecnt < j - k) mode = x, modecnt = j - k;
        // else if(modecnt == j - k) mode = (rev[mode] < rev[x] ? x : mode)
      }
      return {rev[mode], modecnt};
    }
    auto [mode, modecnt] = large[lblock][rblock];
    for(int i = l; i < lblock * sq; i++){
      int x = val[i], k = rank[i];
      if(k + modecnt >= pos[x].size() || pos[x][k + modecnt] >= r) continue;
      int j = k + modecnt;
      while(j < pos[x].size() && pos[x][j] < r) j++;
      if(modecnt < j - k) mode = x, modecnt = j - k;
      // else if(modecnt == j - k) mode = (rev[mode] < rev[x] ? x : mode)
    }
    for(int i = r - 1; i >= rblock * sq; i--){
      int x = val[i], k = rank[i];
      if(k - modecnt < 0 || pos[x][k - modecnt] < l) continue;
      int j = k - modecnt;
      while(j >= 0 && pos[x][j] >= l) j--;
      if(modecnt < k - j) mode = x, modecnt = k - j;
      // else if(modecnt == j - k) mode = (rev[mode] < rev[x] ? x : mode)
    }
    return {rev[mode], modecnt};
  }
};
#endif