#ifndef _HORN_SAT_H_
#define _HORN_SAT_H_
#include <vector>
#include <unordered_set>
#include <cassert>
#include <queue>
// 各リテラルについてnotが付くものが最大1つまで
struct horn_sat{
  int n;
  std::vector<bool> answer;
private:
  struct clause{
    int neg;
    std::unordered_set<int> pos;
    clause(const std::vector<std::pair<int, bool>> &u): neg(-1){
      for(auto [i, f] : u){
        if(!f){
          assert(neg == -1 || neg == i); // notが2つ以上ある
          neg = i;
        }else{
          pos.insert(i);
        } 
      }
    }
    int size(){
      return (neg != -1) + pos.size();
    }
  };
  std::vector<clause> v;
  std::vector<std::pair<int, bool>> tmp;
public:
  horn_sat(int n): n(n), answer(n, 1){}
  // addをしていき, 次のリテラルに移りたい場合はend_clauseを呼ぶ
  void add(int i, bool f){
    assert(0 <= i && i < n);
    tmp.push_back({i, f});
  }
  void end_clause(){
    if(tmp.empty()) return;
    v.push_back(clause(tmp));
    tmp.clear();
  }
  bool satisfiable(){
    if(!tmp.empty()) end_clause();
    std::vector<std::vector<int>> pos_idx(n), neg_idx(n);
    int m = v.size();
    for(int i = 0; i < m; i++){
      if(v[i].neg != -1) neg_idx[v[i].neg].push_back(i);
      for(int x : v[i].pos) pos_idx[x].push_back(i);
    }
    std::queue<int> one; // 残り要素数1
    for(int i = 0; i < m; i++) if(v[i].size() == 1) one.push(i);
    while(!one.empty()){
      int i = one.front();
      one.pop();
      if(v[i].neg != -1){ // not1つのみ
        int x = v[i].neg;
        for(int j : pos_idx[x]){
          v[j].pos.erase(x);
          if(v[j].size() == 1) one.push(j);
          if(v[j].size() == 0) return false; // 矛盾
        }
        pos_idx[x].clear();
        answer[x] = 0;
      }else{ // pos1つのみ
        int x = *v[i].pos.begin();
        for(int j : neg_idx[x]){
          v[j].neg = -1;
          if(v[j].size() == 1) one.push(j);
          if(v[j].size() == 0) return false; // 矛盾
        }
        neg_idx[x].clear();
      }
    }
    return true;
  }
};
#endif