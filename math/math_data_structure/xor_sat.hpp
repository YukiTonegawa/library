#ifndef _XOR_SAT_H_
#define _XOR_SAT_H_
#include <vector>
#include <unordered_set>
#include <cassert>
#include "../matrix/matrix_mod2.hpp"

// 各リテラルの条件を満たす個数の偶奇を指定する
struct xor_sat{
  int n;
  std::vector<bool> answer;
private:
  struct clause{
    std::vector<std::pair<int, bool>> u;
    bool f;
    clause(const std::vector<std::pair<int, bool>> &u, bool f): u(u), f(f){}
  };
  std::vector<clause> v;
  std::vector<std::pair<int, bool>> tmp;
public:
  xor_sat(int n): n(n), answer(n){assert(n);}
  // addをしていき, 次のリテラルに移りたい場合はend_clauseを呼ぶ
  void add(int i, bool f){
    assert(0 <= i && i < n);
    tmp.push_back({i, f});
  }
  // 与えたリテラルのうちfが1なら奇数個, 0なら偶数個が条件を満たす
  void end_clause(bool f){
    if(tmp.empty()) return;
    v.push_back(clause(tmp, f));
    tmp.clear();
  }
  bool satisfiable(){
    int h = v.size();
    if(h == 0){
      std::fill(answer.begin(), answer.end(), true);
      return true;
    }
    matrix_mod2 mat(h, n);
    dynamic_bitset bit(h, 0);
    for(int i = 0; i < h; i++){
      bool e = v[i].f;
      for(auto [j, f] : v[i].u){
        bool g = mat.get(i, j);
        if(!f) e = !e;
        mat.set(i, j, !g);
      }
      bit.set(i, e);
    }
    auto [d, ans] = mat.system_of_linear_equations(bit);
    if(ans.size() == 0) return false; // 解なし
    for(int i = 0; i < n; i++) answer[i] = ans.get(i);
    return true;
  }
};
#endif