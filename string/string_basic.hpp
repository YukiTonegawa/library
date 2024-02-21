#ifndef _STRING_BASIC_H_
#define _STRING_BASIC_H_
#include <vector>
#include <string>

std::vector<char> stovector(const std::string &s){
  int n = s.size();
  std::vector<char> res(n);
  for(int i = 0; i < n; i++) res[i] = s[i];
  return res;
}
// cを0としてintに直す
std::vector<int> stovector_int(const std::string &s, char c = 'a'){
  int n = s.size();
  std::vector<int> res(n);
  for(int i = 0; i < n; i++) res[i] = s[i] - c;
  return res;
}
#include <sstream>
std::vector<std::string> split(const std::string &s, char delim){
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)){
    if(!item.empty()) elems.push_back(item);
  }
  return elems;
}

template<typename T>
std::vector<std::pair<T, int>> runlength(const std::vector<T> &v){
  std::vector<std::pair<T, int>> ret;
  int n = v.size();
  for(int i = 0; i < n; i++){
    if(ret.empty() || ret.back().first != v[i]){
      ret.push_back({v[i], 1});
    }else{
      ret.back().second++;
    }
  }
  return ret;
}
// 全ての要素が [0, NUM_VAL)
template<int NUM_VAL = 26>
struct string_processor{
private:
  static constexpr int s = 64, sdiv = 6, smod = 63;
  using Z = uint64_t;
  int N;
  std::vector<std::array<int, NUM_VAL>> B;
  std::vector<std::array<Z, NUM_VAL>> S;
  std::vector<std::vector<int>> V;
public:
  string_processor(){}
  string_processor(const std::vector<int> &v): N(v.size()){
    int M = (N + s - 1) / s;
    V.resize(NUM_VAL, std::vector<int>());
    std::array<int, NUM_VAL> pop;
    std::array<Z, NUM_VAL> pop_small;
    pop.fill(0);
    pop_small.fill(0);
    B.resize(M + 1, pop);
    S.resize(M, pop_small);
    for(int i = 0, t = 0, sz = 0; i < N; i++, t++){
      int x = v[i];
      assert(0 <= x && x < NUM_VAL);
      V[x].push_back(i);
      pop[x]++, pop_small[x] |= (Z(1) << t);
      if(t == s - 1 || i == N - 1){
        for(int j = 0; j < NUM_VAL; j++){
          if(j) pop[j] += pop[j - 1], pop_small[j] |= pop_small[j - 1];
          B[sz + 1][j] = pop[j] + B[sz][j];
          S[sz][j] = pop_small[j];
        }
        pop.fill(0);
        pop_small.fill(0);
        t = -1;
        sz++;
      }
    }
  }
  // count c, i < r, O(1)
  int rank(int r, int c){
    if(c == 0) return rank_lower(r, c);
    assert(0 <= r && r <= N);
    assert(0 <= c && c < NUM_VAL);
    int rq = r >> sdiv, rm = r & smod;
    int ret = B[rq][c] - B[rq][c - 1];
    if(rm) ret += __builtin_popcountll((S[rq][c] ^ S[rq][c - 1]) << (s - rm));
    return ret;
  }
  // count [0, c], i < r, O(1)
  int rank_lower(int r, int c){
    assert(0 <= r && r <= N);
    assert(0 <= c && c < NUM_VAL);
    int rq = r >> sdiv, rm = r & smod;
    int ret = B[rq][c];
    if(rm) ret += __builtin_popcountll(S[rq][c] << (s - rm));
    return ret;
  }
  // k番目のc, ない場合は-1 O(1)
  int select(int k, int c){
    assert(0 <= c && c < NUM_VAL);
    return (k < 0 || V[c].size() <= k) ? -1 : V[c][k];
  }
  // k番目のc以下, ない場合は-1 O(logN)
  int select_lower(int k, int c){
    assert(0 <= c && c < NUM_VAL);
    int l = 0, r = N + 1;
    while(r - l > 1){
      int mid = (l + r) >> 1;
      if(rank_lower(mid, c) <= k) l = mid;
      else r = mid;
    }
    return r == N + 1 ? -1 : l;
  }
  // [i, n)で最も左のc
  int find_next(int i, int c){
    return select(rank(i, c), c);
  }
  // [0, i]で最も右のc
  int find_prev(int i, int c){
    return select(rank(i + 1, c) - 1, c);
  }
};
#endif