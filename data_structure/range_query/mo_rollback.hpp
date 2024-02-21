#ifndef _MO_ROLLBACK_H_
#define _MO_ROLLBACK_H_
#include <vector>
#include <cassert>
#include <array>
#include <math.h>
#include <algorithm>

template<typename mo_struct>
struct mo_rollback{
private:
  struct query{
    int i, b, l, r, label;
    query(int i, int b, int l, int r, int label) : i(i), b(b), l(l), r(r), label(label){}
  };
  std::vector<query> small_q, large_q;
  const int block_size;
  int pq, pl, pr, n, border;
  bool is_shot;
  mo_struct &_st;
  void build(){
    std::sort(large_q.begin(), large_q.end(),  [](query &a, query &b){
      if(a.b != b.b) return a.b < b.b;
      return a.r < b.r;
    });
  }
public:
  mo_rollback(int n, mo_struct &_st) : block_size(round(sqrt(n))), pq(0), pl(0), pr(0), n(n), is_shot(false), _st(_st){}
  inline int insert(int l, int r, int label = -1){
    assert(0 <= l && r <= n && l <= r);
    int qid = large_q.size() + small_q.size();
    if(r - l > block_size) large_q.emplace_back(query(qid, l / block_size,  l, r, label));
    else small_q.emplace_back(query(qid, l / block_size,  l, r, label));
    return qid;
  }
  // すでに全てのクエリを処理し終わっている場合{-1, -1, 0}
  // {クエリの追加された順番, ラベル, 区間長が√N以上のクエリか}
  inline std::tuple<int, int, bool> process(){
    if(pq == 0) build();
    if(pq < large_q.size()){
      query qi = large_q[pq];
      if(!pq || large_q[pq - 1].b < qi.b){
        _st.reset();
        pr = pl = border = std::min(n, (qi.b + 1) * block_size);
        is_shot = false;
      }
      pq++;
      if(is_shot) _st.rollback(), pl = border, is_shot = false;
      while(pr < qi.r) _st.add_right(pr++);
      _st.snapshot();
      is_shot = true;
      while(pl > qi.l) _st.add_left(--pl);
      return {qi.i, qi.label, 1};
    }else if(pq < large_q.size() + small_q.size()){
      query qi = small_q[pq - (int)large_q.size()];
      pq++;
      _st.naive_solve(qi.l, qi.r);
      return {qi.i, qi.label, 0};
    }else return {-1, -1, 0};
  }
  // l, rが欲しい場合
  inline std::tuple<int, int, bool> process2(){
    if(pq == 0) build();
    if(pq < large_q.size()){
      query qi = large_q[pq];
      if(!pq || large_q[pq - 1].b < qi.b){
        _st.reset();
        pr = pl = border = std::min(n, (qi.b + 1) * block_size);
        is_shot = false;
      }
      pq++;
      if(is_shot) _st.rollback(), pl = border, is_shot = false;
      while(pr < qi.r) _st.add_right(pl, pr++);
      _st.snapshot();
      is_shot = true;
      while(pl > qi.l) _st.add_left(--pl, pr);
      return {qi.i, qi.label, 1};
    }else if(pq < large_q.size() + small_q.size()){
      query qi = small_q[pq - (int)large_q.size()];
      pq++;
      _st.naive_solve(qi.l, qi.r);
      return {qi.i, qi.label, 0};
    }else return {-1, -1, 0};
  }
};

/*
struct example{
  void snapshot(){

  }
  void rollback(){

  }
  void add_left(const int i){

  }
  void add_right(const int i){

  }
  // 呼ばれる回数は高々√N回
  void reset(){

  }
  void naive_solve(int l, int r){

  }
};
*/

/*

struct example{
  std::vector<int> A;
  unsafe_deque_rollback<int, 200000> ppq;
  void snapshot(){
    ppq.bookmark();
  }
  void rollback(){
    ppq.rollback_bookmark();
  }
  void add_left(const int i){
    if(ppq.size() <= 1){
      ppq.push_front(A[i]);
      return;
    }
    int y = ppq.front();
    if(A[i] >= y && y >= ppq.get(1)){
      ppq.set(0, A[i]);
    }else{
      ppq.push_front(A[i]);
    }
  }
  void add_right(const int i){
    int sz = ppq.size();
    if(sz <= 1){
      ppq.push_back(A[i]);
      return;
    }
    int y = ppq.back();
    if(ppq.get(sz - 2) >= y && y >= A[i]){
      ppq.set(sz - 1, A[i]);
    }else{
      ppq.push_back(A[i]);
    }
  }
  // 呼ばれる回数は高々√N回
  void reset(){
    ppq.reset();
  }
  int len;
  void naive_solve(int l, int r){
    int r1, r2;
    len = 0;
    for(int i = l; i < r; i++){
      if(len == 0){
        r1 = r2 = A[i], len = 1;
      }else if(len == 1){
        r2 = r1, r1 = A[i], len = 2;
      }else if(r2 >= r1){
        if(r1 >= A[i]){
          r1 = A[i];
        }else{
          r2 = r1, r1 = A[i], len++;
        }
      }else{
        r2 = r1, r1 = A[i], len++;
      }
    }
  }
};

*/
#endif