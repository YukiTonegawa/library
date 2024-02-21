#ifndef _MO_ALGORITHM_H_
#define _MO_ALGORITHM_H_
#include <vector>
#include <math.h>

template<typename mo_struct>
struct mo_algorithm{
private:
  struct query{
    int i, b, l, r, label;
    query(int i, int b, int l, int r, int label) : i(i), b(b), l(l), r(r), label(label){}
  };
  std::vector<query> q;
  int n, block_size, pq, pl, pr;
  mo_struct &_st;
  void build(){
    // block_size: https://nyaannyaan.github.io/library/misc/mo.hpp.html
    int qsz = std::max(1, (int)q.size());
    block_size = std::max(1, (int)(n * sqrt(1.5 / (double)qsz)));
    for(int i = 0; i < q.size(); i++) q[i].b /= block_size;
    std::sort(q.begin(), q.end(),  [](query &a, query &b){
      if(a.b != b.b) return a.b < b.b;
      return a.b & 1 ? a.r > b.r : a.r < b.r;
    });
  }
public:
  mo_algorithm(int n, mo_struct &_st) : n(n), pq(0), pl(0), pr(0), _st(_st){}
  inline int insert(int l, int r, int label = -1){
    int qsz = q.size();
    q.emplace_back(query(qsz, l,  l, r, label));
    return qsz;
  }
  // すでに全てのクエリを処理し終わっている場合は{-1, -1}
  // 今見ているクエリが追加された順番, ラベルを返す
  inline std::pair<int, int> process(){
    if(pq == 0) build();
    if(pq == q.size()) return {-1, -1};
    query &qi = q[pq];
    while(pl > qi.l) _st.add_left(--pl);
    while(pr < qi.r) _st.add_right(pr++);
    while(pl < qi.l) _st.del_left(pl++);
    while(pr > qi.r) _st.del_right(--pr);
    pq++;
    return {qi.i, qi.label};
  }
  // [l, r)が欲しい場合
  inline std::pair<int, int> process2(){
    if(pq == 0) build();
    if(pq == q.size()) return {-1, -1};
    query &qi = q[pq];
    while(pl > qi.l) _st.add_left(--pl, pr);
    while(pr < qi.r) _st.add_right(pl, pr++);
    while(pl < qi.l) _st.del_left(pl++, pr);
    while(pr > qi.r) _st.del_right(pl, --pr);
    pq++;
    return {qi.i, qi.label};
  }
};

/*
struct range_inversion{
  binary_indexed_tree<int> b;
  std::vector<int> a;
  int n;
  long long sum = 0;
  range_inversion(int n): n(n){}

  void add_left(int i){
    sum += b.query(0, a[i]);
    b.update(a[i], 1);
  }
  void del_left(int i){
    sum -= b.query(0, a[i]);
    b.update(a[i], -1);
  }
  void add_right(int i){
    sum += b.query(a[i] + 1, n);
    b.update(a[i], 1);
  }
  void del_right(int i){
    sum -= b.query(a[i] + 1, n);
    b.update(a[i], -1);
  }
};
*/
#endif