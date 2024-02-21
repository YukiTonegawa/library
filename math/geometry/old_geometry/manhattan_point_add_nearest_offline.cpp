#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
.
template<typename Val>
struct bit_min{
  int M;
  const Val INF = std::numeric_limits<Val>::max();
  std::vector<Val> v;
  bit_min(int _n): M(_n), v(M+1, INF){}
  void update(int k, Val x){
    for(int i=k+1;i<=M;i+=(i&(-i))){
      if(v[i] <= x) return;
      v[i] = x;
    }
  }
  Val min(int r){
    Val ret = INF;
    for(int k=r;k>0;k-=(k&(-k))) ret = std::min(ret, v[k]);
    return ret;
  }
};
template<typename Idx, typename Val = Idx>
struct manhattan_point_add_nearest{
private:
  using query = std::tuple<int, Idx, Idx, Val>;//index, x, y, value
  using query2 = std::tuple<Idx, int, int, Idx, int, Val>;//x, query_type, index, y, y_idx, value
  std::vector<query> q;
  std::vector<int> question{0};
  int qs = 0;
  void solve(int l, int r, std::vector<Val> &ans){
    if(r - l < 2) return;
    int mid = (l + r) >> 1;
    solve(l, mid, ans);
    solve(mid, r, ans);
    int left_insert = (mid - l) - (question[mid] - question[l]);
    int right_question = question[r] - question[mid];
    if(left_insert == 0 || right_question == 0) return;
    std::vector<Idx> y;
    for(int i=l;i<mid;i++) if(std::get<0>(q[i]) == -1) y.push_back(std::get<2>(q[i]));//point add
    std::sort(y.begin(), y.end());
    y.erase(std::unique(y.begin(), y.end()), y.end());
    bit_min<Val> pp(y.size()), pm(y.size()), mp(y.size()), mm(y.size());
    std::vector<query2> tmp_q;
    for(int i=l;i<mid;i++){
      if(std::get<0>(q[i]) == -1){
        auto [idx, x_tmp, y_tmp, val] = q[i];
        int y_idx = std::lower_bound(y.begin(), y.end(), y_tmp) - y.begin();
        tmp_q.push_back(std::make_tuple(x_tmp, 2, -1, y_tmp, y_idx, val));
        tmp_q.push_back(std::make_tuple(-x_tmp, 3, -1, y_tmp, y_idx, val));
      }
    }
    for(int i=mid;i<r;i++){
      if(std::get<0>(q[i]) != -1){
        auto [idx, x_tmp, y_tmp, _] = q[i];
        int y_idx = std::lower_bound(y.begin(), y.end(), y_tmp) - y.begin();
        tmp_q.push_back(std::make_tuple(x_tmp, 0, idx, y_tmp, y_idx, 0));
        tmp_q.push_back(std::make_tuple(-x_tmp, 1, idx, y_tmp, y_idx, 0));
      }
    }
    std::sort(tmp_q.begin(), tmp_q.end());
    std::reverse(tmp_q.begin(), tmp_q.end());
    for(auto [x_tmp, query_type, idx, y_tmp, y_idx, val]:tmp_q){
      if(query_type == 0){
        Val v = pp.min(y.size() - y_idx);
        if(v != pp.INF) ans[idx] = std::min(ans[idx], v - x_tmp - y_tmp);
        v = pm.min(y_idx);
        if(v != pm.INF) ans[idx] = std::min(ans[idx], v - x_tmp + y_tmp);
      }else if(query_type == 1){
        Val v = mp.min(y.size() - y_idx);
        if(v != mp.INF) ans[idx] = std::min(ans[idx], v - x_tmp - y_tmp);
        v = mm.min(y_idx);
        if(v != mm.INF) ans[idx] = std::min(ans[idx], v - x_tmp + y_tmp);
      }else if(query_type == 2){
        pp.update(y.size() - y_idx - 1, x_tmp + y_tmp + val);
        pm.update(y_idx, x_tmp - y_tmp + val);
      }else if(query_type == 3){
        mp.update(y.size() - y_idx - 1, x_tmp + y_tmp + val);
        mm.update(y_idx, x_tmp - y_tmp + val);
      }
    }
  }
public:
  manhattan_point_add_nearest(){}
  void add_point(Idx x, Idx y, Val z = 0){
    q.push_back(std::make_tuple(-1, x, y, z));
    question.push_back(0);
  }
  void min_dist(Idx x, Idx y){
    q.push_back(std::make_tuple(qs++, x, y, 0));
    question.push_back(1);
  }
  std::vector<Val> solve(){
    std::vector<Val> ret(qs, std::numeric_limits<Val>::max());
    for(int i=1;i<question.size();i++) question[i] += question[i-1];
    solve(0, q.size(), ret);
    return ret;
  }
};
