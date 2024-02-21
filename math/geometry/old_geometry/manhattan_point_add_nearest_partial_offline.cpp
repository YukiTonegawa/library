#include <vector>
#include <tuple>
#include <algorithm>
.
template<typename Val>
struct bit_min{
  int M;
  Val INF = std::numeric_limits<Val>::max();
  std::vector<Val> v;
  bit_min(){}
  bit_min(int _n): M(_n), v(M+1, INF){}
  bit_min(const std::vector<Val> &val): M(val.size()), v(1){
    v.insert(v.begin()+1, val.begin(), val.end());
    for(int i=1;i<=val.size();i++){
      int nxt = i + (i&(-i));
      if(nxt<=M) v[nxt] = std::min(v[nxt], v[i]);
    }
  }
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
template<typename Idx = int, typename Val = Idx>
struct manhattan_point_add_nearest{
private:
  int M = 1;
  std::vector<Idx> X;
  std::vector<std::vector<Idx>> Y, Y_rev;
  std::vector<bit_min<Val>> pp, pm, mp, mm;
  using point = std::tuple<Idx, Idx, Val>;
public:
  manhattan_point_add_nearest(){}
  manhattan_point_add_nearest(std::vector<point> v){
    int n = v.size();
    std::sort(v.begin(), v.end(), [](const point &a, const point &b){return std::get<1>(a) < std::get<1>(b);});
    for(int i=0;i<n;i++) X.push_back(std::get<0>(v[i]));
    std::sort(X.begin(), X.end());
    X.erase(std::unique(X.begin(), X.end()), X.end());
    M = (int)X.size();
    std::vector<std::vector<Val>> tmp_mp(M+1), tmp_mm(M+1), tmp_pp(M+1), tmp_pm(M+1);
    pp.resize(M+1), pm.resize(M+1), mp.resize(M+1), mm.resize(M+1);
    Y.resize(M+1), Y_rev.resize(M+1);

    for(int i=0;i<n;i++){
      auto [x, y, z] = v[i];
      int k = lower_bound(X.begin(), X.end(), x) - X.begin();
      for(int j=k+1;j<=M;j+=(j&(-j))){
        if(Y[j].empty()||Y[j].back()!=y){
          Y[j].push_back(y);
          tmp_mp[j].push_back(-x+y+z);
          tmp_mm[j].push_back(-x-y+z);
        }else{
          tmp_mp[j].back() = std::min(tmp_mp[j].back(), -x+y+z);
          tmp_mm[j].back() = std::min(tmp_mm[j].back(), -x-y+z);
        }
      }
      for(int j=X.size()-k;j<=M;j+=(j&(-j))){
        if(Y_rev[j].empty()||Y_rev[j].back()!=y){
          Y_rev[j].push_back(y);
          tmp_pp[j].push_back(x+y+z);
          tmp_pm[j].push_back(x-y+z);
        }else{
          tmp_pp[j].back() = std::min(tmp_pp[j].back(), x+y+z);
          tmp_pm[j].back() = std::min(tmp_pm[j].back(), x-y+z);
        }
      }
    }
    for(int i=0;i<=M;i++){
      std::reverse(tmp_mp[i].begin(), tmp_mp[i].end());
      std::reverse(tmp_pp[i].begin(), tmp_pp[i].end());
      mp[i] = bit_min(tmp_mp[i]);
      mm[i] = bit_min(tmp_mm[i]);
      pp[i] = bit_min(tmp_pp[i]);
      pm[i] = bit_min(tmp_pm[i]);
    }
  }
  void add_point(Idx x, Idx y, Val z = 0){
    int xidx = lower_bound(X.begin(), X.end(), x) - X.begin();
    for(int i=xidx+1;i<=M;i+=(i&(-i))){
      int yidx = lower_bound(Y[i].begin(), Y[i].end(), y) - Y[i].begin();
      mp[i].update(Y[i].size() - yidx - 1, -x+y+z);
      mm[i].update(yidx, -x-y+z);
    }
    for(int i=X.size()-xidx;i<=M;i+=(i&(-i))){
      int yidx = lower_bound(Y_rev[i].begin(), Y_rev[i].end(), y) - Y_rev[i].begin();
      pp[i].update(Y_rev[i].size() - yidx - 1, x+y+z);
      pm[i].update(yidx, x-y+z);
    }
  }
  Val min_dist(Idx x, Idx y){
    Val ret = std::numeric_limits<Val>::max();
    int xidx = lower_bound(X.begin(), X.end(), x) - X.begin();
    for(int i=xidx;i>0;i-=(i&(-i))){
      int yidx = lower_bound(Y[i].begin(), Y[i].end(), y) - Y[i].begin();
      Val v = mm[i].min(yidx);
      if(v != mm[i].INF) ret = std::min(ret, v + x + y);
      v = mp[i].min(Y[i].size() - yidx);
      if(v != mp[i].INF) ret = std::min(ret, v + x - y);
    }
    for(int i=(X.size()-xidx);i>0;i-=(i&(-i))){
      int yidx = lower_bound(Y_rev[i].begin(), Y_rev[i].end(), y) - Y_rev[i].begin();
      Val v = pm[i].min(yidx);
      if(v != pm[i].INF) ret = std::min(ret, v - x + y);
      v = pp[i].min(Y_rev[i].size() - yidx);
      if(v != pp[i].INF) ret = std::min(ret, v - x - y);
    }
    return ret;
  }
};
