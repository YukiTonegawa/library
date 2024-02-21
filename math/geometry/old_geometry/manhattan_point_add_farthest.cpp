#include <vector>
#include <set>
#include <array>
#include <tuple>
#include <queue>

template<typename Idx, int dim = 2>
struct manhattan_point_add_farthest{
private:
  const static int manhattan_dim = dim;
  const static int chebyshev_dim = 1 << (dim - 1);
  using manhattan_point = std::array<Idx, dim>;
  using chebyshev_point = std::array<Idx, chebyshev_dim>;
  std::vector<std::pair<Idx, Idx>> max_min;
  inline chebyshev_point manhattan2chebyshev(manhattan_point &p){
    if(dim == 2) return {p[0] + p[1], p[0] - p[1]};
    if(dim == 3) return {p[0] + p[1] + p[2], p[0] + p[1] - p[2], p[0] - p[1] + p[2], p[0] - p[1] - p[2]};
    chebyshev_point ret;
    std::queue<std::tuple<Idx, int, int>> q;
    q.push({p[0], 1, 0});
    while(!q.empty()){
      auto [val, k, idx] = q.front();
      q.pop();
      if(k == manhattan_dim){
        ret[idx] = val;
        continue;
      }
      q.push({val + p[k], k + 1, idx});
      q.push({val - p[k], k + 1, idx + (1 << (manhattan_dim - k - 1))});
    }
    return ret;
  }
public:
  manhattan_point_add_farthest(){}
  void add_point(manhattan_point p){
    if(max_min.empty()){
      max_min.resize(chebyshev_dim, {std::numeric_limits<Idx>::min(), std::numeric_limits<Idx>::max()});
    }
    chebyshev_point q = manhattan2chebyshev(p);
    for(int i=0;i<chebyshev_dim;i++){
      max_min[i].first = std::max(max_min[i].first, q[i]);
      max_min[i].second = std::min(max_min[i].second, q[i]);
    }
  }
  Idx max_dist(manhattan_point p){
    assert(!max_min.empty());
    chebyshev_point q = manhattan2chebyshev(p);
    Idx ret = std::numeric_limits<Idx>::min();
    for(int i=0;i<chebyshev_dim;i++){
      ret = std::max(ret, abs(max_min[i].first - q[i]));
      ret = std::max(ret, abs(max_min[i].second - q[i]));
    }
    return ret;
  }
};

template<typename Idx, int dim = 2>
struct manhattan_point_add_farthest_erasable{
private:
  const static int manhattan_dim = dim;
  const static int chebyshev_dim = 1 << (dim - 1);
  using manhattan_point = std::array<Idx, dim>;
  using chebyshev_point = std::array<Idx, chebyshev_dim>;
  std::vector<std::multiset<Idx>> max_min;
  inline chebyshev_point manhattan2chebyshev(manhattan_point &p){
    if(dim == 2) return {p[0] + p[1], p[0] - p[1]};
    if(dim == 3) return {p[0] + p[1] + p[2], p[0] + p[1] - p[2], p[0] - p[1] + p[2], p[0] - p[1] - p[2]};
    chebyshev_point ret;
    std::queue<std::tuple<Idx, int, int>> q;
    q.push({p[0], 1, 0});
    while(!q.empty()){
      auto [val, k, idx] = q.front();
      q.pop();
      if(k == manhattan_dim){
        ret[idx] = val;
        continue;
      }
      q.push({val + p[k], k + 1, idx});
      q.push({val - p[k], k + 1, idx + (1 << (manhattan_dim - k - 1))});
    }
    return ret;
  }
public:
  manhattan_point_add_farthest_erasable(){}
  void add_point(manhattan_point &p){
    if(max_min.empty()){
      max_min.resize(chebyshev_dim);
    }
    chebyshev_point q = manhattan2chebyshev(p);
    for(int i=0;i<chebyshev_dim;i++){
      max_min[i].insert(q[i]);
    }
  }
  Idx erase(manhattan_point p){
    assert(!max_min.empty());
    chebyshev_point q = manhattan2chebyshev(p);
    for(int i=0;i<chebyshev_dim;i++){
      max_min[i].erase(max_min[i].find(q[i]));
    }
  }
  Idx max_dist(manhattan_point p){
    assert(!max_min.empty());
    chebyshev_point q = manhattan2chebyshev(p);
    Idx ret = std::numeric_limits<Idx>::min();
    for(int i=0;i<chebyshev_dim;i++){
      ret = std::max(ret, abs(*max_min[i].begin() - q[i]));
      ret = std::max(ret, abs(*(--max_min[i].end()) - q[i]));
    }
    return ret;
  }
};
