#ifndef _SQRT_DECOMPOSITION_H_
#define _SQRT_DECOMPOSITION_H_
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>

template<typename Val>
struct sqrtd_range_sum{
  int B, M, N;
  std::vector<Val> large, small;
  sqrtd_range_sum(int n): B(std::max(1, (int)sqrt(n / 2))), M((n + B - 1) / B), N(M * B), large(M, 0), small(N, 0){}
  // Ai <- Ai + x  O(1)
  void add(int i, Val x){
    large[i / B] += x;
    small[i] += x;
  }
  // Ai <- x  O(1)
  void set(int i, Val x){
    large[i / B] += x - small[i];
    small[i] = x;
  }
  // O(1)
  Val get(int i){
    return small[i];
  }
  // sum [l, r) O(√N)
  Val query(int l, int r){
    int lb = l / B, rb = r / B;
    Val res = 0;
    if(lb == rb){
      for(int i = l; i < r; i++) res += small[i];
      return res;
    }
    for(int i = lb; i < rb; i++) res += large[i];
    for(int i = lb * B; i < l; i++) res -= small[i];
    for(int i = rb * B; i < r; i++) res += small[i];
    return res;
  }
  // sum [0, r) O(√N)
  Val query(int r){
    int rb = r / B;
    Val res = 0;
    for(int i = 0; i < rb; i++) res += large[i];
    for(int i = rb * B; i < r; i++) res += small[i];
    return res;
  }
  // 全ての要素が0以上の時のみ使える
  // sum[l, r]がx以上になる最左の{r, sum}, ない場合は{-1, [l, n)のsum} O(√N)
  std::pair<int, Val> lower_bound(int l, Val x){
    Val s = 0;
    int m = (l + B - 1) / B * B;
    while(l < m){
      s += small[l];
      if(s >= x) return {l, s};
      l++;
    }
    for(l /= B; l < M; l++){
      if(s + large[l] >= x){
        for(int i = l * B;; i++){
          s += small[i];
          if(s >= x) return {i, s};
        }
      }
      s += large[l];
    }
    return {-1, s};
  }
  // 全ての要素が0か1の時のみ使える
  // l以降(l含む)で初めて現れる0の位置 O(√N)
  int next_zero(int l){
    Val s = 0;
    int m = (l + B - 1) / B * B;
    while(l < m){
      if(small[l] == 0) return {l, s};
      s += small[l];
      l++;
    }
    for(l /= B; l < M; l++){
      if(large[l] < B){        
        for(int i = l * B;; i++){
          if(small[i] == 0) return {i == N ? -1 : i, s};
          s += small[i];
        }
      }
      s += large[l];
    }
    return {-1, s};
  }
};

template<typename Val>
struct sqrtd_range_add_range_sum{
  int B, M, N, _n;
  std::vector<Val> large, large_slope, small, small_slope;
  sqrtd_range_add_range_sum(int n): B(std::max(1, (int)sqrt(n / 2))), M((n + B - 1) / B), N(M * B), _n(n), large(M, 0), large_slope(M, 0), small(N, 0), small_slope(N, 0){}
  // Ai <- Ai + x  O(1)
  void add(int i, Val x){
    large[i / B] += x;
    small[i] += x;
  }
  // Ai <- Ai + x (i in [l, r))  O(1)
  void add(int l, int r, Val x){
    int lb = l / B;
    large_slope[lb] += x;
    small_slope[l] += x;
    large[lb] -= x * l;
    small[l] -= x * l;
    if(r < N){
      int rb = r / B;
      large_slope[rb] -= x;
      small_slope[r] -= x;
      large[rb] += x * r;
      small[r] += x * r;
    }
  }
  // sum [l, r) O(√N)
  Val query(int l, int r){
    assert(l <= r);
    int lb = l / B, rb = r / B;
    Val res = 0, ssum = 0;
    for(int i = 0; i < lb; i++) ssum += large_slope[i];
    Val tmpssum = 0;
    for(int i = lb * B; i < l; i++){
      tmpssum += small_slope[i];
      res -= small[i];
    }
    res -= (ssum + tmpssum) * l;
    for(int i = lb; i < rb; i++){
      res += large[i];
      ssum += large_slope[i];
    }
    for(int i = rb * B; i < r; i++){
      res += small[i];
      ssum += small_slope[i];
    }
    return res + ssum * r;
  }
  // sum [0, r) O(√N)
  Val query(int r){
    int rb = r / B;
    Val res = 0, ssum = 0;
    for(int i = 0; i < rb; i++){
      res += large[i];
      ssum += large_slope[i];
    }
    for(int i = rb * B; i < r; i++){
      res += small[i];
      ssum += small_slope[i];
    }
    return res + ssum * r;
  }
  // O(n)
  std::vector<Val> to_list(){
    std::vector<Val> res(_n, 0);
    Val sum = 0, ssum = 0;
    for(int i = 0; i < _n; i++){
      sum += small[i];
      ssum += small_slope[i];
      res[i] = sum + ssum * (i + 1);
    }
    for(int i = _n - 1; i > 0; i--) res[i] -= res[i - 1];
    return res;
  }
  // 全ての要素が0以上の時のみ使える
  // sum[l, r]がx以上になる最左の{r, sum}, ない場合は{-1, [l, n)のsum} O(√N)
  int lower_bound(int l, Val x){
    Val ssum = 0;
    int lb = l / B;
    for(int i = 0; i < lb; i++) ssum += large_slope[i];
    Val tmpsum = 0;
    for(int i = lb * B; i < l; i++){
      tmpsum += small_slope[i];
      x += small[i];
    }
    x += (ssum + tmpsum) * l;
    int r = (lb + 1) * B;
    for(int i = lb; i < M; i++){
      if((x - large[i]) <= (ssum + large_slope[i]) * r){
        for(int j = lb * B; ; j++){
          x -= small[j];
          ssum += small_slope[j];
          if(x <= ssum * (j + 1)) return (j >= _n ? -1: j);
        }
      }
      x -= large[i];
      ssum += large_slope[i];
    }
    return -1;
  }
};
template<typename Val>
struct sqrtd_plus_minus_one_range_min{
  static constexpr Val inf = std::numeric_limits<Val>::max() / 2;
  int B, M, N;
  struct node{
    Val x;
    int cnt;
    node *prev, *next;
    node(Val _x, int _cnt): x(_x), cnt(_cnt), prev(nullptr), next(nullptr){}
  };
  std::vector<Val> large;
  std::vector<node*> small;
  sqrtd_plus_minus_one_range_min(const std::vector<Val> &v): B(std::max(1, (int)sqrt(v.size() / 2))), M((v.size() + B - 1) / B), N(M * B), large(M, inf), small(N){
    for(int i = 0; i < M; i++){
      int l = i * B, r = (i + 1) * B;
      std::vector<std::pair<Val, int>> tmp;
      for(int j = l; j < r; j++){
        Val x = (j < v.size() ? v[j] : inf);
        large[i] = std::min(large[i], x);
        tmp.push_back({x, j});
      }
      std::sort(tmp.begin(), tmp.end());
      node *pre = nullptr;
      for(int j = 0; j < r - l; j++){
        if(!j || tmp[j].first > tmp[j - 1].first){
          auto u = new node(tmp[j].first, 1);
          small[tmp[j].second] = u;
          u->prev = pre;
          if(pre) pre->next = u;
          pre = u;
        }else{
          pre->cnt++;
          small[tmp[j].second] = pre;
        }
      }
    }
  }
  // Ai <- Ai + 1  O(1)
  void add_one(int i){
    auto u = small[i];
    int b = i / B;
    if(u->x == large[b] && u->cnt == 1) large[b]++;
    if(u->cnt == 1){
      if(u->next && u->next->x == u->x + 1){
        u->next->prev = u->prev;
        if(u->prev) u->prev->next = u->next;
        small[i] = u->next;
        small[i]->cnt++;
        delete u;
      }else{
        u->x++;
      }
    }else{
      u->cnt--;
      if(u->next && u->next->x == u->x + 1){
        small[i] = u->next;
        small[i]->cnt++;
      }else{
        small[i] = new node(u->x + 1, 1);
        small[i]->next = u->next;
        if(u->next) u->next->prev = small[i];
        u->next = small[i];
        small[i]->prev = u;
      }
    }
  }
  // Ai <- Ai - 1  O(1)
  void sub_one(int i){
    auto u = small[i];
    int b = i / B;
    if(u->x == large[b]) large[b]--;
    if(u->cnt == 1){
      if(u->prev && u->prev->x == u->x - 1){
        u->prev->next = u->next;
        if(u->next) u->next->prev = u->prev;
        small[i] = u->prev;
        small[i]->cnt++;
        delete u;
      }else{
        u->x--;
      }
    }else{
      u->cnt--;
      if(u->prev && u->prev->x == u->x - 1){
        small[i] = u->prev;
        small[i]->cnt++;
      }else{
        small[i] = new node(u->x - 1, 1);
        small[i]->prev = u->prev;
        if(u->prev) u->prev->next = small[i];
        u->prev = small[i];
        small[i]->next = u;
      }
    }
  }
  // min [l, r) O(√N)
  Val query(int l, int r){
    assert(l < r);
    int lb = (l + B - 1) / B, rb = r / B;
    Val res = inf;
    if(lb >= rb){
      for(int i = l; i < r; i++) res = std::min(res, small[i]->x);
      return res;
    }
    for(int i = lb; i < rb; i++) res = std::min(res, large[i]);
    for(int i = l; i < lb * B; i++) res = std::min(res, small[i]->x);
    for(int i = rb * B; i < r; i++) res = std::min(res, small[i]->x);
    return res;
  }
  // min [0, r) O(√N)
  Val query(int r){
    assert(0 < r);
    int rb = r / B;
    Val res = inf;
    for(int i = 0; i < rb; i++) res = std::min(res, large[i]);
    for(int i = rb * B; i < r; i++) res = std::min(res, small[i]->x);
    return res;
  }
  std::pair<int, Val> query_idx(int l, int r){
    assert(l < r);
    int lb = (l + B - 1) / B, rb = r / B;
    Val res = inf;
    int res_idx = 0;
    if(lb >= rb){
      for(int i = l; i < r; i++) if(res > small[i]->x) res = small[i]->x, res_idx = i;
      return {res_idx, res};
    }
    for(int i = lb; i < rb; i++) if(res > large[i]) res = large[i], res_idx = N + i;
    for(int i = l; i < lb * B; i++) if(res > small[i]->x) res = small[i]->x, res_idx = i;
    for(int i = rb * B; i < r; i++) if(res > small[i]->x) res = small[i]->x, res_idx = i;
    if(res_idx >= N){
      for(int i = (res_idx - N) * B;; i++){
        if(small[i]->x == res) return {i, res};
      }
    }
    return {res_idx, res};
  }
  // min [l, r]がx以下になる最左の{r, min}, ない場合は{-1, [l, n)のmin} O(√N)
  std::pair<int, Val> lower_bound(int l, Val x){
    Val s = inf;
    int m = (l + B - 1) / B * B;
    while(l < m){
      s = std::min(s, small[l]->x);
      if(s <= x) return {l, s};
      l++;
    }
    for(l /= B; l < M; l++){
      if(large[l] <= x){
        for(int i = l * B;; i++){
          s = std::min(s, small[i]);
          if(s <= x) return {i, s};
        }
      }
      s = std::min(s, large[l]);
    }
    return {-1, s};
  }
};
template<typename Val>
struct plus_minus_one_all_min{
  static constexpr Val inf = std::numeric_limits<Val>::max() / 2;
  int N;
  struct node{
    Val x;
    int cnt;
    node *prev, *next;
    node(Val _x, int _cnt): x(_x), cnt(_cnt), prev(nullptr), next(nullptr){}
  };
  Val large;
  std::vector<node*> small;
  plus_minus_one_all_min(const std::vector<Val> &v): N(v.size()), large(inf), small(N){
    std::vector<std::pair<Val, int>> tmp;
    for(int j = 0; j < N; j++){
      Val x = v[j];
      large = std::min(large, x);
      tmp.push_back({x, j});
    }
    std::sort(tmp.begin(), tmp.end());
    node *pre = nullptr;
    for(int j = 0; j < N; j++){
      if(!j || tmp[j].first > tmp[j - 1].first){
        auto u = new node(tmp[j].first, 1);
        small[tmp[j].second] = u;
        u->prev = pre;
        if(pre) pre->next = u;
        pre = u;
      }else{
        pre->cnt++;
        small[tmp[j].second] = pre;
      }
    }
  }
  // Ai <- Ai + 1  O(1)
  void add_one(int i){
    auto u = small[i];
    if(u->x == large && u->cnt == 1) large++;
    if(u->cnt == 1){
      if(u->next && u->next->x == u->x + 1){
        u->next->prev = u->prev;
        if(u->prev) u->prev->next = u->next;
        small[i] = u->next;
        small[i]->cnt++;
        delete u;
      }else{
        u->x++;
      }
    }else{
      u->cnt--;
      if(u->next && u->next->x == u->x + 1){
        small[i] = u->next;
        small[i]->cnt++;
      }else{
        small[i] = new node(u->x + 1, 1);
        small[i]->next = u->next;
        if(u->next) u->next->prev = small[i];
        u->next = small[i];
        small[i]->prev = u;
      }
    }
  }
  // Ai <- Ai - 1  O(1)
  void sub_one(int i){
    auto u = small[i];
    if(u->x == large) large--;
    if(u->cnt == 1){
      if(u->prev && u->prev->x == u->x - 1){
        u->prev->next = u->next;
        if(u->next) u->next->prev = u->prev;
        small[i] = u->prev;
        small[i]->cnt++;
        delete u;
      }else{
        u->x--;
      }
    }else{
      u->cnt--;
      if(u->prev && u->prev->x == u->x - 1){
        small[i] = u->prev;
        small[i]->cnt++;
      }else{
        small[i] = new node(u->x - 1, 1);
        small[i]->prev = u->prev;
        if(u->prev) u->prev->next = small[i];
        u->prev = small[i];
        small[i]->next = u;
      }
    }
  }
  // min O(1)
  Val query(){
    return large;
  }
};
#endif
