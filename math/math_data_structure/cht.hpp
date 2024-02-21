#ifndef _CHT_H_
#define _CHT_H_
#include <algorithm>
#include <vector>
/*
制限
  x_high + 1がオーバーフローしない
  x_low + x_highがオーバーフローしない
  f(x_high), f(x_low)がオーバーフローしない
*/
template<typename Val>
struct lichao_tree{
  static constexpr Val inf = std::numeric_limits<Val>::max();
  Val x_low, x_high;
  struct line{
    Val a, b;
    line(Val a, Val b): a(a), b(b){}
    inline Val get(Val x){return a * x + b;}
    bool operator != (const line &r){return a != r.a || b != r.b;}
  };
private:
  struct node{
    line x;
    node *l, *r;
    node(const line &x) : x(x), l(nullptr), r(nullptr){}
  };
  node *root;
  node *add_line(node *v, line &x, Val l, Val r, Val x_l, Val x_r){
    if(!v) return new node(x);
    Val t_l = v->x.get(l), t_r = v->x.get(r);
    if(l + 1 == r){
      if(x_l < t_l) v->x = x;
      return v;
    }else if(t_l <= x_l && t_r <= x_r){
      return v;
    }else if(t_l >= x_l && t_r >= x_r){
      v->x = x;
      return v;
    }else{
      Val m = (l + r) / 2;
      Val t_m = v->x.get(m), x_m = x.get(m);
      if(t_m > x_m){
        std::swap(v->x, x);
        if(x_l >= t_l) v->l = add_line(v->l, x, l, m, t_l, t_m);
        else v->r = add_line(v->r, x, m, r, t_m, t_r);
      }else{
        if(t_l >= x_l) v->l = add_line(v->l, x, l, m, x_l, x_m);
        else v->r = add_line(v->r, x, m, r, x_m, x_r);
      }
      return v;
    }
  }
  node *add_segment(node *v, line &x, Val a, Val b, Val l, Val r, Val x_l, Val x_r){
    if(r <= a || b <= l) return v;
    if(a <= l && r <= b){
      line y(x);
      return add_line(v, y, l, r, x_l, x_r);
    }
    if(v){
      Val t_l = v->x.get(l), t_r = v->x.get(r);
      if(t_l <= x_l && t_r <= x_r) return v;
    }else{
      v = new node(line(0, inf));
    }
    Val m = (l + r) / 2;
    Val x_m = x.get(m);
    v->l = add_segment(v->l, x, a, b, l, m, x_l, x_m);
    v->r = add_segment(v->r, x, a, b, m, r, x_m, x_r);
    return v;
  }
  Val min(node *v, Val l, Val r, Val x){
    if(!v) return inf;
    if(l + 1 == r) return v->x.get(x);
    Val m = (l + r) / 2;
    if(x < m) return std::min(v->x.get(x), min(v->l, l, m, x));
    else return std::min(v->x.get(x), min(v->r, m, r, x));
  }
  line min2(node *v, Val l, Val r, Val x){
    if(!v) return line(0, inf);
    if(l + 1 == r) return v->x;
    Val m = (l + r) / 2;
    if(x < m){
      line res = min2(v->l, l, m, x);
      return v->x.get(x) <= res.get(x) ? v->x : res;
    }else{
      line res = min2(v->r, m, r, x);
      return v->x.get(x) <= res.get(x) ? v->x : res;
    }
  }
  // 定数倍高速化できそうだが, addがボトルネックになりがちなのでとりあえずこれで
  void enumerate(node *v, Val l, Val r, std::vector<line> &L, std::vector<std::pair<Val, line>> &res){
    L.push_back(v->x);
    Val m = (l + r) / 2;

    auto calc = [&](Val lx, Val rx)->void{
      while(true){
        line a = min2(lx);
        if(res.empty() || res.back().second != a) res.push_back({lx, a});
        if(a.b == inf) return;
        Val intersection = rx; // min(ceil(交点))
        line next_line = line(0, inf);
        for(int i = 0; i < L.size(); i++){
          if(L[i].b == inf || L[i].a >= a.a) continue;
          // (a.a - L[i].a) x = L[i].b - a.b
          Val diff_a = a.a - L[i].a;
          Val ceil_x = (L[i].b - a.b + diff_a - 1) / diff_a;
          if(lx < ceil_x && ceil_x < intersection){
            assert(lx < ceil_x);
            next_line = L[i];
            intersection = ceil_x;
          }
        }
        if(intersection == rx) break;
        lx = intersection;
      }
    };
    if(v->l && v->r){
      enumerate(v->l, l, m, L, res);
      enumerate(v->r, m, r, L, res);
    }else if(!v->l && !v->r){
      calc(l, r);
    }else if(!v->l){
      calc(l, m);
      enumerate(v->r, m, r, L, res);
    }else{
      enumerate(v->l, l, m, L, res);
      calc(m, r);
    }
    L.pop_back();
  }
public:
  lichao_tree(Val x_low, Val x_high): x_low(x_low), x_high(x_high), root(nullptr){}
  // 直線ax + bを追加
  void add_line(Val a, Val b){
    line x(a, b);
    root = add_line(root, x, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
  }
  // 線分ax + b, [l, r)を追加
  void add_segment(Val l, Val r, Val a, Val b){
    line x(a, b);
    root = add_segment(root, x, l, r, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
  }
  // f(x)の最小値
  Val min(Val x){
    return min(root, x_low, x_high + 1, x);
  }
  // f(x)の最小値をとる線
  line min2(Val x){
    return min2(root, x_low, x_high + 1, x);
  }
  // {l, f}, fはl, (次のl)において最小値を取る
  std::vector<std::pair<Val, line>> enumerate(){
    if(!root) return {{x_low, line(0, inf)}};
    std::vector<std::pair<Val, line>> res;
    std::vector<line> L;
    enumerate(root, x_low, x_high + 1, L, res);
    return res;
  }
};

template<typename Val>
struct compressed_lichao_tree{
  static constexpr Val inf = std::numeric_limits<Val>::max();
  struct line{
    Val a, b;
    line(Val a, Val b): a(a), b(b){}
    inline Val get(Val x){return a * x + b;}
  };
private:
  int N, M;
  Val x_high;
  std::vector<line> L;
  std::vector<Val> sample_points;

  void add_line(int k, line &x, int lx, int rx, Val x_l, Val x_r){
    if(k >= 2 * M - 1) return;
    Val l = sample_points[lx], r = sample_points[rx];
    Val t_l = L[k].get(l), t_r = L[k].get(r);
    if(M - 1 <= k){
      if(x_l < t_l) L[k] = x;
      return;
    }else if(t_l <= x_l && t_r <= x_r){
      return;
    }else if(t_l >= x_l && t_r >= x_r){
      L[k] = x;
      return;
    }else{
      int mx = (lx + rx) / 2;
      Val m = sample_points[mx];
      Val t_m = L[k].get(m), x_m = x.get(m);
      if(t_m > x_m){
        std::swap(L[k], x);
        if(x_l >= t_l) add_line(k * 2 + 1, x, lx, mx, t_l, t_m);
        else add_line(k * 2 + 2, x, mx, rx, t_m, t_r);
      }else{
        if(t_l >= x_l) add_line(k * 2 + 1, x, lx, mx, x_l, x_m);
        else add_line(k * 2 + 2, x, mx, rx, x_m, x_r);
      }
    }
  }
  void add_segment(int k, line &x, int ax, int bx, int lx, int rx, Val x_l, Val x_r){
    if(2 * M - 1 <= k || rx < ax || bx < lx) return;
    if(ax <= lx && rx <= bx){
      line y(x);
      return add_line(k, y, lx, rx, x_l, x_r);
    }
    if(L[k].get(sample_points[lx]) <= x_l && L[k].get(sample_points[rx]) <= x_r) return;
    int mx = (lx + rx) / 2;
    Val m = sample_points[mx];
    Val x_m = x.get(m);
    add_segment(k * 2 + 1, x, ax, bx, lx, mx, x_l, x_m);
    add_segment(k * 2 + 2, x, ax, bx, mx, rx, x_m, x_r);
  }
  Val min(int k, int lx, int rx, Val x){
    if(M - 1 <= k) return L[k].get(x);
    int mx = (lx + rx) / 2;
    Val m = sample_points[mx];
    if(x < m) return std::min(L[k].get(x), min(k * 2 + 1, lx, mx, x));
    else return std::min(L[k].get(x), min(k * 2 + 2, mx, rx, x));
  }
  line min2(int k, Val lx, Val rx, Val x){
    if(M - 1 <= k) return L[k];
    int mx = (lx + rx) / 2;
    Val m = sample_points[mx];
    if(x < m){
      line res = min2(k * 2 + 1, lx, mx, x);
      return L[k].get(x) <= res.get(x) ? L[k] : res;
    }else{
      line res = min2(k * 2 + 2, mx, rx, x);
      return L[k].get(x) <= res.get(x) ? L[k] : res;
    }
  }
  static constexpr int ceil_pow2(int n){
    int m = 1;
    while(m < n) m <<= 1;
    return m;
  }
public:
  // x: クエリを飛ばす点の集合(線分の端点は追加しなくていい)
  // xはソート済みかつユニーク
  compressed_lichao_tree(const std::vector<Val> &x):
  N(x.size()), M(ceil_pow2(N)), L(2 * M - 1, line(0, inf)), sample_points(x){
    x_high = (x.empty() ? 0 : x.back() + 1);
    sample_points.resize(M + 1, x_high);
  }
  // 直線ax + bを追加
  void add_line(Val a, Val b){
    line x(a, b);
    add_line(0, x, 0, M, x.get(sample_points[0]), x.get(x_high));
  }
  // 線分ax + b, [l, r)を追加
  void add_segment(Val l, Val r, Val a, Val b){
    line x(a, b);
    int lx = std::lower_bound(sample_points.begin(), sample_points.end(), l) - sample_points.begin();
    int rx = std::lower_bound(sample_points.begin(), sample_points.end(), r) - sample_points.begin();
    add_segment(0, x, lx, rx, 0, M, x.get(sample_points[0]), x.get(x_high));
  }
  // f(x)の最小値
  Val min(Val x){
    return min(0, 0, M, x);
  }
  // f(x)の最小値をとる線
  line min2(Val x){
    return min2(0, 0, M, x);
  }
};

template<typename Val>
struct lichao_tree_range_min{
  static constexpr Val inf = std::numeric_limits<Val>::max();
  Val x_low, x_high;
  struct line{
    Val a, b;
    line(Val a, Val b): a(a), b(b){}
    inline Val get(Val x){return a * x + b;}
    bool operator != (const line &r){return a != r.a || b != r.b;}
  };
private:
  struct node{
    line x;
    node *l, *r;
    Val min;
    node(const line &x, Val m) : x(x), l(nullptr), r(nullptr), min(m){}
  };
  node *root;
  node *add_line(node *v, line &x, Val l, Val r, Val x_l, Val x_r){
    if(!v) return new node(x, std::min(x_l, x_r));
    Val t_l = v->x.get(l), t_r = v->x.get(r);
    if(l + 1 == r){
      if(x_l < t_l) v->x = x, v->min = x_l;
      return v;
    }else if(t_l <= x_l && t_r <= x_r){
      return v;
    }else if(t_l >= x_l && t_r >= x_r){
      v->x = x;
      v->min = std::min(v->min, std::min(x_l, x_r));
      return v;
    }else{
      Val m = (l + r) / 2;
      if(m == r) --m;
      Val t_m = v->x.get(m), x_m = x.get(m);
      if(t_m > x_m){
        std::swap(v->x, x);
        v->min = std::min(v->min, std::min(x_l, x_r));
        if(x_l >= t_l) v->l = add_line(v->l, x, l, m, t_l, t_m);
        else v->r = add_line(v->r, x, m + 1, r, t_m + x.a, t_r);
      }else{
        if(t_l >= x_l) v->l = add_line(v->l, x, l, m, x_l, x_m);
        else v->r = add_line(v->r, x, m + 1, r, x_m + x.a, x_r);
      }
      if(v->l && v->min > v->l->min) v->min = v->l->min;
      if(v->r && v->min > v->r->min) v->min = v->r->min;
      return v;
    }
  }
  node *add_segment(node *v, line &x, Val a, Val b, Val l, Val r, Val x_l, Val x_r){
    if(r < a || b < l) return v;
    if(a <= l && r <= b){
      line y(x);
      return add_line(v, y, l, r, x_l, x_r);
    }
    if(v){
      Val t_l = v->x.get(l), t_r = v->x.get(r);
      if(t_l <= x_l && t_r <= x_r) return v;
    }else{
      v = new node(line(0, inf), inf);
    }
    Val m = (l + r) / 2;
    if(m == r) --m;
    Val x_m = x.get(m);
    v->l = add_segment(v->l, x, a, b, l, m, x_l, x_m);
    v->r = add_segment(v->r, x, a, b, m + 1, r, x_m + x.a, x_r);
    if(v->l && v->min > v->l->min) v->min = v->l->min;
    if(v->r && v->min > v->r->min) v->min = v->r->min;
    return v;
  }
  Val min(node *v, Val x, Val l, Val r){
    if(!v) return inf;
    if(l == r) return v->x.get(x);
    Val m = (l + r) / 2;
    if(m == r) --m;
    if(x <= m) return std::min(v->x.get(x), min(v->l, x, l, m));
    else return std::min(v->x.get(x), min(v->r, x, m + 1, r));
  }
  line min2(node *v, Val x, Val l, Val r){
    if(!v) return line(0, inf);
    if(l == r) return v->x;
    Val m = (l + r) / 2;
    if(m == r) --m;
    if(x <= m){
      line res = min2(v->l, x, l, m);
      return v->x.get(x) <= res.get(x) ? v->x : res;
    }else{
      line res = min2(v->r, x, m + 1, r);
      return v->x.get(x) <= res.get(x) ? v->x : res;
    }
  }
  Val range_min(node *v, Val a, Val b, Val l, Val r){
    if(!v || r < a || b < l) return inf;
    if(a <= l && r <= b){
      return v->min;
    }
    Val m = (l + r) / 2;
    if(m == r) --m;
    Val min1 = std::min(v->x.get(std::max(a, l)), v->x.get(std::min(b, r)));
    Val min2 = std::min(range_min(v->l, a, b, l, m), range_min(v->r, a, b, m + 1, r));
    return std::min(min1, min2);
  }
  // 定数倍高速化できそうだが, addがボトルネックになりがちなのでとりあえずこれで
  void enumerate(node *v, Val l, Val r, std::vector<line> &L, std::vector<std::pair<Val, line>> &res){
    L.push_back(v->x);
    Val m = (l + r) / 2;
    if(m == r) --m;

    // [lx, rx]
    auto calc = [&](Val lx, Val rx)->void{
      while(true){
        line a = min2(lx);
        if(res.empty() || res.back().second != a) res.push_back({lx, a});
        if(a.b == inf) return;
        Val intersection = rx + 1; // min(ceil(交点))
        line next_line = line(0, inf);
        for(int i = 0; i < L.size(); i++){
          if(L[i].b == inf || L[i].a >= a.a) continue;
          // (a.a - L[i].a) x = L[i].b - a.b
          Val diff_a = a.a - L[i].a;
          assert(diff_a != 0);
          Val ceil_x = (L[i].b - a.b + diff_a - 1) / diff_a;
          if(lx < ceil_x && ceil_x < intersection){
            assert(lx < ceil_x);
            next_line = L[i];
            intersection = ceil_x;
          }
        }
        if(intersection > rx) break;
        lx = intersection;
      }
    };
    if(v->l && v->r){
      enumerate(v->l, l, m, L, res);
      enumerate(v->r, m + 1, r, L, res);
    }else if(!v->l && !v->r){
      calc(l, r);
    }else if(!v->l){
      calc(l, m);
      enumerate(v->r, m + 1, r, L, res);
    }else{
      enumerate(v->l, l, m, L, res);
      calc(m + 1, r);
    }
    L.pop_back();
  }
public:
  lichao_tree_range_min(Val x_low, Val x_high): x_low(x_low), x_high(x_high), root(nullptr){}
  // 直線ax + bを追加
  void add_line(Val a, Val b){
    line x(a, b);
    root = add_line(root, x, x_low, x_high, x.get(x_low), x.get(x_high));
  }
  // 線分ax + b, [l, r)を追加
  void add_segment(Val l, Val r, Val a, Val b){
    assert(l <= r);
    line x(a, b);
    root = add_segment(root, x, l, r - 1, x_low, x_high, x.get(x_low), x.get(x_high));
  }
  // f(x)の最小値
  Val min(Val x){
    return min(root, x, x_low, x_high);
  }
  // f(x)の最小値をとる線
  line min2(Val x){
    return min2(root, x, x_low, x_high);
  }
  // l <= x < rを満たす整数xでのf(x)の最小値
  Val range_min(Val l, Val r){
    return range_min(root, l, r - 1, x_low, x_high);
  }
  // {l, f}, fはl, (次のl)において最小値を取る
  std::vector<std::pair<Val, line>> enumerate(){
    if(!root) return {{x_low, line(0, inf)}};
    std::vector<std::pair<Val, line>> res;
    std::vector<line> L;
    enumerate(root, x_low, x_high, L, res);
    return res;
  }
};
template<typename Val>
struct lichao_tree_range_add{
  static constexpr Val inf = std::numeric_limits<Val>::max();
  Val x_low, x_high;
  struct line{
    Val a, b;
    line(Val a, Val b): a(a), b(b){}
    inline Val get(Val x){return a * x + b;}
    bool operator != (const line &r){return a != r.a || b != r.b;}
  };
private:
  struct node{
    line x;
    node *l, *r;
    Val lazy;
    node(const line &x) : x(x), l(nullptr), r(nullptr), lazy(0){}
  };
  node *root;
  
  void push_lazy(node *v){
    if(v->lazy == 0) return;
    if(v->l) propagate(v->l, v->lazy);
    if(v->r) propagate(v->r, v->lazy);
    v->lazy = 0;
  }
  void propagate(node *v, Val lazy){
    v->lazy += lazy;
    if(v->x.b != inf) v->x.b += lazy;
  }
  // vの直線をv->l, v->rにpushする
  void push_line(node *v, Val l, Val r){
    Val m = (l + r) / 2;
    if(m == r) --m;
    Val xl = v->x.get(l), xm = v->x.get(m), xr = v->x.get(r);
    v->l = add_line(v->l, v->x, l, m, xl, xm);
    v->r = add_line(v->r, v->x, m + 1, r, xm + v->x.a, xr);
    v->x = {0, inf};
  }
  node *add_line(node *v, line &x, Val l, Val r, Val x_l, Val x_r){
    if(!v) return new node(x);
    push_lazy(v);
    Val t_l = v->x.get(l), t_r = v->x.get(r);
    if(l + 1 == r){
      if(x_l < t_l) v->x = x;
      return v;
    }else if(t_l <= x_l && t_r <= x_r){
      return v;
    }else if(t_l >= x_l && t_r >= x_r){
      v->x = x;
      return v;
    }else{
      Val m = (l + r) / 2;
      if(m == r) --m;
      Val t_m = v->x.get(m), x_m = x.get(m);
      if(t_m > x_m){
        std::swap(v->x, x);
        if(x_l >= t_l) v->l = add_line(v->l, x, l, m, t_l, t_m);
        else v->r = add_line(v->r, x, m + 1, r, t_m + x.a, t_r);
      }else{
        if(t_l >= x_l) v->l = add_line(v->l, x, l, m, x_l, x_m);
        else v->r = add_line(v->r, x, m + 1, r, x_m + x.a, x_r);
      }
      return v;
    }
  }
  node *add_segment(node *v, line &x, Val a, Val b, Val l, Val r, Val x_l, Val x_r){
    if(r < a || b < l) return v;
    if(a <= l && r <= b){
      line y(x);
      return add_line(v, y, l, r, x_l, x_r);
    }
    if(v){
      Val t_l = v->x.get(l), t_r = v->x.get(r);
      if(t_l <= x_l && t_r <= x_r) return v;
      push_lazy(v);
    }else{
      v = new node(line(0, inf));
    }
    Val m = (l + r) / 2;
    if(m == r) --m;
    Val x_m = x.get(m);
    v->l = add_segment(v->l, x, a, b, l, m, x_l, x_m);
    v->r = add_segment(v->r, x, a, b, m + 1, r, x_m + x.a, x_r);
    return v;
  }
  void add(node *v, Val x, Val a, Val b, Val l, Val r){
    if(!v || r < a || b < l) return;
    if(a <= l && r <= b){
      propagate(v, x);
      return;
    }
    push_lazy(v);
    push_line(v, l, r);
    Val m = (l + r) / 2;
    if(m == r) --m;
    add(v->l, x, a, b, l, m);
    add(v->r, x, a, b, m + 1, r);
  }
  Val min(node *v, Val x, Val l, Val r){
    if(!v) return inf;
    if(l == r) return v->x.get(x);
    push_lazy(v);
    Val m = (l + r) / 2;
    if(m == r) --m;
    if(x <= m) return std::min(v->x.get(x), min(v->l, x, l, m));
    else return std::min(v->x.get(x), min(v->r, x, m + 1, r));
  }
  line min2(node *v, Val x, Val l, Val r){
    if(!v) return line(0, inf);
    if(l == r) return v->x;
    push_lazy(v);
    Val m = (l + r) / 2;
    if(m == r) --m;
    if(x <= m){
      line res = min2(v->l, x, l, m);
      return v->x.get(x) <= res.get(x) ? v->x : res;
    }else{
      line res = min2(v->r, x, m + 1, r);
      return v->x.get(x) <= res.get(x) ? v->x : res;
    }
  }
public:
  lichao_tree_range_add(Val x_low, Val x_high): x_low(x_low), x_high(x_high), root(nullptr){}
  // 直線ax + bを追加
  void add_line(Val a, Val b){
    line x(a, b);
    root = add_line(root, x, x_low, x_high, x.get(x_low), x.get(x_high));
  }
  // 線分ax + b, [l, r)を追加
  void add_segment(Val l, Val r, Val a, Val b){
    assert(l <= r);
    line x(a, b);
    root = add_segment(root, x, l, r - 1, x_low, x_high, x.get(x_low), x.get(x_high));
  }
  // [l, r)にxを足す
  void add(Val l, Val r, Val x){
    add(root, x, l, r - 1, x_low, x_high);
  }
  // f(x)の最小値
  Val min(Val x){
    return min(root, x, x_low, x_high);
  }
  // f(x)の最小値をとる線
  line min2(Val x){
    return min2(root, x, x_low, x_high);
  }
};
/*
template<typename Val>
struct lichao_tree_quadratic{
  static constexpr Val inf = std::numeric_limits<Val>::max();
  Val x_low, x_high;
  struct line{
    Val a, b, c;
    line(Val a, Val b, Val c): a(a), b(b), c(c){}
    inline Val get(Val x){return a * x * x + b * x + c;}
    bool operator != (const line &r){return a != r.a || b != r.b;}
  };
private:
  struct node{
    line x;
    node *l, *r;
    node(const line &x) : x(x), l(nullptr), r(nullptr){}
  };
  node *root;
  node *add_line(node *v, line &x, Val l, Val r, Val x_l, Val x_r){
    if(!v) return new node(x);
    Val t_l = v->x.get(l), t_r = v->x.get(r);
    if(l + 1 == r){
      if(x_l < t_l) v->x = x;
      return v;
    }else if(t_l <= x_l && t_r <= x_r){
      return v;
    }else if(t_l >= x_l && t_r >= x_r){
      v->x = x;
      return v;
    }else{
      Val m = (l + r) / 2;
      Val t_m = v->x.get(m), x_m = x.get(m);
      if(t_m > x_m){
        std::swap(v->x, x);
        if(x_l >= t_l) v->l = add_line(v->l, x, l, m, t_l, t_m);
        else v->r = add_line(v->r, x, m, r, t_m, t_r);
      }else{
        if(t_l >= x_l) v->l = add_line(v->l, x, l, m, x_l, x_m);
        else v->r = add_line(v->r, x, m, r, x_m, x_r);
      }
      return v;
    }
  }
  node *add_segment(node *v, line &x, Val a, Val b, Val l, Val r, Val x_l, Val x_r){
    if(r <= a || b <= l) return v;
    if(a <= l && r <= b){
      line y(x);
      return add_line(v, y, l, r, x_l, x_r);
    }
    if(v){
      Val t_l = v->x.get(l), t_r = v->x.get(r);
      if(t_l <= x_l && t_r <= x_r) return v;
    }else{
      v = new node(line(0, 0, inf));
    }
    Val m = (l + r) / 2;
    Val x_m = x.get(m);
    v->l = add_segment(v->l, x, a, b, l, m, x_l, x_m);
    v->r = add_segment(v->r, x, a, b, m, r, x_m, x_r);
    return v;
  }
  Val min(node *v, Val l, Val r, Val x){
    if(!v) return inf;
    if(l + 1 == r) return v->x.get(x);
    Val m = (l + r) / 2;
    if(x < m) return std::min(v->x.get(x), min(v->l, l, m, x));
    else return std::min(v->x.get(x), min(v->r, m, r, x));
  }
  line min2(node *v, Val l, Val r, Val x){
    if(!v) return line(0, inf);
    if(l + 1 == r) return v->x;
    Val m = (l + r) / 2;
    if(x < m){
      line res = min2(v->l, l, m, x);
      return v->x.get(x) <= res.get(x) ? v->x : res;
    }else{
      line res = min2(v->r, m, r, x);
      return v->x.get(x) <= res.get(x) ? v->x : res;
    }
  }
  // 定数倍高速化できそうだが, addがボトルネックになりがちなのでとりあえずこれで
  void enumerate(node *v, Val l, Val r, std::vector<line> &L, std::vector<std::pair<Val, line>> &res){
    L.push_back(v->x);
    Val m = (l + r) / 2;

    auto calc = [&](Val lx, Val rx)->void{
      while(true){
        line a = min2(lx);
        if(res.empty() || res.back().second != a) res.push_back({lx, a});
        if(a.b == inf) return;
        Val intersection = rx; // min(ceil(交点))
        line next_line = line(0, inf);
        for(int i = 0; i < L.size(); i++){
          if(L[i].b == inf || L[i].a >= a.a) continue;
          // (a.a - L[i].a) x = L[i].b - a.b
          Val diff_a = a.a - L[i].a;
          Val ceil_x = (L[i].b - a.b + diff_a - 1) / diff_a;
          if(lx < ceil_x && ceil_x < intersection){
            assert(lx < ceil_x);
            next_line = L[i];
            intersection = ceil_x;
          }
        }
        if(intersection == rx) break;
        lx = intersection;
      }
    };
    if(v->l && v->r){
      enumerate(v->l, l, m, L, res);
      enumerate(v->r, m, r, L, res);
    }else if(!v->l && !v->r){
      calc(l, r);
    }else if(!v->l){
      calc(l, m);
      enumerate(v->r, m, r, L, res);
    }else{
      enumerate(v->l, l, m, L, res);
      calc(m, r);
    }
    L.pop_back();
  }
public:
  lichao_tree_quadratic(Val x_low, Val x_high): x_low(x_low), x_high(x_high), root(nullptr){}
  // 線分ax^2 + bx + cを追加
  void add_line(Val a, Val b, Val c){
    line x(a, b, c);
    if(a){
      Val mid = -(b / (2 * a));
      if(mid < x_low || x_high < mid){
        root = add_segment(root, x, x_low, x_high + 1, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
      }else{
        root = add_segment(root, x, x_low, mid + 1, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
        root = add_segment(root, x, mid, x_high + 1, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
      }
    }else{
      root = add_segment(root, x, x_low, x_high + 1, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
    }
  }
  // 線分ax^2 + bx + c [l, r)を追加
  void add_segment(Val l, Val r, Val a, Val b, Val c){
    line x(a, b, c);
    if(a){
      Val mid = -(b / (2 * a));
      if(mid < l || r <= mid){
        root = add_segment(root, x, l, r, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
      }else{
        root = add_segment(root, x, l, mid + 1, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
        root = add_segment(root, x, mid, r, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
      }
    }else{
      root = add_segment(root, x, l, r, x_low, x_high + 1, x.get(x_low), x.get(x_high + 1));
    }
  }
  // f(x)の最小値
  Val min(Val x){
    return min(root, x_low, x_high + 1, x);
  }
  // f(x)の最小値をとる線
  line min2(Val x){
    return min2(root, x_low, x_high + 1, x);
  }
  // {l, f}, fはl, (次のl)において最小値を取る
  std::vector<std::pair<Val, line>> enumerate(){
    if(!root) return {{x_low, line(0, inf)}};
    std::vector<std::pair<Val, line>> res;
    std::vector<line> L;
    enumerate(root, x_low, x_high + 1, L, res);
    return res;
  }
};
*/
#endif
