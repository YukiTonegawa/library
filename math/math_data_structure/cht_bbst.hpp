#ifndef _CHT_BBST_H_
#define _CHT_BBST_H_
#include <vector>
#include <cassert>

template<typename Val, typename ValSum>
struct cht_bbst{
  static constexpr Val inf = std::numeric_limits<Val>::max() / 2;
  struct line{
    Val a, b;
    line(Val a, Val b): a(a), b(b){}
    inline Val get(Val x){return a * x + b;}
  };
private:
  // ak + b のl <= k < rにおける和
  // a / 2 * {r * (r - 1) - l * (l - 1)} + b * (r - l)
  static constexpr ValSum line_sum(Val a, Val b, Val l, Val r){
    return a * ((ValSum)r * (r - 1) - (ValSum)l * (l - 1)) / 2 + (ValSum)b * (r - l);
  }
  struct node{
    node *l, *r, *p;
    Val lx, rx, ly, ry;
    line x;
    Val min_sub, max_sub;
    ValSum sum, sum_sub;
    node(Val a, Val b, Val _lx, Val _rx): l(nullptr), r(nullptr), p(nullptr), lx(_lx), rx(_rx), x(a, b){
      ly = x.get(lx);
      ry = x.get(rx);
      min_sub = std::min(ly, ry);
      max_sub = std::max(ly, ry);
      sum = sum_sub = line_sum(a, b, lx, rx + 1);
    }
  };
  Val xlow, xhigh;
  node *root;
  int size(node *v){
    return v ? v->sz : 0;
  }
  void update(node *v){
    v->min_sub = std::min(v->ly, v->ry);
    v->max_sub = std::max(v->ly, v->ry);
    v->sum_sub = v->sum;
    if(v->l){
      v->min_sub = std::min(v->min_sub, v->l->min_sub);
      v->max_sub = std::max(v->max_sub, v->l->max_sub);
      v->sum_sub += v->l->sum_sub;
    }
    if(v->r){
      v->min_sub = std::min(v->min_sub, v->r->min_sub);
      v->max_sub = std::max(v->max_sub, v->r->max_sub);
      v->sum_sub += v->r->sum_sub;
    }
  }
  void rotate_right(node *v){
    node *p = v->p, *pp = p->p;
    if((p->l = v->r)) v->r->p = p;
    v->r = p, p->p = v;
    update(p), update(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      update(pp);
    }
  }
  void rotate_left(node *v){
    node *p = v->p, *pp = p->p;
    if((p->r = v->l)) v->l->p = p;
    v->l = p, p->p = v;
    update(p), update(v);
    if((v->p = pp)){
      if(pp->l == p) pp->l = v;
      if(pp->r == p) pp->r = v;
      update(pp);
    }
  }
  void splay(node *v){
    while(v->p){
      node *p = v->p;
      if(!p->p){
        if(p->l == v) rotate_right(v);
        else rotate_left(v);
      }else{
        node *pp = p->p;
        if(pp->l == p){
          if(p->l == v) rotate_right(p);
          else rotate_left(v);
          rotate_right(v);
        }else{
          if(p->r == v) rotate_left(p);
          else rotate_right(v);
          rotate_left(v);
        }
      }
    }
  }
  node *shrink_left(node *v, Val a, Val b){
    node *ret = nullptr;
    while(v){
      // 傾きが同じ直線がある場合
      if(v->x.a == a){
        // 何もせず終了
        if(v->x.b <= b) return v;
        v = v->l;
        continue;
      }
      // ceil(交点)
       Val x = (b - v->x.b + v->x.a - a - 1) / (v->x.a - a);
       if(v->lx < x){
        if(x <= v->rx){
          v->rx = x - 1;
          v->ry = v->x.get(v->rx);
          v->sum = line_sum(v->x.a, v->x.b, v->lx, v->rx + 1);
          update(v);
          return v;
        }
        ret = v;
        v = v->r;
       }else{
        v = v->l;
       }
    }
    return ret;
  }
  node *shrink_right(node *v, Val a, Val b){
    node *ret = nullptr;
    while(v){
      Val x = (v->x.b - b + a - v->x.a - 1) / (a - v->x.a);
      if(x <= v->rx){
        if(v->lx <= x){
          v->lx = x;
          v->ly = v->x.get(v->lx);
          v->sum = line_sum(v->x.a, v->x.b, v->lx, v->rx + 1);
          update(v);
          return v;
        }
        ret = v;
        v = v->l;
      }else{
        v = v->r;
      }
    }
    return ret;
  }
  void add_line_inner(Val a, Val b){
    node *v = root, *u = nullptr, *tmp;
    while(v){
      tmp = v;
      if(v->x.a > a){
        u = v;
        v = v->r;
      }else if(v->x.a < a){
        v = v->l;
      }else{
        u = v;
        break;
      }
    }
    splay(tmp); // 最後に訪れたノードでsplay

    // u := 傾きがa以下の直線のノードの根
    // w := 傾きがaより大きい直線のノードの根
    // u側とw側でsplit
    node *w = nullptr;
    if(u){
      splay(u);
      w = u->r;
      u->r = nullptr;
      if(w) w->p = nullptr;
    }else{
      w = tmp;
    }
    Val l = xlow, r = xhigh;
    // uおよび左側の部分木は残る, 右側は消える
    if(u){
      u = shrink_left(u, a, b);
      if(u){
        splay(u);
        l = max(l, u->rx + 1);
        u->r = nullptr;
        update(u);
      }
    }
    // wおよび右側の部分木は残る, 左側は消える
    if(w){
      w = shrink_right(w, a, b);
      if(w){
        splay(w);
        r = std::min(r, w->lx - 1);
        w->l = nullptr;
        update(w);
      }
    }
    if(l > r){
      if(u){
        if((u->r = w)) w->p = u;
        update(u);
        root = u;
      }else{
        if((w->l = u)) u->p = w;
        update(w);
        root = w;
      }
      return;
    }
    v = new node(a, b, l, r);
    // root <- {u, v, w}
    if((v->l = u)) u->p = v;
    if((v->r = w)) w->p = v;
    update(v);
    root = v;
  }
  line query_inner(Val x){
    node *v = root;
    while(true){
      if(x < v->lx) v = v->l;
      else if(x > v->rx) v = v->r;
      else{
        splay(v);
        root = v;
        return v->x;
      }
    }
  }
  std::pair<Val, Val> find_range_inner(Val x){
    node *v = root;
    while(true){
      if(x < v->lx) v = v->l;
      else if(x > v->rx) v = v->r;
      else{
        splay(v);
        root = v;
        return std::make_pair(v->lx, v->rx + 1);
      }
    }
  }
  // {左側, [l, r)に完全に含まれる区間, 右側}
  std::tuple<node*, node*, node*> split_range(Val l, Val r){
    node *v = root, *left = nullptr, *tmp;
    while(v){
      tmp = v;
      if(v->lx < l){
        left = v;
        v = v->r;
      }else{
        v = v->l;
      }
    }
    splay(tmp); // 最後に訪れたノードでsplay
    if(left){
      splay(left);
      v = left->r;
      left->r = nullptr;
      if(v) v->p = nullptr;
      update(left);
    }else v = tmp;

    node *right = nullptr;
    while(v){
      tmp = v;
      if(r <= v->rx){
        right = v;
        v = v->l;
      }else{
        v = v->r;
      }
    }
    splay(tmp);
    if(right){
      splay(right);
      v = right->l;
      right->l = nullptr;
      if(v) v->p = nullptr;
      update(right);
    }else v = tmp;
    return {left, v, right};
  }
  void merge(node *left, node *v, node *right){
    if(left){
      if((left->r = v)) v->p = left;
      update(left);
      v = left;
    }
    if(right){
      if((right->l = v)) v->p = right;
      update(right);
      v = right;
    }
    root = v;
  }

public:
  // [xlow, xhigh]
  cht_bbst(Val xlow, Val xhigh): xlow(xlow), xhigh(xhigh), root(new node(0, inf, xlow, xhigh)){}
  
  void add_line(Val a, Val b){
    add_line_inner(a, b);
  }
  // xでの最小値
  Val min(Val x){
    assert(xlow <= x && x <= xhigh);
    return query_inner(x).get(x);
  }
  // xでの最小値をとる直線
  line min2(Val x){
    assert(xlow <= x && x <= xhigh);
    return query_inner(x);
  }
  // 座標xで最小値をとる直線が最小値である範囲[l, r)(整数値)
  std::pair<Val, Val> find_range(Val x){
    assert(xlow <= x && x <= xhigh);
    return find_range_inner(x);
  }
  // [l, r)の最小値の最小値
  Val range_min(Val l, Val r){
    assert(xlow <= l && r <= xhigh + 1);
    auto [a, b, c] = split_range(l, r);
    Val ret = b ? b->min_sub : inf;
    if(a && l <= a->rx){
      ret = std::min(ret, std::min(a->ry, a->x.get(l)));
    }
    if(c && c->lx < r){
      ret = std::min(ret, std::min(c->ly, c->x.get(r - 1)));
    }
    merge(a, b, c);
    return ret;
  }
  // [l, r)の最小値の最大値　
  Val range_max(Val l, Val r){
    assert(xlow <= l && r <= xhigh + 1);
    auto [a, b, c] = split_range(l, r);
    Val ret = b ? b->max_sub : -inf;
    if(a && l <= a->rx){
      ret = std::max(ret, std::max(a->ry, a->x.get(l)));
    }
    if(c && c->lx < r){
      ret = std::max(ret, std::max(c->ly, c->x.get(r - 1)));
    }
    merge(a, b, c);
    return ret;
  }
  // [l, r)の最小値の和(整数座標におけるy座標の和)
  ValSum range_sum(Val l, Val r){
    assert(xlow <= l && r <= xhigh + 1);
    auto [a, b, c] = split_range(l, r);
    ValSum ret = b ? b->sum_sub : 0;
    if(a && l <= a->rx){
      ret += line_sum(a->x.a, a->x.b, l, a->rx + 1);
    }
    if(c && c->lx < r){
      ret += line_sum(c->x.a, c->x.b, c->lx, r);
    }
    merge(a, b, c);
    return ret;
  }
};
#endif