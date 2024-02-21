/*
#include ".lib/template.hpp"
#include ".lib/math/mod.hpp"
#include ".lib/fast_io.hpp"

template<typename Idx, typename Val>
struct offline_static_slope_add_rectangle_sum{
private:
  struct bit6{
    int N;
    std::vector<std::array<Val, 6>> sum;
    bit6(int N): N(N), sum(N + 1, {0, 0, 0, 0, 0, 0}){}
    void update(Idx l, Idx r, int lc, int rc, Val z1, Val z2, Val z3){
      Val a = l * z1, b = r * z1, c = l * z2, d = r * z2, e = l * z3, f = r * z3;
      for(int i = lc + 1; i <= N; i += (i & (-i))){
        sum[i][0] -= a;
        sum[i][1] += z1;
        sum[i][2] -= c;
        sum[i][3] += z2;
        sum[i][4] -= e;
        sum[i][5] += z3;
      }
      for(int i = rc + 1; i <= N; i += (i & (-i))){
        sum[i][0] += b;
        sum[i][1] -= z1;
        sum[i][2] += d;
        sum[i][3] -= z2;
        sum[i][4] += f;
        sum[i][5] -= z3;
      }
    }
    std::tuple<Val, Val, Val> query(Idx r, int rc){
      Val a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
      for(int i = rc; i > 0; i -= (i & (-i))){
        a += sum[i][0];
        b += sum[i][1];
        c += sum[i][2];
        d += sum[i][3];
        e += sum[i][4];
        f += sum[i][5];
      }
      return {a + (b * r), c + (d * r), e + (f * r)};
    }
    std::tuple<Val, Val, Val> query(Idx l, Idx r, int lc, int rc){
      auto [cr, dxr, dx2r] = query(r, rc);
      auto [cl, dxl, dx2l] = query(l, lc);
      return {cr - cl, dxr - dxl, dx2r - dx2l};
    }
  };
  struct Update{
    Idx lx, rx, ly, ry;
    Val z1, z2;
    int lyc, ryc;
    Update(Idx lx, Idx rx, Idx ly, Idx ry, Val z1, Val z2, int lyc = 0, int ryc = 0): lx(lx), rx(rx), ly(ly), ry(ry), z1(z1), z2(z2), lyc(lyc), ryc(ryc){}
  };
  struct Query{
    Idx lx, rx, ly, ry;
    int id, lyc, ryc;
    Query(Idx lx, Idx rx, Idx ly, Idx ry, int id, int lyc = 0, int ryc = 0): lx(lx), rx(rx), ly(ly), ry(ry), id(id), lyc(lyc), ryc(ryc){}
  };
  std::vector<Update> U;
  std::vector<Query> Q;

  void solve(std::vector<Val> &ans){
    int N = U.size(), M = Q.size();
    std::vector<Idx> Y;
    for(int i = 0; i < N; i++){
      Y.push_back(U[i].ly);
      Y.push_back(U[i].ry);
    }
    std::sort(Y.begin(), Y.end());
    Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
    auto lb = [&](Idx y) -> int { return std::lower_bound(Y.begin(), Y.end(), y) - Y.begin();};
    for(int i = 0; i < N; i++){      
      int lyc = lb(U[i].ly), ryc = lb(U[i].ry);
      U[i].lyc = lyc, U[i].ryc = ryc;
      U.push_back(Update(U[i].rx, 0, U[i].ly, U[i].ry, -U[i].z1, -U[i].z2, lyc, ryc));
    }
    for(int i = 0; i < M; i++){
      int lyc = lb(Q[i].ly), ryc = lb(Q[i].ry);
      Q[i].lyc = lyc, Q[i].ryc = ryc;
      Q.push_back(Query(Q[i].rx, 0, Q[i].ly, Q[i].ry, Q[i].id + M, lyc, ryc));
    }
    std::sort(U.begin(), U.end(), [](const Update &a, const Update &b){return a.lx < b.lx;});
    std::sort(Q.begin(), Q.end(), [](const Query &a, const Query &b){return a.lx < b.lx;});
    assert(U.size() == 2 * N && Q.size() == 2 * M);
    bit6 bit(Y.size());
    int uid = 0, qid = 0;
    Val i2 = (Val)1 / 2; // use modint
    while(qid < 2 * M){
      if(uid < 2 * N && U[uid].lx < Q[qid].lx){
        Val a = -U[uid].z1 * U[uid].lx - U[uid].z2 * U[uid].lx * (U[uid].lx + 1) * i2;
        Val b = U[uid].z1 + U[uid].z2 * i2;
        Val c = U[uid].z2 * i2;
        bit.update(U[uid].ly, U[uid].ry, U[uid].lyc, U[uid].ryc, a, b, c);
        uid++;
      }else{
        auto [a, b, c] = bit.query(Q[qid].ly, Q[qid].ry, Q[qid].lyc, Q[qid].ryc);
        int id = Q[qid].id;
        if(id >= M){
          ans[id - M] += a + Q[qid].lx * b + (Val)Q[qid].lx * Q[qid].lx * c;
        }else{
          ans[id    ] -= a + Q[qid].lx * b + (Val)Q[qid].lx * Q[qid].lx * c;
        }
        qid++;
      }
    }
  }
public:
  offline_static_slope_add_rectangle_sum(){}
  // [lx, rx) × [ly, ry)にax + bを足す
  void update(Idx lx, Idx rx, Idx ly, Idx ry, Val a, Val b){
    U.push_back(Update(lx, rx, ly, ry, b, a));
  }
  // sum([lx, rx) × [ly, ry))
  void query(Idx lx, Idx rx, Idx ly, Idx ry){
    Q.push_back(Query(lx, rx, ly, ry, Q.size()));
  }
  std::vector<Val> solve(){
    std::vector<Val> ans(Q.size(), 0);
    solve(ans);
    return ans;
  }
};

using mint = modint998244353;

int main(){
  int n, q;
  n = io.in();
  q = io.in();
  offline_static_slope_add_rectangle_sum<int, mint> rect;
  range(i, 0, n){
    int l, d, r, u, w;
    l = io.in();
    d = io.in();
    r = io.in();
    u = io.in();
    w = io.in();
    rect.update(l, r, d, u, 0, w);
  }
  range(i, 0, q){
    int l, d, r, u;
    l = io.in();
    d = io.in();
    r = io.in();
    u = io.in();
    rect.query(l, r, d, u);
  }
  auto ans = rect.solve();
  for(int i = 0; i < ans.size(); i++) io.out(ans[i].val(), '\n');
}
*/
