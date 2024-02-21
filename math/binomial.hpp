#ifndef _BINOMIAL_H_
#define _BINOMIAL_H_
#include "mod.hpp"
#include "../data_structure/range_query/mo_algorithm.hpp"
// f(i, j) := comb(i, j) * comb(n - i, m - j), {0 <= i <= n, 0 <= j <= m}
// h(i, j) := ∑{0 <= k <= j} f(i, j) 　<-jについての和　
// f(i, j)は(n - m) * mグリッドにおける点(i - j, j)を通る最短経路{左上(0, 0)と右下(n, m)の最短経路}の数
// f(i, j)とf(i, j + 1)は最短経路の数え上げにおいて独立{(i - j, j)と(i - j - 1, j + 1)を両方通ることが不可能}
// h(i, j)は(i - k, k) {0 <= k <= i}を通る場合の数(ナナメ)である

// 遷移
// h(i, j + 1) = h(i, j) + f(i, j + 1)
// h(i, j - 1) = h(i, j) - f(i, j)
// h(i, j) -> h(i + 1, j) : 点(i - j, j)および点(i - j, j + 1)を通る最短経路の数を引く -comb(i, j) * comb(n - i, m - j)
// h(i, j) -> h(i - 1, j) : 点(i - j - 1, j)および点(i - j - 1, j + 1)を通る最短経路の数を足す +comb(i - 1, j) * comb(n - i + 1, m - j)

/* #を通る最短経路の数
h(3, 1)
s..........
...........
.#.........
#.........t

h(3, 1)の言い換え
s..........
...........
##.........
..........t

s..........
...........
.#.........
.#........t
*/


template<typename mint>
std::vector<mint> diagonal_path_sum(int n, int m, std::vector<std::pair<int, int>> Q){
  struct diagonal_path_sum_struct{
    int n, m;
    modcomb<mint> mcb;
    // f(0, 0) = comb(n, m)
    mint x;

    diagonal_path_sum_struct(int n, int m): n(n), m(m), mcb(n + 1), x(mcb.comb(n, m)){}

    void add_left(int l, int r){
      x -= mcb.comb(r - 1, l + 1) * mcb.comb(n - r + 1, m - l - 1);
    }
    void add_right(int l, int r){
      r--;
      x -= mcb.comb(r, l) * mcb.comb(n - r - 1, m - l - 1);
    }
    void del_left(int l, int r){
      x += mcb.comb(r - 1, l + 1) * mcb.comb(n - r + 1, m - l - 1);
    }
    void del_right(int l, int r){
      r--;
      x += mcb.comb(r, l) * mcb.comb(n - r - 1, m - l - 1);
    }
  };
  diagonal_path_sum_struct dp(n, m);
  mo_algorithm mo(n + 1, dp);
  int q = Q.size();
  std::vector<mint> ans(q);
  for(int i = 0; i < q; i++) mo.insert(Q[i].first, Q[i].second + 1);
  for(int i = 0; i < q; i++) ans[mo.process2().first] = dp.x;
  return ans;
}

/*
// f(i, j) = ∑{0 <= k <= j} nCi * (n - i)Ck
// n個のものを{i, j, n - i - j}の3つに分ける場合の数のjに関する和
// f(i, j) = nCi * (n - i + 1)C(j + 1)
template<typename mint>
mint split_3_row(int n, int i, int j){
  static modcomb<mint> mcb;
  mcb.recalc(n + 1);
  return mcb.comb(n, i) * mcb.comb(n - i + 1, j + 1);
}
*/


// f(i, j) =  iCj + iC(j - 2) + iC(j - 4) ...
template<typename mint>
std::vector<mint> stripe_sum(std::vector<std::pair<int, int>> Q){
  struct stripe_sum_struct{
    modcomb<mint> mcb;
    mint i2;
    // f(0, 0) = 1
    mint x, y, tmp;
    stripe_sum_struct(): i2(mint(2).inv()), x(1), y(0){}

    void add_left(int l, int r){
      tmp = x - mcb.comb(r - 1, l + 1);
      x = y;
      y = tmp;
    }
    void del_left(int l, int r){
      tmp = y + mcb.comb(r - 1, l + 1);
      y = x;
      x = tmp;
    }
    void add_right(int l, int r){
      tmp = x + y;  // ∑{0 <= k <= j}iCk
      x = tmp;
      y = tmp - mcb.comb(r, l);
    }
    void del_right(int l, int r){
      r--;
      tmp = x + y; // ∑{0 <= k <= j}iCk
      tmp = (tmp + mcb.comb(r, l)) * i2; // ∑{0 <= k <= j}(i - 1)Ck
      tmp = (tmp + mcb.comb(r - 1, l)) * i2; // ∑{0 <= k <= j}(i - 2)Ck
      x = tmp;
      y = tmp - mcb.comb(r, l);
    }
  };
  int max_n = 0;
  for(auto [i, j] : Q) max_n = std::max(max_n, i); // nCrの最大のn
  stripe_sum_struct st;
  st.mcb.recalc(max_n + 1);
  mo_algorithm mo(max_n + 1, st);
  int q = Q.size();
  std::vector<mint> ans(q);
  for(int i = 0; i < q; i++){
    if(Q[i].first < Q[i].second) Q[i].second = Q[i].first;
    mo.insert(Q[i].second, Q[i].first + 1);
  }
  for(int i = 0; i < q; i++){
    auto [idx, id] = mo.process2();
    ans[idx] = st.x;
  }
  return ans;
}



// f(i, j) := ∑{0 <= k <= j} comb(i, k)
// パスカルの三角形におけるi行目[0, j]列目の和

// 遷移
// f(i, j + 1) = f(i, j) + comb(i, j + 1)を足す
// f(i, j - 1) = f(i, j) - comb(i, j)を引く
// f(i + 1, j) = f(i, j) + f(i, j - 1)
//   = 2 * f(i, j) - comb(i, j) = f(i + 1, j)
// f(i, j) -> f(i - 1, j)も同様

/*
i < jの場合i = jにする
1
1 1
1 2 1
1 3 3 1
1 4 6 4 1
*/
template<typename mint>
std::vector<mint> pascal_sum_row(std::vector<std::pair<int, int>> Q){
  struct pascal_sum_row_struct{
    modcomb<mint> mcb;
    mint i2 = mint(2).inv();
    mint x = i2;

    pascal_sum_row_struct(){}

    void add_left(int l, int r){
      x -= mcb.comb(r - 1, l + 1);
    }
    void del_left(int l, int r){
      x += mcb.comb(r - 1, l + 1);
    }
    void add_right(int l, int r){
      r--;
      x = 2 * x - mcb.comb(r, l);
    }
    void del_right(int l, int r){
      r--;
      x = (x + mcb.comb(r, l)) * i2;
    }
  };
  int max_n = 0;
  for(auto [i, j] : Q) max_n = std::max(max_n, i); // nCrの最大のn
  pascal_sum_row_struct pas;
  pas.mcb.recalc(max_n + 1);
  mo_algorithm mo(max_n + 1, pas);
  int q = Q.size();
  std::vector<mint> ans(q);
  for(int i = 0; i < q; i++){
    if(Q[i].first < Q[i].second) Q[i].second = Q[i].first;
    mo.insert(Q[i].second, Q[i].first + 1);
  }
  for(int i = 0; i < q; i++){
    auto [idx, id] = mo.process2();
    ans[idx] = pas.x;
  }
  return ans;
}


// 等差数列の重み付き
// g(i, j) := ∑{0 <= k <= j} (ak + b) * comb(i, k)

// 遷移
// g(i, j) -> g(i, j + 1): (a(j + 1) + b) * comb(i, j + 1)を足す
// g(i, j) -> g(i, j - 1): (aj + b) * comb(i, j)を引く
// g(i, j) -> g(i + 1, j):
// 各項(ak + b) * comb(i + 1, k)
// = (ak + b) * (comb(i, k) + comb(i, k - 1))
// = (ak + b) * comb(i, k) + (a(k - 1) + b) * comb(i, k - 1) + a * comb(i, k - 1)
// より
// g(i + 1, j) = g(i, j) + g(i, j - 1) + a * f(i, j - 1)
//             = 2 * g(i, j) - (aj + b) * comb(i, j) + a * f(i, j) - a * comb(i, j)
//             = 2 * g(i, j) + a * f(i, j) - (aj + a + b) * comb(i, j)
// fはpascal_sum_row, fとgを同時に計算する
// g(i - 1, j) = (g(i, j) - a * f(i - 1, j) + (aj + a + b) * comb(i - 1, j)) / 2


template<typename mint>
std::vector<mint> pascal_sum_row_arithmetic_weighted(mint a, mint b, std::vector<std::pair<int, int>> Q){
  struct pascal_sum_row_arithmetic_struct{
    modcomb<mint> mcb;
    mint i2 = mint(2).inv();
    // [0, 1)が合うように初期値を設定
    // f(0, 0) = 1
    // g(0, 0) = b
    mint a, b, g, f = i2;

    pascal_sum_row_arithmetic_struct(mint _a, mint _b): a(_a), b(_b), g(b * i2 - a * i2 * i2){}

    void add_left(int l, int r){
      l++;
      g -= (a * l + b) * mcb.comb(r - 1, l);
      f -= mcb.comb(r - 1, l);
    }
    void del_left(int l, int r){
      l++;
      g += (a * l + b) * mcb.comb(r - 1, l);
      f += mcb.comb(r - 1, l);
    }
    void add_right(int l, int r){
      r--;
      g = 2 * g + a * f - (a * l + a + b) * mcb.comb(r, l);
      f = 2 * f - mcb.comb(r, l);
    }
    void del_right(int l, int r){
      r--;
      f = (f + mcb.comb(r, l)) * i2;
      g = (g - a * f + (a * l + a + b) * mcb.comb(r - 1, l)) * i2;
    }
  };
  int max_n = 0;
  for(auto [i, j] : Q) max_n = std::max(max_n, i); // nCrの最大のn
  pascal_sum_row_arithmetic_struct pas(a, b);
  pas.mcb.recalc(max_n + 1);
  mo_algorithm mo(max_n + 1, pas);
  int q = Q.size();
  std::vector<mint> ans(q);
  for(int i = 0; i < q; i++){
    if(Q[i].first < Q[i].second) Q[i].second = Q[i].first;
    mo.insert(Q[i].second, Q[i].first + 1);
  }
  for(int i = 0; i < q; i++){
    auto [idx, id] = mo.process2();
    ans[idx] = pas.g;
  }
  return ans;
}

// 等比級数の重み付き
// f(i, j) := ∑{0 <= k <= j} ab^k * comb(i, k)
// 遷移
// f(i, j) -> f(i, j + 1): ab^(j + 1) * comb(i, j + 1) を足す
// f(i, j) -> f(i, j - 1): ab^j * comb(i, j)を引く
// f(i, j) -> f(i + 1, j)
// f(i + 1, j) = f(i, j) + b * f(i, j - 1)
//             = (1 + b) * f(i, j) - b * ab^j * comb(i, j)
//
// f(i - 1, j) = (f(i, j) + b * ab^j * comb(i - 1, j)) / (1 + b)
template<typename mint>
std::vector<mint> pascal_sum_row_geometric_weighted(mint a, mint b, std::vector<std::pair<int, int>> Q){
  struct pascal_sum_row_geometric_struct{
    modcomb<mint> mcb;
    // [0, 1)が合うように初期値を設定
    // f(0, 0) = a
    const mint a, b, binv, b_plus_oneinv;
    mint a_b_pow_l, f;

    pascal_sum_row_geometric_struct(mint _a, mint _b): a(_a), b(_b), binv(b.inv()), b_plus_oneinv((b + 1).inv()), a_b_pow_l(a), f(a * b_plus_oneinv){}

    void add_left(int l, int r){
      l++;
      f -= a_b_pow_l * mcb.comb(r - 1, l);
      a_b_pow_l *= binv;
    }
    void del_left(int l, int r){
      l++;
      a_b_pow_l *= b;
      f += a_b_pow_l * mcb.comb(r - 1, l);
    }
    void add_right(int l, int r){
      r--;
      f = (1 + b) * f - b * a_b_pow_l * mcb.comb(r, l);
    }
    void del_right(int l, int r){
      r--;
      f = (f + b * a_b_pow_l * mcb.comb(r - 1, l)) * b_plus_oneinv;
    }
  };
  int max_n = 0;
  for(auto [i, j] : Q) max_n = std::max(max_n, i); // nCrの最大のn
  pascal_sum_row_geometric_struct pas(a, b);
  pas.mcb.recalc(max_n + 1);
  mo_algorithm mo(max_n + 1, pas);
  int q = Q.size();
  std::vector<mint> ans(q);
  for(int i = 0; i < q; i++){
    if(Q[i].first < Q[i].second) Q[i].second = Q[i].first;
    mo.insert(Q[i].second, Q[i].first + 1);
  }
  for(int i = 0; i < q; i++){
    auto [idx, id] = mo.process2();
    ans[idx] = pas.f;
  }
  return ans;
}


// f(i, j) := ∑{0 <= k <= i} comb(k, j)
template<typename mint>
mint pascal_sum_column(int i, int j){
  static modcomb<mint> mcb;
  mcb.recalc(i + 1);
  return mcb.comb(i + 1, j + 1);
}
// f(i, j) := ∑{0 <= k <= i} (ak + b) * comb(k, j)
// ガリガリ計算するとこの式になる
template<typename mint>
mint pascal_sum_column_arithmetic_weighted(mint a, mint b, int i, int j){
  static modcomb<mint> mcb;
  mcb.recalc(i + 1);
  return (a * i + b) * mcb.comb(i + 1, j + 1) - a * mcb.comb(i + 1, j + 2);
}

// ガリガリ計算すると
// f(i, j + 1) = (b * f(i, j) - a * b^(i + 1) * comb(i + 1, j + 1)) / (1 - b)
// f(i, j - 1) = ((1 - b) * f(i, j) + a * b ^ (i + 1) * comb(i + 1, j)) / b
// になる. moじゃなくてもできそうな気がする
/*
template<typename mint>
mint pascal_sum_column_geometric_weighted(mint a, mint b, int i, int j){

}
*/

// f(i, j) = ∑{0 <= k <= i}∑{0 <= l <= j} comb(k, l)
template<typename mint>
std::vector<mint> pascal_sum(std::vector<std::pair<int, int>> Q){
  int q = Q.size();
  for(int i = 0; i < q; i++){
    Q[i].first++;
    Q[i].second++;
  }
  auto res = pascal_sum_row<mint>(Q);
  for(int i = 0; i < q; i++) res[i] -= 1;
  return res;
}
#endif