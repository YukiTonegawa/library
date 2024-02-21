#ifndef _MATRIX_MOD_H_
#define _MATRIX_MOD_H_
#include "../mod.hpp"
#include <vector>

template<typename mint>
struct matrix_mod{
  int n, m;
  using _mint = mint;
private:
  using vec = std::vector<mint>;
  using matrix = matrix_mod<mint>;
  // n × k 行列と k × m 行列の積(n × m行列)
  // K == 0だと壊れる
  static matrix __mul_mat(const matrix &vl, const matrix &vr){
    int N = vl.n, K = vl.m, M = vr.m;
    assert(K == vr.n);
    assert(K);
    if(N == 0) return matrix(0, M);
    if(M == 0) return matrix(N, 0);
    auto vr_t = vr.t();
    matrix ret(N, M, 0);
    for(int i = 0; i < N; i++){
      for(int j = 0; j < M; j++){
        __int128_t S = 0;
        for(int k = 0; k < K; k++){
          S += (long long)vl.val[i][k].val() * vr_t[j][k].val();
        }
        ret[i][j] = S % mint::mod();
      }
    }
    return ret;
  }
  // n × m 行列と n × m 行列の和(n × m行列)
  static void __add_mat_inplace(matrix &vl, const matrix &vr){
    assert(vl.n == vr.n && vl.m == vr.m);
    int N = vl.n, M = vl.m;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < M; j++){
        vl[i][j] += vr[i][j];
      }
    }
  }
  // n × m 行列と n × m 行列の差(n × m行列)
  static void __sub_mat_inplace(matrix &vl, const matrix &vr){
    assert(vl.n == vr.n && vl.m == vr.m);
    int N = vl.n, M = vl.m;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < M; j++){
        vl[i][j] -= vr[i][j];
      }
    }
  }
  static void __mul_val_inplace(matrix &vl, mint vr){
    int N = vl.n, M = vl.m;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < M; j++){
        vl[i][j] *= vr;
      }
    }
  }
  static void __add_val_inplace(matrix &vl, mint vr){
    int N = vl.n, M = vl.m;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < M; j++){
        vl[i][j] += vr;
      }
    }
  }
  static void __sub_val_inplace(matrix &vl, mint vr){
    int N = vl.n, M = vl.m;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < M; j++){
        vl[i][j] -= vr;
      }
    }
  }
  std::vector<vec> val;
public:
  matrix_mod(): n(0), m(0){}
  matrix_mod(int _n, int _m, mint x = mint(0)) : n(_n), m(_m), val(_n, vec(_m, x)){}
  matrix_mod(const matrix_mod &v) : n(v.n), m(v.m), val(v.val){}
  matrix_mod(const vec &v): n(1), m(v.size()), val(1, vec(v.size())){val[0] = v;}
  matrix_mod(const std::vector<vec> &v): n(v.size()), m(v[0].size()), val(v){}

  matrix_mod operator +  (const matrix_mod &vr){matrix_mod tmp(*this); return tmp += vr;}
  matrix_mod operator -  (const matrix_mod &vr){matrix_mod tmp(*this); return tmp -= vr;}
  matrix_mod operator *  (const matrix_mod &vr){return __mul_mat(*this, vr);}
  matrix_mod operator ^  (const long long vr){return pow(vr);}
  matrix_mod operator *  (const mint vr){matrix_mod tmp(*this); return tmp *= vr;}
  matrix_mod operator += (const matrix_mod &vr){__add_mat_inplace(*this, vr); return *this;}
  matrix_mod operator -= (const matrix_mod &vr){__sub_mat_inplace(*this, vr); return *this;}
  matrix_mod operator *= (const matrix_mod &vr){return (*this) = __mul_mat(*this, vr);}
  matrix_mod operator ^= (const long long vr){return (*this) = pow(vr);}
  matrix_mod operator *= (const mint vr){__mul_val_inplace(*this, vr); return *this;}
  vec& operator [] (const int i){return val[i];}

  // n次の単位行列
  static matrix_mod eye(int n){
    matrix_mod ret(n, n, 0);
    for(int i = 0; i < n; i++) ret[i][i] = mint(1);
    return ret;
  }
  void print()const{
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        std::cout << val[i][j] << (j == m - 1 ? '\n' : ' ');
      }
    }
  }
  matrix_mod pow(long long k){
    assert(n && m && n == m); // 正方行列でなければならない
    matrix_mod ret = eye(n); // k == 0の場合単位行列を返す
    matrix_mod m(*this);
    while(k){
      if(k & 1) ret *= m;
      m *= m;
      k >>= 1;
    }
    return ret;
  }
  // 転置
  matrix_mod t()const{
    matrix_mod ret(m, n, 0);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        ret[j][i] = val[i][j];
      }
    }
    return ret;
  }
  //掃き出し法で上三角行列を作る, {変形後の行列、行のスワップ回数}を返す O(NM^2)
  std::pair<matrix_mod, int> gaussian_elimination(){
    matrix_mod v(*this);
    int row = 0;//確定していない行
    int swp = 0;
    for(int i = 0; i < m && row < n; i++){
      //i列目が0でない行を探す
      int r = -1;
      for(int j = row; j < n; j++){
        if(v[j][i].val()){
          r = j;
          break;
        }
      }
      if(r == -1) continue;
      if(r != row){
        swp++;
        std::swap(v[r], v[row]);
      }
      //i列目が0でない行の処理
      for(int j = row + 1; j < n; j++){
        if(v[j][i].val() == 0) continue;
        mint x = v[j][i] / v[row][i];
        for(int k = i; k < m; k++){
          v[j][k] -= x * v[row][k];
        }
      }
      row++;
    }
    return {v, swp};
  }
  //掃き出し法で上三角行列を作る, {変形後の行列、行のスワップ回数}を返す O(NM^2 * log mod)
  std::pair<matrix_mod, int> gaussian_elimination_arbitrary_mod(){
    matrix_mod v(*this);
    int row = 0;//確定していない行
    int swp = 0;
    for(int i = 0; i < m && row < n; i++){
      //i列目が0でない行を探す
      int r = -1;
      for(int j = row; j < n; j++){
        if(v[j][i].val()){
          r = j;
          break;
        }
      }
      if(r == -1) continue;
      if(r != row){
        swp++;
        std::swap(v[r], v[row]);
      }
      //i列目が0でない行の処理
      for(int j = row + 1; j < n; j++){
        while(v[j][i].val() != 0){
          if(v[row][i].val() > v[j][i].val()){
            swp++;
            std::swap(v[row], v[j]);
          }
          int x = v[j][i].val() / v[row][i].val();
          for(int k = i; k < m; k++){
            v[j][k] -= x * v[row][k];
          }
        }
      }
      row++;
    }
    return {v, swp};
  }
  //すでに上三角行列になっていることが前提
  int rank(){
    int cnt = 0;
    for(int i = 0; i < n; i++, cnt++){
      bool f = false;
      for(int j = i; j < m; j++){
        if(val[i][j].val()){
          f = true;
          break;
        }
      }
      if(!f) break;
    }
    return cnt;
  }
  // 行列式 O(N^3)
  mint det(){
    assert(n == m); // 正方行列のみ
    auto [tmp, swp] = gaussian_elimination();
    mint res(1);
    for(int i = 0; i < n; i++) res *= tmp[i][i];
    return swp & 1 ? -res : res;
  }
  // 行列式 O(N^3 * log mod)
  mint det_arbitrary_mod(){
    assert(n == m); // 正方行列のみ
    auto [tmp, swp] = gaussian_elimination_arbitrary_mod();
    mint res(1);
    for(int i = 0; i < n; i++) res *= tmp[i][i];
    return swp & 1 ? -res : res;
  }
  // (n, m) + (n, l) -> (n, m + l)　横に結合
  matrix_mod concat_horizontal(matrix_mod vr){
    assert(n == vr.n);
    matrix_mod res(*this);
    for(int i = 0; i < n; i++){
      res[i].insert(res[i].end(), vr[i].begin(), vr[i].end());
    }
    res.m += vr.m;
    return res;
  }
  // (n, m) + (l, m) -> (n + l, m) 縦に結合
  matrix_mod concat_vertical(matrix_mod vr){
    assert(m == vr.m);
    matrix_mod res(*this);
    for(int i = 0; i < vr.n; i++) res.val.push_back(vr[i]);
    res.n += vr.n;
    return res;
  }
  // (n, m) -> (n, k), (n, m - k)
  std::pair<matrix_mod, matrix_mod> split_horizontal(int k){
    assert(0 <= k && k <= m);
    matrix_mod a(n, k), b(n, m - k);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < k; j++){
        a[i][j] = val[i][j];
      }
    }
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m - k; j++){
        b[i][j] = val[i][j + k];
      }
    }
    return {a, b};
  }
  // (n, m) -> (k, m), (n - k, m)
  std::pair<matrix_mod, matrix_mod> split_vertical(int k){
    assert(0 <= k && k <= n);
    matrix_mod a(k, m), b(n - k, m);
    for(int i = 0; i < k; i++){
      for(int j = 0; j < m; j++){
        a[i][j] = val[i][j];
      }
    }
    for(int i = 0; i < n - k; i++){
      for(int j = 0; j < m; j++){
        b[i][j] = val[k + i][j];
      }
    }
    return {a, b};
  }
  matrix_mod inv(){
    assert(n == m);
    auto [tmp, swp] = concat_horizontal(eye(n)).gaussian_elimination();
    for(int i = 0; i < n; i++){
      mint x = tmp[i][i];
      if(!x.val()) return matrix_mod{};// 存在しない
      x = x.inv();
      for(int j = 0; j < 2 * n; j++) tmp[i][j] *= x;
    }
    for(int i = n - 1; i >= 0; i--){
      for(int j = i + 1; j < n; j++){
        if(!tmp[i][j].val()) continue;
        mint c = tmp[i][j];
        for(int k = j; k < 2 * n; k++){
          tmp[i][k] -= c * tmp[j][k];
        }
      }
    }
    return tmp.split_horizontal(n).second;
  }
  // https://ja.wikipedia.org/wiki/LU%E5%88%86%E8%A7%A3
  // キャッシュのためにuを転置して実装
  std::pair<matrix_mod, matrix_mod> lu_decomposition(){
    matrix_mod l = eye(n), u(n, n);
    for(int i = 0; i < n; i++){
      // u[i][i]を決定
      u[i][i] = val[i][i];
      for(int j = 0; j < i; j++) u[i][i] -= l[i][j] * u[i][j];
      if(u[i][i].val() == 0) return {matrix_mod{}, matrix_mod{}}; // 不可能
      mint iuii = u[i][i].inv();
      // l[0, n)[i]を決定
      for(int j = i + 1; j < n; j++){
        l[j][i] = val[j][i];
        for(int k = 0; k < i; k++) l[j][i] -= l[j][k] * u[i][k];
        l[j][i] *= iuii;
      }
      // u[i][0, n)を決定
      for(int j = i + 1; j < n; j++){
        u[j][i] = val[i][j];
        for(int k = 0; k < i; k++) u[j][i] -= l[i][k] * u[j][k];
      }
    }
    u = u.t();
    return {l, u};
  }
  // Ax = b
  // (n, m) * (m, 1) -> (n, 1)
  // を満たす連立方程式を解く, 解空間の次元、(rank*変数)の基底を返す
  // 解空間の基底は任意のt_iについてA * (v1t1 + v2t2 ...) = 0を満たす
  // つまり plus + res[0]t_0 + res[1]t_1 + res[2]t_2...は全て解を満たす
  // 解が存在しない場合解空間の次元として-1を返す
  std::tuple<int, matrix_mod, vec> system_of_linear_equations(const vec &vr){
    assert(vr.size() == n);
    matrix_mod tmp = concat_horizontal(matrix_mod(vr).t()).gaussian_elimination().first;
    //解空間の次元 = 変数の数 - 階数
    int r = tmp.rank();
    std::vector<int> fc(r, -1);//各行に初めて非零要素が現れる列
    for(int i = 0; i < r; i++){
      mint tmp_inv;
      bool f = false;
      for(int j = i; j < tmp.m; j++){
        if(tmp[i][j].val() == 0) continue;
        if(j == tmp.m - 1 && !f){
          return {-1, matrix_mod{}, vec{}}; // 解なし
        }
        if(!f){
          tmp_inv = tmp[i][j].inv();
          fc[i] = j;
          f = true;
        }
        tmp[i][j] = tmp[i][j] * tmp_inv;
      }
    }
    int d = tmp.m - 1 - r, v = tmp.m - 1;
    vec plus(v, 0);
    for(int i = r - 1; i >= 0; i--){
      int idx = fc[i];
      assert(idx != -1);
      plus[idx] = tmp[i][v];
      for(int j = idx + 1; j < v; j++){
        plus[idx] -= plus[j] * tmp[i][j];
      }
    }
    matrix_mod res(d, v, 0);
    std::vector<bool> not_fc(v, true);
    for(int i = 0; i < r; i++) not_fc[fc[i]] = false;
    for(int i = 0, j = 0; i < v; i++) if(not_fc[i]) res[j++][i] = 1;
    for(int i = r - 1; i >= 0; i--){ //各行に1つまだ確定していない変数が現れる
      int col = fc[i];
      assert(col != -1);
      assert(tmp[i][col].val() == 1);
      for(int k = 0; k < d; k++){ // 次元
        for(int j = col + 1; j < v; j++){ // すでに確定した要素
          res[k][col] -= res[k][j] * tmp[i][j];
        }
      }
    }
    return {d, res, plus};
  }
};
#endif