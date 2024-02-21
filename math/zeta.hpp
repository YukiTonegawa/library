#ifndef _ZETA_H_
#define _ZETA_H_

#include <vector>
#include <functional>
#include <limits>
#include <cassert>

namespace zeta{
  namespace subset{
    // z[S] = f(v[T], T⊂S)   O(N * 2^N)
    template<typename T>
    std::vector<T> zeta_subset(int n, std::vector<T> dp, const std::function<T(T, T)> &f){
      assert(dp.size() == (1 << n));
      int m = 1 << n;
      for(int i = 1; i < m; i <<= 1){
        for(int j = 0; j < m; j++){
          if(!(j & i)) dp[j | i] = f(dp[j | i], dp[j]);
        }
      }
      return dp;
    }
    // z[S] = sum(v[T], T⊂S)
    template<typename T>
    std::vector<T> zeta_subset_sum(int n, const std::vector<T> &v){
      return zeta_subset<T>(n, v, [](T l, T r){return l + r;});
    }
    // 包除原理(sumの逆変換)
    template<typename T>
    std::vector<T> mobius_subset_sum(int n, const std::vector<T> &v){
      return zeta_subset<T>(n, v, [](T l, T r){return l - r;});
    }
    // z[S] = max(v[T], T⊂S)
    template<typename T>
    std::vector<T> zeta_subset_max(int n, const std::vector<T> &v){
      return zeta_subset<T>(n, v, [](T l, T r){return max(l, r);});
    }
  };
  namespace superset{
    // z[S] = f(v[T], T⊃S)    O(N * 2^N)
    template<typename T>
    std::vector<T> zeta_superset(int n, std::vector<T> dp, const std::function<T(T, T)> &f){
      assert(dp.size() == (1 << n));
      int m = 1 << n;
      for(int i = 1; i < m; i <<= 1){
        for(int j = 0; j < m; j++){
          if(!(j & i)) dp[j] = f(dp[j], dp[j | i]);
        }
      }
      return dp;
    }
    // z[S] = sum(v[T], T⊃S)
    template<typename T>
    std::vector<T> zeta_superset_sum(int n, const std::vector<T> &v){
      return zeta_superset<T>(n, v, [](T l, T r){return l + r;});
    }
    // 包除原理(sumの逆変換)
    template<typename T>
    std::vector<T> mobius_superset_sum(int n, const std::vector<T> &v){
      return zeta_superset<T>(n, v, [](T l, T r){return l - r;});
    }
    // z[S] = max(v[T], T⊃S)
    template<typename T>
    std::vector<T> zeta_superset_max(int n, const std::vector<T> &v){
      return zeta_superset<T>(n, v, [](T l, T r){return max(l, r);});
    }
  };
  namespace intersect{
    // z[S] = f(v[T], (T&S)!=φ)    O(N * 2^N)
    // できる演算が限定されるかも？ sum, max, minはできそう
    // 最初の部分で
    // v[111] = 5 ->
    // A[111] = 5
    // A[011] = A[101] = A[110] = -5
    // A[100] = A[010] = A[001] = 5としたくなるが
    // 5 + inv(5) + inv(5) + inv(5) + 5 + 5 + 5 = 5 = v[111]が成り立つ
    // maxも逆元は定義できないが inv(任意の数) = 0 として
    // max({5, inv(5), inv(5), inv(5), 5, 5, 5}) = 5 = v[111]が成り立つためできるはず
    // つまり f(Tの立っているbitが奇数 ? A[T] : inv(A[T]), T⊂S) = v[S]ならできそう？

    template<typename T>
    std::vector<T> zeta_intersect(int n, std::vector<T> dp, const std::function<T(T, T)> &f, const std::function<T(T)> &inv){
      assert(dp.size() == (1 << n));
      int m = 1 << n;
      for(int i = 1; i < m; i <<= 1){
        for(int j = 0; j < m; j++){
          if(!(j & i)) dp[j] = f(dp[j], dp[j | i]);
        }
      }
      for(int j=0;j<m;j++){
        if(!(__builtin_popcount(j) & 1)) {
          dp[j] = inv(dp[j]);
        }
      }
      dp[0] = 0;
      for(int i = 1; i < m; i <<= 1){
        for(int j = 0; j < m; j++){
          if(!(j & i)) dp[j | i] = f(dp[j | i], dp[j]);
        }
      }
      return dp;
    }
    // z[S] = sum(v[T], (T&S)!=φ)
    template<typename T>
    std::vector<T> zeta_intersect_sum(int n, const std::vector<T> &v){
      return zeta_intersect<T>(n, v, [](T l, T r){return l + r;}, [](T l){return -l;});
    }
    // z[S] = max(v[T], (T&S)!=φ)
    template<typename T>
    std::vector<T> zeta_intersect_max(int n, const std::vector<T> &v){
      return zeta_intersect<T>(n, v, [](T l, T r){return std::max(l, r);},
                      [](T l){return std::numeric_limits<T>::min();});
    }
    // z[S] = min(v[T], (T&S)!=φ)
    template<typename T>
    std::vector<T> zeta_intersect_min(int n, const std::vector<T> &v){
      return zeta_intersect<T>(n, v, [](T l, T r){return std::min(l, r);},
                      [](T l){return std::numeric_limits<T>::max();});
    }
  };
  namespace disjoint{
    // z[S] = f(v[T], (T&S)==φ)    O(N * 2^N)
    // Z[S] = f(v[T], T⊂(~S))
    template<typename T>
    std::vector<T> zeta_disjoint(int n, const std::vector<T> &v, const std::function<T(T, T)> &f){
      int m = 1 << n;
      auto z = subset::zeta_subset<T>(n, v, f);
      std::vector<T> ret(m);
      for(int i = 0; i < m;i++) ret[i] = z[(m - 1) ^ i];
      return ret;
    }
    template<typename T>
    std::vector<T> mobius_disjoint(int n, const std::vector<T> &v, const std::function<T(T, T)> &f){
      int m = 1 << n;
      std::vector<T> ret(m);
      for(int i = 0; i < m; i++) ret[i] = v[(m - 1) ^ i];
      return subset::zeta_subset<T>(n, ret, f);
    }
    // z[S] = sum(v[T], (T&S)==φ)
    template<typename T>
    std::vector<T> zeta_disjoint_sum(int n, const std::vector<T> &v){
      return zeta_disjoint<T>(n, v, [](T l, T r){return l + r;});
    }
    // sumの逆変換
    template<typename T>
    std::vector<T> mobius_disjoint_sum(int n, const std::vector<T> &v){
      return mobius_disjoint<T>(n, v, [](T l, T r){return l - r;});
    }
    // z[S] = max(v[T], (T&S)==φ)
    template<typename T>
    std::vector<T> zeta_disjoint_max(int n, const std::vector<T> &v){
      return zeta_disjoint<T>(n, v, [](T l, T r){return std::max(l, r);});
    }
  };
  namespace multiple{
    //1-indexed, z[i] = f(v[j], 1<=i<=j<=n and j%i==0) O(NlogN)
    template<typename T>
    std::vector<T> zeta_multiple(int n, const std::vector<T> &v, const std::function<T(T, T)> &f, T e){
      std::vector<T> ret(n, e);
      for(int i = 0; i < v.size(); i++) ret[i] = v[i];
      for(int i = 1; i <= n; i++){
        for(int j = 2 * i; j <= n; j += i){
          ret[i - 1] = f(ret[i - 1], ret[j - 1]);
        }
      }
      return ret;
    }
    template<typename T>
    std::vector<T> mobius_multiple(int n, const std::vector<T> &v, const std::function<T(T, T)> &f, T e){
      std::vector<T> ret(n, e);
      for(int i = 0; i < v.size(); i++) ret[i] = v[i];
      for(int i = n; i >= 1; i--){
        for(int j = 2 * i; j <= n; j += i){
          ret[i - 1] = f(ret[i - 1], ret[j - 1]);
        }
      }
      return ret;
    }
  };
  namespace divisor{
    //1-indexed, z[i] = f(v[j], 1<=j<=i<=n and i%j==0) O(NlogN)
    template<typename T>
    std::vector<T> zeta_divisor(int n, const std::vector<T> &v, const std::function<T(T, T)> &f, T e){
      std::vector<T> ret(n, e);
      for(int i = 0; i < v.size(); i++) ret[i] = v[i];
      for(int i = n; i >= 1; i--){
        for(int j = 2 * i; j <= n; j += i){
          ret[j - 1] = f(ret[j - 1], ret[i - 1]);
        }
      }
      return ret;
    }
    template<typename T>
    std::vector<T> mobius_divisor(int n, const std::vector<T> &v, const std::function<T(T, T)> &f, T e){
      std::vector<T> ret(n, e);
      for(int i = 0; i < v.size(); i++) ret[i] = v[i];
      for(int i = 1; i <= n; i++){
        for(int j = 2 * i; j <= n; j += i){
          ret[j - 1] = f(ret[j - 1], ret[i - 1]);
        }
      }
      return ret;
    }
  };
};

namespace convolution{
  template<typename T, int max_log2 = 20>
  std::vector<T> convolution_subset(int n, const std::vector<T> &A, const std::vector<T> &B){
    assert(n <= max_log2);
    assert(A.size() == (1 << n));
    int m = 1 << n;
    std::vector<std::array<std::pair<T, T>, max_log2 + 1>> C(m);
    for(int i = 0; i < m; i++){
      C[i].fill({0, 0});
      C[i][__builtin_popcount(i)] = {A[i], B[i]};
    }
    for(int i = 1; i < m; i <<= 1){
      for(int j = 0; j < m; j++){
        int max_size = __builtin_popcount(j | i);
        if(!(j & i)){
          for(int k = 0; k <= max_size; k++){
            C[j | i][k].first += C[j][k].first;
            C[j | i][k].second += C[j][k].second;
          }
        }
      }
    }
    for(int j = 0; j < m; j++){
      for(int k = 1; k <= n; k++){
        C[j][k].first += C[j][k - 1].first;
        C[j][k].second += C[j][k - 1].second;
      }
    }
    for(int i = 0; i < m; i++){
      for(int j = n; j >= 0; j--){
        T cur = 0;
        for(int k = 0; k <= j; k++) cur += C[i][j - k].first * C[i][k].second;
        C[i][j].first = cur;
      }
    }
    for(int i = 1; i < m; i <<= 1){
      for(int j = 0; j < m; j++){
        if(!(j & i)){
          for(int k = 0; k <= n; k++) C[j | i][k].first -= C[j][k].first;
        }
      }
    }
    std::vector<T> res(m);
    for(int i = 0; i < m; i++){
      int cnt = __builtin_popcount(i);
      res[i] = cnt ? C[i][cnt].first - C[i][cnt - 1].first : C[i][0].first;
    }
    return res;
  }
  // 2つの配列の畳み込み
  // ・or, and, xor(アダマール変換), subset
  // ・gcd, lcm
  // ・+(FFT, NTT), *(原始根を使ったNTT)
  using namespace zeta;
  // z[k] = ∑{i, j, (i|j)==k} A[i]*B[j] O(N * 2^N)
  template<typename T>
  std::vector<T> convolution_or(int n, const std::vector<T> &A, const std::vector<T> &B){
    std::vector<T> Az = subset::zeta_subset_sum(n, A);
    std::vector<T> Bz = subset::zeta_subset_sum(n, B);
    for(int i = 0; i < (1 << n); i++) Az[i] *= Bz[i];
    //return Az;  z[k] = ∑{i, j, (i|j)⊂k} A[i]*B[j] が欲しい場合
    return subset::mobius_subset_sum(n, Az);
  }
  // z[k] = ∑{i, j, (i&j)==k} A[i]*B[j] O(N * 2^N)
  template<typename T>
  std::vector<T> convolution_and(int n, const std::vector<T> &A, const std::vector<T> &B){
    std::vector<T> Az = superset::zeta_superset_sum(n, A);
    std::vector<T> Bz = superset::zeta_superset_sum(n, B);
    for(int i=0;i < (1 << n); i++) Az[i] *= Bz[i];
    //return Az;  z[k] = ∑{i, j, (i&j)⊂k} A[i]*B[j] が欲しい場合
    return superset::mobius_superset_sum(n, Az);
  }
  // z[k] = ∑{i, j, (i^j)==k} A[i]*B[j] O(N * 2^N)
  template<typename T>
  std::vector<T> convolution_xor(int n, std::vector<T> A, std::vector<T> B){
    assert(A.size() == (1 << n) && B.size() == (1 << n));
    int m = 1 << n;
    for(int i = 1; i < m; i <<= 1){
      for(int j = 0; j < m; j++){
        if(!(j&i)){
          T l = A[j], r = A[j|i];
          A[j] = l + r, A[j|i] = l - r;
          T L = B[j], R = B[j|i];
          B[j] = L + R, B[j|i] = L - R;
        }
      }
    }
    //T i2 = (998244353 + 1) / 2;
    for(int i = 0; i < m; i++) A[i] *= B[i];
    for(int i = 1; i < m; i <<= 1){
      for(int j = 0; j < m; j++){
        if(!(j & i)){
          T l = A[j], r = A[j|i];
          A[j] = (l + r) / 2, A[j | i] = (l - r) / 2;
          //A[j] = (l + r) * i2, A[j | i] = (l - r) * i2; // 素数modの場合
        }
      }
    }
    return A;
  }
  // z[k] = ∑{i, j, min(i, j) == k} A[i]*B[j] O(N)
  template<typename T>
  std::vector<T> convolution_min(const std::vector<T> &A, const std::vector<T> &B){
    int N = std::min(A.size(), B.size());
    std::vector<T> Az(N, 0), Bz(N, 0);
    for(int i = N - 1; i >= 0; i--){
      if(i != N - 1) Az[i] += Az[i + 1];
      Az[i] += A[i];
    }
    for(int i = N - 1; i >= 0; i--){
      if(i != N - 1) Bz[i] += Bz[i + 1];
      Bz[i] += B[i];
    }
    std::vector<T> Cz(N);
    for(int i = 0; i < N; i++) Cz[i] = Az[i] * Bz[i];
    for(int i = 0; i < N - 1; i++) Cz[i] -= Cz[i + 1];
    return Cz;
  }
  // z[k] = ∑{i, j, max(i, j) == k} A[i]*B[j] O(N)
  template<typename T>
  std::vector<T> convolution_max(const std::vector<T> &A, const std::vector<T> &B){
    int N = max(A.size(), B.size());
    std::vector<T> Az(N, 0), Bz(N, 0);
    for(int i = 0; i < N; i++){
      if(i) Az[i] += Az[i - 1];
      if(i < A.size()) Az[i] += A[i];
    }
    for(int i = 0; i < N; i++){
      if(i) Bz[i] += Bz[i - 1];
      if(i < B.size()) Bz[i] += B[i];
    }
    std::vector<T> Cz(N);
    for(int i = 0; i < N; i++) Cz[i] = Az[i] * Bz[i];
    for(int i = N - 1; i >= 1; i--) Cz[i] -= Cz[i - 1];
    return Cz;
  }
  // 1-indexed, z[k] = ∑{i, j, gcd(i, j)==k} A[i]*B[j]
  // gcd(i, j)は最大でmin(|A|, |B|)になる    O(max(NlogN, MlogM))
  template<typename T>
  std::vector<T> convolution_gcd(const std::vector<T> &A, const std::vector<T> &B){
    int n = A.size();
    int m = B.size();
    std::vector<T> Az = multiple::zeta_multiple<T>(n, A, [](T l, T r){return l + r;}, 0);
    std::vector<T> Bz = multiple::zeta_multiple<T>(m, B, [](T l, T r){return l + r;}, 0);
    std::vector<T> Cz(std::max(n, m), 0);
    for(int i = 0; i < std::min(n, m); i++) Cz[i] = Az[i] * Bz[i];
    return multiple::mobius_multiple<T>(std::max(n, m), Cz, [](T l, T r){return l - r;}, 0);
  }
  // 1-indexed, z[k] = ∑{i, j, lcm(i, j)==k} A[i]*B[j]
  // lcm(i, j)は最大で|A|*|B|になる     O(NMlog(NM))
  // https://atcoder.jp/contests/agc038/submissions/18538058　
  // のようにgcd畳み込みを変形させるパターンの方が多いかも
  template<typename T>
  std::vector<T> convolution_lcm(const std::vector<T> &A, const std::vector<T> &B, int limit_dim){
    int N = std::min((long long)A.size() * (long long)B.size(), (long long)limit_dim);
    std::vector<T> Az = divisor::zeta_divisor<T>(N, A, [](T l, T r){return l + r;}, 0);
    std::vector<T> Bz = divisor::zeta_divisor<T>(N, B, [](T l, T r){return l + r;}, 0);
    std::vector<T> Cz(N);
    for(int i = 0; i < N; i++) Cz[i] = Az[i] * Bz[i];
    return divisor::mobius_divisor<T>(N, Cz, [](T l, T r){return l - r;}, 0);
  }
};
#endif