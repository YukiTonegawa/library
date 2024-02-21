#ifndef _SUBSET_SUM_H_
#define _SUBSET_SUM_H_

#include <vector>
#include <cmath>
#include <algorithm>

// 1, 2...nからそれぞれ最大1個使ってSを作る, 作れない場合は空配列, O(√S)
std::vector<long long> subset_sum_linear(long long n, long long S){
  if(S <= 0 || n * (n + 1) / 2 < S) return {};
  std::vector<long long> ret;
  for(int i = n; S && i >= 1; i--){
    if(i <= S){
      ret.push_back(i);
      S -= i;
    }
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}
// 1, 2...nからそれぞれ最大1個, 合計k個使ってSを作る, 作れない場合は空配列, O(√S)
std::vector<long long> subset_sum_linear_k(long long n, long long S, int k){
  auto calc_max = [](long long N, long long K){
    return N * (N + 1) / 2 - (N - K) * (N - K + 1) / 2;
  };
  long long min_s = (long long)k * (k + 1) / 2;
  if(S < min_s || calc_max(n, k) < S) return {};
  std::vector<long long> ret;
  for(int i = n; S && i >= 1; i--){
    if(calc_max(i - 1, k) < S){
      ret.push_back(i);
      S -= i;
      k--;
    }
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}
// l, l + 1...r - 1だけを使ってSを作る, 作れない場合は空配列　, O(√S)
std::vector<long long> subset_sum_range(long long l, long long r, long long S){
  // k個使うと決め打つと 1...r-1-lからk個使ってS-klを作る問題になる
  // k個使うときの最小値　
  auto calc_min = [l](long long k){
    return k * l + k * (k + 2) / 2;
  };
  // k個使うときの最大値
  auto calc_max = [l, r](long long k){
    assert(r - l >= k);
    return r * (r - 1) / 2 - (r - k) * (r - k - 1) / 2;
  };
  long long L = 0, R = r - l + 1;
  while(R - L > 1){
    long long M = (L + R) / 2;
    if(calc_min(M) > S) R = M;
    else L = M;
  }
  long long Rmin = R; // これ以上多くの要素を使うと最小値がSを超える
  L = -1, R = r - l + 1;
  while(R - L > 1){
    long long M = (L + R) / 2;
    if(calc_max(M) >= S) R = M;
    else L = M;
  }
  long long Lmax = R; // これより小さいと最大値がS未満になる
  if(Lmax > r - l || Rmin <= Lmax) return {};
  // Lmax個使う
  auto ret = subset_sum_linear_k(r - l, S - Lmax * (l - 1), Lmax);
  for(auto &v : ret) v += l - 1;
  return ret;
}
// 2C2, 3C2, 4C2...をそれぞれ最大1個足してSを作る, 作れない場合は空配列, O(S^(1/3)))
// S > 35の時必ず構築できる
std::vector<long long> subset_sum_nc2(long long S){
  static std::vector<std::vector<long long>> low_table;
  static constexpr int threshold = 35;
  if(low_table.empty()){
    std::vector<long long> sq{0, 1, 3, 6, 10, 15, 21, 28};
    low_table.resize(threshold + 1);
    for(int i = 1; i <= threshold; i++){
      for(int j = 1; j < sq.size(); j++){
        int k = i - sq[j];
        if(k < 0) break;
        if(k == 0 || (!low_table[k].empty() && low_table[k].back() < j + 1)){
          low_table[i] = low_table[k];
          low_table[i].push_back(j + 1);
          break;
        }
      }
    }
  }
  std::vector<long long> ret;
  // xC2が_S以下になる最大のx
  auto find_x = [](long long _S){
    long long l = 0, r = _S;
    while(r - l > 1){
      long long m = (l + r) / 2;
      if(m * (m - 1) > 2 * _S) r = m;
      else l = m;
    }
    return l;
  };
  while(S > threshold){
    long long x = find_x(S);
    for(; ; x--){
      long long T = S - (x * (x - 1) / 2);
      if(T > threshold || T == 0 || (!low_table[T].empty() && low_table[T].back() < x)){
        ret.push_back(x);
        S = T;
        break;
      }
    }
  }
  std::reverse(ret.begin(), ret.end());
  if(S) ret.insert(ret.begin(), low_table[S].begin(), low_table[S].end());
  return ret;
}

// 1^2, 2^2...をそれぞれ最大1個足してSを作る, 作れない場合は空配列, O(S^(1/3)))
// S > 128の時必ず構築できる
std::vector<long long> subset_sum_square(long long S){
  static std::vector<std::vector<long long>> low_table;
  static constexpr int threshold = 128;
  if(low_table.empty()){
    std::vector<long long> sq{1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121};
    low_table.resize(threshold + 1);
    for(int i = 1; i <= threshold; i++){
      for(int j = 0; j < sq.size(); j++){
        int k = i - sq[j];
        if(k < 0) break;
        if(k == 0 || (!low_table[k].empty() && low_table[k].back() < j + 1)){
          low_table[i] = low_table[k];
          low_table[i].push_back(j + 1);
          break;
        }
      }
    }
  }
  std::vector<long long> ret;
  while(S > threshold){
    long long x = sqrtl(S);
    for(; ; x--){
      long long T = S - x * x;
      if(T > threshold || T == 0 || (!low_table[T].empty() && low_table[T].back() < x)){
        ret.push_back(x);
        S = T;
        break;
      }
    }
  }
  std::reverse(ret.begin(), ret.end());
  if(S) ret.insert(ret.begin(), low_table[S].begin(), low_table[S].end());
  return ret;
}
#endif