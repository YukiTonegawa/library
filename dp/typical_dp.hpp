#ifndef _TYPICAL_DP_H_
#define _TYPICAL_DP_H_
#include <string>
#include <cassert>

// 桁dpの種類

// [1] 文字列を整数とみなして[L, R)かつ条件を満たすようなものを数え上げる
// [1]では, leading-zeroと他の位置での0を区別する必要がある場合がある

// n桁に揃える. ascii codeで'0'の前の'/'をleading-zeroとする
std::string padding_slash(const std::string &s, int n){
  int m = s.size();
  assert(m <= n);
  return std::string(n - m, '/') + s;
}
void keta_dp_range_template(std::string l, const std::string &r){
  assert(l.size() <= r.size());
  int n = r.size();
  l = padding_slash(l, n);
  for(int i = 0; i < n; i++){
    char lc = l[i], rc = r[i];
    // 0 := 現時点でl, rと一致している
    for(char c = lc; c <= rc; c++){
      if(lc == '/' && c == '0') continue; // leading-zero
      // 遷移先
      // int flag = (int)(lc < c) + 2 * (int)(c < rc);
      // dp[i + 1][flag] += dp[i][0];
    }
    // 1 := 現時点でLより大きく, Rと一致している
    for(char c = '0'; c <= rc; c++){
      // 遷移先
      // int flag = 1 + 2 * (int)(c < rc);
      // dp[i + 1][flag] += dp[i][1];
    }
    // 2 := 現時点でLと一致, Rより小さい
    for(char c = lc; c <= '9'; c++){
      if(lc == '/' && c == '0') continue;// leading-zero
      // 遷移先
      // int flag = (int)(lc < c) + 2;
      // dp[i + 1][flag] += dp[i][2];
    }
    // 3 := 現時点でLより大きく, Rより小さい
    for(char c = '0'; c <= '9'; c++){
      // 遷移先
      // int flag = 3;
      //dp[i + 1][flag] += dp[i][3];
    }
  }
}

//1 + ax^k (mod P) を掛ける O(N)
template<typename mint>
void multiply_binom_inplace(mint a, int k, std::vector<mint> &v){
  assert(k >= 0);
  int n = v.size();
  for(int i = n - 1; i >= k; i--) v[i] += a * v[i - k];
}
template<typename mint>
std::vector<mint> multiply_binom(mint a, int k, const std::vector<mint> &v){
  std::vector<mint> res = v;
  multiply_binom_inplace<mint>(a, k, res);
  return res;
}
//1 + ax^k (mod P)で割る O(N)
//1 - ax^k + a^2x^2k - a^3x^3k...を掛ける事と同じ
template<typename mint>
void divide_binom_inplace(mint a, int k, std::vector<mint> &v){
  a *= -1;
  assert(k);
  int n = v.size();
  std::vector<std::vector<mint>> sum(k, std::vector<mint>(n / k + 1, 0)); //kで割ったあまり, 累積和
  for(int i = 0; i < n; i++){
    int r = i % k, idx = i / k;
    if(idx == 0) sum[r][idx] = v[i];
    else sum[r][idx] = v[i] + a * sum[r][idx - 1];
    v[i] = sum[r][idx];
  }
}
template<typename mint>
std::vector<mint> divide_binom(mint a, int k, const std::vector<mint> &v){
  std::vector<mint> res = v;
  divide_binom_inplace<mint>(a, k, res);
  return res;
}
#endif