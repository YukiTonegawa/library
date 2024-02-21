#include "../../template.hpp"
#include "../../math/matrix/matrix_mod2.hpp"

// 長さn(bit)のaの末尾に長さmのbを結合する
void concat_bitset64(int n, int m, std::vector<unsigned long long> &a, const std::vector<unsigned long long> &b){
  static constexpr int bitlen = 64;
  static constexpr int bitlen_shift = 6;
  static constexpr int bitlen_mod = 63;
  int N = (n + bitlen_mod) >> bitlen_shift;
  int M = (m + bitlen_mod) >> bitlen_shift;
  assert(a.size() == N && b.size() == M);
  if((n & bitlen_mod) == 0){
    a.insert(a.end(), b.begin(), b.end());
    return;
  }
  int K = (n + m + bitlen_mod) >> bitlen_shift; // 結合後のvectorのサイズ
  a.resize(K, 0);
  int s1 = n & bitlen_mod, s2 = bitlen - s1;
  for(int i = N - 1, j = 0; j <= M; i++, j++){
    if(j < M) a[i] |= b[j] << s1;
    if(j) a[i] |= b[j - 1] >> s2;
  }
}

// [k, n)をコピーして返す
std::vector<unsigned long long> copy_bitset64(int n, const std::vector<unsigned long long> &a, int k){
  using ull = unsigned long long;
  static constexpr int bitlen = 64;
  static constexpr int bitlen_shift = 6;
  static constexpr int bitlen_mod = 63;
  assert(0 <= k && k <= n);
  int N = (n + bitlen_mod) >> bitlen_shift;
  int M = (n - k + bitlen_mod) >> bitlen_shift;
  std::vector<ull> ret(M, 0);
  if((k & bitlen_mod) == 0){
    std::copy(a.begin() + (k >> bitlen_shift), a.end(), ret.begin());
    return ret;
  }
  int s1 = k & bitlen_mod, s2 = bitlen - s1;
  for(int i = (k >> bitlen_shift), j = 0; j < M; i++, j++){
    ret[j] = a[i] >> s1;
    if(i + 1 < N) ret[j] |= a[i + 1] << s2;
  }
  return ret;
}

/*
// [l, r)をコピーして返す
std::vector<unsigned long long> copy_range_bitset64(int n, const std::vector<unsigned long long> &a, int l, int r){
  using ull = unsigned long long;
  static constexpr int bitlen_shift = 6;
  static constexpr int bitlen_mod = 63;
  assert(k <= n);
  int N = (n + bitlen_mod) >> bitlen_shift;
  int M = (n - k + bitlen_mod) >> bitlen_shift;
  std::vector<ull> ret(M, 0);
  if((k & bitlen_mod) == 0){
    std::copy(a.begin() + (k >> bitlen_shift), a.end(), ret.begin());
    return ret;
  }
  int s1 = k & bitlen_mod;
  int s2 = n - s1;
  for(int i = (k >> bitlen_shift), j = 0; j < M; i++, j++){
    ret[j] = a[i] >> s1;
    if(i + 1 < N) ret[j] |= a[i + 1] << s2;
  }
  return ret;
}
*/

// [0, k), [k, n)
std::pair<std::vector<unsigned long long>, std::vector<unsigned long long>> split_bitset64(int n, std::vector<unsigned long long> &a, int k){
  static constexpr int bitlen_shift = 6;
  static constexpr int bitlen_mod = 63;
  assert(0 <= k && k <= n);
  auto b = copy_bitset64(n, a, k);
  a.resize((k + bitlen_mod) >> bitlen_shift);
  if((k & bitlen_mod) != 0){
    a.back() &= ((unsigned long long)1 << (k & bitlen_mod)) - 1;
  }
  return {a, b};
}

int main(){
  matrix_mod2 A(5, 4, 1);
  auto a = A.gaussian_elimination().first;
  a.print();
}


/*
template<typename T>
std::vector<bool> convolution_bitset_ntt(const std::vector<T> &a, const std::vector<T> &b){
  static constexpr int mod = 998244353;
  using mint = static_modint<mod>;
  int n = a.size(), m = b.size();
  if (!n || !m) return {1};
  std::vector<mint> a2(n), b2(m);
  for (int i = 0; i < n; i++) a2[i] = mint(a[i] & 1);
  for (int i = 0; i < m; i++) b2[i] = mint(b[i] & 1);
  auto c2 = _convolution<mint>(move(a2), move(b2));
  std::vector<bool> c(n + m - 1);
  for(int i = 0; i < n + m - 1; i++) c[i] = c2[i].val() > 0;
  return c;
}
template<typename T>
std::vector<bool> convolution_bitset_rk(std::vector<T> a, std::vector<T> b){
  int n = a.size(), m = b.size();
  if (!n || !m) return {1};
  std::vector<int> pos, pos2;
  for(int i = 0; i < n; i++) if(a[i] & 1) pos.push_back(i);
  for(int i = 0; i < m; i++) if(b[i] & 1) pos2.push_back(i);
  if(pos.size() > pos2.size()){
    std::swap(n, m);
    std::swap(a, b);
    std::swap(pos, pos2);
  }
  std::vector<bool> c(n + m - 1);
  int try_cnt = 0;
  for(int i = 0, j = 0; i < n + m - 1; i++){
    if(j < pos.size() && pos[j] == i) j++;
    bool f = false;
    for(int k = j - 1; !f && k >= 0 && i - pos[k] < m; k--, try_cnt++){
      f |= (b[i - pos[k]] & 1);
    }
    c[i] = f;
    if(try_cnt > 4 * (n + m - 1)) return convolution_bitset_ntt<T>(a, b);
  }
  return c;
}
*/