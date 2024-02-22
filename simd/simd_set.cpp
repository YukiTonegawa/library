/*
#pragma GCC optimize ("O3")
#pragma GCC target ("avx2, avx512")

#include <immintrin.h>
#include ".lib/template.hpp"
#include ".lib/fast_io.hpp"
#include ".lib/data_structure/range_query/mo_algorithm.hpp"

// [0, 262144)
struct simd_set{
  static constexpr int D = 6;
  static constexpr int split = 8;
  static constexpr int split_mod = 7; // k % split = k & split_mod
  static constexpr int split_div = 3; // k / split = k >> split_div
  static constexpr int inf = 1 << (split_div * D); // 8^D
  static constexpr int node_count = (inf - 1) / split_mod; // 1 + 8 + 8^2...8^(D-1)
  using i256 = __m256i;
  std::array<i256, split> one; // i番目のみ1 (0 <= i < split)
  std::array<i256, split> mask; // [0, i)番目のみ111...1(0 <= i < split)
  std::array<i256, node_count> val;
  simd_set(){
    one.fill(_mm256_set1_epi64x(0));
    one[0] = _mm256_insert_epi32(one[0], 1, 0);
    one[1] = _mm256_insert_epi32(one[1], 1, 1);
    one[2] = _mm256_insert_epi32(one[2], 1, 2);
    one[3] = _mm256_insert_epi32(one[3], 1, 3);
    one[4] = _mm256_insert_epi32(one[4], 1, 4);
    one[5] = _mm256_insert_epi32(one[5], 1, 5);
    one[6] = _mm256_insert_epi32(one[6], 1, 6);
    one[7] = _mm256_insert_epi32(one[7], 1, 7);
    mask[0] = _mm256_set1_epi64x(0);
    mask[1] = _mm256_insert_epi32(mask[0], ~0, 0);
    mask[2] = _mm256_insert_epi32(mask[1], ~0, 1);
    mask[3] = _mm256_insert_epi32(mask[2], ~0, 2);
    mask[4] = _mm256_insert_epi32(mask[3], ~0, 3);
    mask[5] = _mm256_insert_epi32(mask[4], ~0, 4);
    mask[6] = _mm256_insert_epi32(mask[5], ~0, 5);
    mask[7] = _mm256_insert_epi32(mask[6], ~0, 6);
  }
  // 必: kが存在しないことが保証されている
  void insert(int k){
    k += node_count;
    while(k--){
      int dir = k & split_mod;
      k >>= split_div;
      val[k] = _mm256_add_epi32(val[k], one[dir]);
    }
  }
  // 必: kが存在することが保証されている
  void erase(int k){
    k += node_count;
    while(k--){
      int dir = k & split_mod;
      k >>= split_div;
      val[k] = _mm256_sub_epi32(val[k], one[dir]);
    }
  }
  // [0, r)の要素数
  i256 rank(int r){
    i256 sum = _mm256_set1_epi64x(0);
    r += node_count;
    while(r--){
      int dir = r & split_mod;
      r >>= split_div;
      sum = _mm256_add_epi32(sum, _mm256_and_si256(val[r], mask[dir]));
    }
    return sum;
  }
  // [0, r)の要素数を返しつつrを追加
  // 必: rが存在しないことが保証されている
  i256 rank_insert(int r){
    i256 sum = _mm256_set1_epi64x(0);
    r += node_count;
    while(r--){
      int dir = r & split_mod;
      r >>= split_div;
      sum = _mm256_add_epi32(sum, _mm256_and_si256(val[r], mask[dir]));
      val[r] = _mm256_add_epi32(val[r], one[dir]);
    }
    return sum;
  }
  // [0, r)の要素数を返しつつrを削除
  // 必: rが存在することが保証されている
  i256 rank_erase(int r){
    i256 sum = _mm256_set1_epi64x(0);
    r += node_count;
    while(r--){
      int dir = r & split_mod;
      r >>= split_div;
      sum = _mm256_add_epi32(sum, _mm256_and_si256(val[r], mask[dir]));
      val[r] = _mm256_sub_epi32(val[r], one[dir]);
    }
    return sum;
  }
  long long sum32(__m256i x){
    return (long long)_mm256_extract_epi32(x, 0) + _mm256_extract_epi32(x, 1) + _mm256_extract_epi32(x, 2) + _mm256_extract_epi32(x, 3) +
                      _mm256_extract_epi32(x, 4) + _mm256_extract_epi32(x, 5) + _mm256_extract_epi32(x, 6) + _mm256_extract_epi32(x, 7);
  }
};

template<typename T>
std::vector<long long> offline_static_range_inversion(const std::vector<T> &V, const std::vector<std::pair<int, int>> &Q){
  struct range_inversion_st{
    simd_set b;
    std::vector<int> a;
    int n, all = 0;
    __m256i sum = _mm256_set1_epi64x(0);
    ll allsum = 0;
    range_inversion_st(int n): n(n){}
    void add_left(int i){
      all++;
      sum = _mm256_add_epi32(sum, b.rank_insert(a[i]));
    }
    void del_left(int i){
      all--;
      sum = _mm256_sub_epi32(sum, b.rank_erase(a[i]));
    }
    void add_right(int i){
      allsum += all++;
      sum = _mm256_sub_epi32(sum, b.rank_insert(a[i]));
    }
    void del_right(int i){
      allsum -= --all;
      sum = _mm256_add_epi32(sum, b.rank_erase(a[i]));
    }
  };
  int n = V.size(), q = Q.size();
  range_inversion_st ri(n);
  std::vector<int> tmp(n);
  std::vector<std::pair<T, int>> z;
  for(int i = 0; i < n; i++) z.push_back({V[i], i});
  std::sort(z.begin(), z.end());
  for(int i = 0; i < n; i++) tmp[z[i].second] = i;
  ri.a = tmp;
  ri.b = simd_set();
  mo_algorithm mo(n, ri);
  for(auto [l, r] : Q) mo.insert(l, r);
  std::vector<long long> ans(q);
  while(true){
    auto [i, id] = mo.process();
    if(i == -1) break;
    ans[i] = ri.allsum + ri.b.sum32(ri.sum);
  }
  return ans;
}

int main(){
  int n, q;
  n = io.in();
  q = io.in();
  vector<int> tmp(n);
  vector<pair<int, int>> Q(q);
  range(i, 0, n) tmp[i] = io.in();
  range(i, 0, q) Q[i].first = io.in(), Q[i].second = io.in();
  for(ll ans : offline_static_range_inversion(tmp, Q)) io.out(ans, '\n');
}
*/