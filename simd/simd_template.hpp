#ifndef _SIMD_TEMPLATE_H_
#define _SIMD_TEMPLATE_H_

/* コード上部にこれを貼る
#pragma GCC optimize ("O3")
#pragma GCC target ("avx2")
#include <immintrin.h>
*/

#include <immintrin.h>
#include <vector>
#include <array>

const __m256i allzero = _mm256_setzero_si256();
const __m256i allone  = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);

/*
std::vector<__m256i> pack_int(const std::vector<int> &v){
  static constexpr int block_size = 8;
  int n = v.size();
  int m = (n + block_size - 1) / block_size;
  std::vector<__m256i> res(m);
  for(int i = 0; i < m; i++){
    std::array<int, block_size> tmp;
    for(int j = 0, k = i * block_size; j < block_size; j++, k++){
      tmp[j] = k < n ? v[k] : 0; // 余った分は0埋め
    }
    res[i] = _mm256_set_epi32(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7]);
  }
  return res;
}
std::vector<__m256i> pack_ll(const std::vector<long long> &v){
  static constexpr int block_size = 4;
  int n = v.size();
  int m = (n + block_size - 1) / block_size;
  std::vector<__m256i> res(m);
  for(int i = 0; i < m; i++){
    std::array<int, block_size> tmp;
    for(int j = 0, k = i * block_size; j < block_size; j++, k++){
      tmp[j] = k < n ? v[k] : 0; // 余った分は0埋め
    }
    res[i] = _mm256_set_epi64x(tmp[0], tmp[1], tmp[2], tmp[3]);
  }
  return res;
}
std::vector<int> unpack_int(const std::vector<__m256i> &v){
  static constexpr int block_size = 8;
  int m = v.size();
  std::vector<int> res(m * block_size);
  for(int i = 0; i < m; i++){
    int j = i * block_size;
    res[j    ] = _mm256_extract_epi32(v[i], 0);
    res[j + 1] = _mm256_extract_epi32(v[i], 1);
    res[j + 2] = _mm256_extract_epi32(v[i], 2);
    res[j + 3] = _mm256_extract_epi32(v[i], 3);
    res[j + 4] = _mm256_extract_epi32(v[i], 4);
    res[j + 5] = _mm256_extract_epi32(v[i], 5);
    res[j + 6] = _mm256_extract_epi32(v[i], 6);
    res[j + 7] = _mm256_extract_epi32(v[i], 7);
  }
  return res;
}
std::vector<long long> unpack_ll(const std::vector<__m256i> &v){
  static constexpr int block_size = 4;
  int m = v.size();
  std::vector<long long> res(m * block_size);
  for(int i = 0; i < m; i++){
    int j = i * block_size;
    res[j    ] = _mm256_extract_epi64(v[i], 0);
    res[j + 1] = _mm256_extract_epi64(v[i], 1);
    res[j + 2] = _mm256_extract_epi64(v[i], 2);
    res[j + 3] = _mm256_extract_epi64(v[i], 3);
  }
  return res;
}
*/
#endif