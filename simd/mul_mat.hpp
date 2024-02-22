/* mul_mat_simdを使う場合mainファイルの行頭に書く
  #pragma GCC optimize ("O3")
  #pragma GCC target ("avx2")
  #include <immintrin.h>
*/
/*
struct montgomery_parallel{
  constexpr static int B = 31;
  long long _MOD;
  __m256i MOD, MOD_CMP, R2, INV, MASK;
  __m256i MOD2_CMP, MOD2, ZERO_CMP = _mm256_set1_epi64x(-1);
  void init(){
    long long _R2, _INV = 0, t = 0, s = 1;
    for(int i = 0; i < B; i++){
      if((t & 1) == 0) t += _MOD, _INV += s;
      t >>= 1, s <<= 1;
    }
    _R2 = (1ULL << B) % _MOD;
    _R2 = (_R2 * _R2) % _MOD;
    MOD = _mm256_set1_epi64x(_MOD);
    MOD_CMP = _mm256_set1_epi64x(_MOD - 1);
    MOD2 = _mm256_set1_epi64x(2 * _MOD);
    MOD2_CMP = _mm256_set1_epi64x(2 * _MOD - 1);
    R2 = _mm256_set1_epi64x(_R2);
    INV = _mm256_set1_epi64x(_INV);
    MASK = _mm256_set1_epi64x((1ULL << B) - 1);
  }
  montgomery_parallel(int _MOD): _MOD(_MOD){init();}

  // [0, 2MOD)
  inline __m256i reduce(__m256i x){
    __m256i y = _mm256_mul_epi32(x, INV);
    y = _mm256_and_si256(y, MASK);
    y = _mm256_mul_epi32(y, MOD);
    y = _mm256_add_epi64(y, x);
    y = _mm256_srli_epi64(y, B);
    return y;
  }
  inline __m256i mul_raw(__m256i a, __m256i b){
    return reduce(_mm256_mul_epi32(reduce(_mm256_mul_epi32(a, b)), R2));
  }
  // [0, 2MOD) -> [0, MOD)
  inline __m256i fix(__m256i x){
    return _mm256_sub_epi64(x, _mm256_and_si256(MOD, _mm256_cmpgt_epi64(x, MOD_CMP)));
  }
  // [0, MOD) -> モンゴメリ表現
  inline __m256i generate(__m256i x){
    return reduce(_mm256_mul_epi32(x, R2));
  }
  // モンゴメリ表現同士の積, [0, 2MOD)
  inline __m256i mul(__m256i ma, __m256i mb){
    return reduce(_mm256_mul_epi32(ma, mb));
  }
  // モンゴメリ表現同士の和, [0, 2MOD)
  inline __m256i add(__m256i ma, __m256i mb){
    ma = _mm256_add_epi64(ma, mb); // [0, 4MOD)
    return _mm256_sub_epi64(ma, _mm256_and_si256(MOD2, _mm256_cmpgt_epi64(ma, MOD2_CMP))); // [0, 2MOD)
  }
  // モンゴメリ表現同士の差, [0, 2MOD)
  inline __m256i sub(__m256i ma, __m256i mb){
    ma = _mm256_sub_epi64(ma, mb); // (-2MOD, 2MOD)
    __m256i neg =  _mm256_andnot_si256(_mm256_cmpgt_epi64(ma, ZERO_CMP), ZERO_CMP);
    return _mm256_add_epi64(ma, _mm256_and_si256(MOD2, neg));
  }
};

template<typename __matrix>
__matrix mul_mat_simd(__matrix &A, __matrix &B){
  using mint = typename __matrix::_mint;
  constexpr int __mod = mint::mod();
  int N = A.n, M = B.n, K = B.m;
  assert(A.m == M);
  montgomery_parallel mr(__mod);
  std::vector<std::vector<__m256i>> Arow(N), Bcol(K);
  std::array<int, 4> tmp;
  for(int i = 0; i < N; i++){
    int k = 0;
    for(int j = 0; j < M / 4; j++){
      long long a, b, c, d;
      a = A[i][k++].val();
      b = A[i][k++].val();
      c = A[i][k++].val();
      d = A[i][k++].val();
      Arow[i].push_back(_mm256_set_epi64x(a, b, c, d));
    }
    if(M % 4){
      tmp.fill(0);
      for(int j = 0; j < M % 4; j++) tmp[j] = A[i][k++].val();
      Arow[i].push_back(_mm256_set_epi64x(tmp[0], tmp[1], tmp[2], tmp[3]));
    }
  }
  for(int i = 0; i < K; i++){
    int k = 0;
    for(int j = 0; j < M / 4; j++){
      long long a, b, c, d;
      a = B[k++][i].val();
      b = B[k++][i].val();
      c = B[k++][i].val();
      d = B[k++][i].val();
      Bcol[i].push_back(_mm256_set_epi64x(a, b, c, d));
    }
    if(M % 4){
      tmp.fill(0);
      for(int j = 0; j < M % 4; j++) tmp[j] = B[k++][i].val();
      Bcol[i].push_back(_mm256_set_epi64x(tmp[0], tmp[1], tmp[2], tmp[3]));
    }
  }
  __matrix ans(N, K, 0);
  
  int M2 = (M + 3) / 4;
  for(int i = 0; i < N; i++) for(int t = 0; t < M2; t++) Arow[i][t] = mr.generate(Arow[i][t]);
  for(int i = 0; i < K; i++) for(int t = 0; t < M2; t++) Bcol[i][t] = mr.generate(Bcol[i][t]);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < K; j++){
      __m256i ans256 = _mm256_setzero_si256();
      for(int t = 0; t < M2; t++){
        ans256 = mr.add(ans256, mr.mul(Arow[i][t], Bcol[j][t]));
      }
      ans256 = mr.reduce(ans256);
      ans[i][j] = mint(ans256[0] + ans256[1] + ans256[2] + ans256[3]);
    }
  }
  return ans;
}
*/