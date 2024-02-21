#ifndef _BERLEKAMP_MASSEY_H_
#define _BERLEKAMP_MASSEY_H_
#include <vector>
// 任意のi>=Lで S[i] = C[1]S[i-1] + C[2]S[i-2]...C[L]S[i-L]を満たすような最小次数のC
// O(N^2)
template<typename mint>
std::vector<mint> berlekamp_massey(const std::vector<mint> &S){
  std::vector<mint> C{1}, B{1};
  int N = S.size();
  int L = 0, m = 1;
  mint b = 1;
  // A(x) <- A(x) - c * B(x)
  auto mul_sub_poly = [](std::vector<mint> &X, const std::vector<mint> &Y, mint z, int yshift) -> void {
    if(X.size() < Y.size() + yshift) X.resize(Y.size() + yshift, 0);
    for(int i = yshift; i < Y.size() + yshift; i++) X[i] -= z * Y[i - yshift];
  };
  for(int n = 0; n < N; n++){
    mint d = 0;
    for(int i = 0; i <= L; i++) d += C[i] * S[n - i];
    if(d == 0){
      m++;
    }else if(2 * L <= n){
      auto T = C;
      mul_sub_poly(C, B, d / b, m);
      L = n + 1 - L;
      B = T;
      b = d;
      m = 1;
    }else{
      mul_sub_poly(C, B, d / b, m);
      m++;
    }
  }
  std::vector<mint> res(L + 1, 0);
  for(int i = 1; i <= L; i++) res[i] = -C[i];
  return res;
}
#endif