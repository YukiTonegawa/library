#ifndef _CEILLOG2_H_
#define _CEILLOG2_H_

unsigned long long ceillog2(unsigned long long x){
  int res = 1;
  while((1ULL << res) < x) res++;
  return res;
}
#endif