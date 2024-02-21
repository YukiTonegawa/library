#ifndef _TERNARY_SEARCH_H_
#define _TERNARY_SEARCH_H_

// 三分探索

// 正整数の範囲で探索するときlが3の倍数かつr == l + 2のとき, l = c1となりwhile(r - l > 1)だと終わらない可能性がある
// (常にc2 < rであるため右側は気にしなくていい)
// f(l), f(l + 1), f(l + 2)を試せばいい
/*
auto f = [s, k](ll x, ll a, ll b){
  return (s + 20 * x + a) * (k - 4 * x + b);
};
auto search = [s, k, f](ll a, ll b){
  ll l = 0, r = k / 4 + 1;
  while(r - l > 2){
    ll c1 = (l * 2 + r) / 3, c2 = (l + r * 2) / 3;
    if(f(c1, a, b) > f(c2, a, b)) r = c2; // 最大値を求める -> 小さい方を切り捨てる
    else l = c1;
  }
  ll ret = max(f(l, a, b), f(l + 1, a, b));
  return ret;
};
*/
#endif
