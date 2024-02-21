#ifndef _ROLLING_HASH_H_
#define _ROLLING_HASH_H_
#include <vector>
#include "string_basic.hpp"
#include "../misc/random_number.hpp"
#include "../math/integer.hpp"
#include "../math/prime.hpp"

// 要素はモンゴメリ表現に変換せずに演算する
struct montgomery_rolling_hash{
  using ull = unsigned long long;
  using ulll = __uint128_t;
  static ull r, rinv, mod;
  static std::vector<ull> rpow, rinvpow, run_inv;
  static montgomery_reduction_64bit mr;
  
  // run関数をO(1)で使いたい場合はuse_run = true
  static void init(int max_n, bool use_run = false){
    assert(0 < max_n || rpow.empty()); // 2回呼び出されると壊れる
    mod = generate_random_prime(1ULL << 58, 1ULL << 63); // 小さすぎず, 足してもullでオーバーフローしない
    mr.set_mod(mod);

    r = random_number(1, mod - 1);
    ull g;
    std::tie(g, rinv) = inv_ext_gcd(r, mod);
    assert(g == 1); // r != 0 && is_prime(mod) -> gcd(r, mod) == 1
    r = mr.generate(r);
    rinv = mr.generate(rinv);
    rpow.resize(max_n + 1);
    rinvpow.resize(max_n + 1);
    run_inv.resize(max_n + 1);
    rpow[0] = rinvpow[0] = run_inv[0] = mr.generate(1);
    ull one = mr.generate(1);

    for(int i = 1; i <= max_n; i++){
      rpow[i] = mr.mul_safe(rpow[i - 1], r);
      rinvpow[i] = mr.mul_safe(rinvpow[i - 1], rinv);
      if(use_run) run_inv[i] = mr.pow_safe(mr.sub_safe(rpow[i], one), mod - 2); // 1 / (rpow[i] - 1)のモンゴメリ表現, i != 0
    }
  }
  template<typename T>
  static ull __mod(T x){
    x %= mod;
    return x < 0 ? mod + x : x;
  }
  static ull hash(const std::string &s){
    assert(s.size() <= rpow.size());
    ull res = 0;
    for(int i = 0; i < s.size(); i++) res = mr.add(res, mr.mul(rpow[i], s[i]));
    return mr.fix(res);
  }
  // 要素aが負の場合やaがmod以上の場合, aとa + modが同一視される
  template<typename T>
  static ull hash(const std::vector<T> &s){
    assert(s.size() <= rpow.size());
    ull res = 0;
    for(int i = 0; i < s.size(); i++) res = mr.add(res, mr.mul(rpow[i], __mod(s[i])));
    return mr.fix(res);
  }
  // table[i] = [0, i)のハッシュテーブル
  template<typename T>
  static std::vector<ull> hash_table(const std::vector<T> &s){
    assert(s.size() <= rpow.size());
    std::vector<ull> res(s.size() + 1);
    res[0] = 0;
    for(int i = 0; i < s.size(); i++){
      res[i + 1] = mr.add_safe(res[i], mr.mul(rpow[i], __mod(s[i])));
    }
    return res;
  }
  // xがk文字目にあった時のハッシュ
  template<typename T>
  static ull get(int k, T x){
    assert(k < rpow.size());
    return mr.mul_safe(rpow[k], __mod(x));
  }
  static ull mul(ull a, ull b){
    return mr.mul_safe(a, b);
  }
  static ull add(ull a, ull b){
    return mr.add_safe(a, b);
  }
  static ull sub(ull a, ull b){
    return mr.sub_safe(a, b);
  }
  // a^b
  static ull pow(ull a, ull b){
    return mr.pow_safe(a, b);
  }
  // 長さalenのハッシュがaの文字の末尾にハッシュがbの文字を結合する
  // O(1)またはO(log alen)
  static ull concat(ull a, ull alen, ull b){
    if(alen < rpow.size()) return add(a, mul(rpow[alen], b));
    return add(a, mul(pow(r, alen), b));
  }
  // 長さalenのハッシュがaの文字をk回繰り返す
  // a + a * rpow[alen] + a * rpow[2 * alen] ...
  // = a * (rpow[k * alen] - 1) / (rpow[alen] - 1)
  static ull run(ull a, ull alen, ull k){
    if(alen == 0) return 0;
    if(k <= 1) return k ? a : 0;
    __uint128_t K = k * alen;
    if(K < (int)rpow.size()) return mul(a, mul(sub(rpow[K], run_inv[0]), run_inv[alen]));

    // rpow[alen]
    unsigned long long r_alen = (alen < rpow.size() ? rpow[alen] : mr.pow_safe(r, alen));

    // rpow[k * alen]
    unsigned long long r_K = mr.pow_safe(r_alen, k);

    // 1 / (rpow[alen] - 1)
    r_alen = mr.pow_safe(sub(r_alen, run_inv[0]), mod - 2);
    return mul(a, mul(sub(r_K, run_inv[0]), r_alen));
  }
};
unsigned long long montgomery_rolling_hash::r = 0;
unsigned long long montgomery_rolling_hash::rinv = 0;
unsigned long long montgomery_rolling_hash::mod = 0;
std::vector<unsigned long long> montgomery_rolling_hash::rpow;
std::vector<unsigned long long> montgomery_rolling_hash::rinvpow;
std::vector<unsigned long long> montgomery_rolling_hash::run_inv;
montgomery_reduction_64bit montgomery_rolling_hash::mr;
// モンゴメリ乗算の方が早い
/* 
struct rolling_hash61{
  using ull = unsigned long long;
  static constexpr ull mod = (1UL << 61) - 1;
  static constexpr ull mod4 = (mod << 2);
  static constexpr ull mask30 = (1UL << 30) - 1;
  static constexpr ull mask31 = (1UL << 31) - 1;
  static constexpr ull mask61 = mod;

  static ull r, rinv;
  static std::vector<ull> rpow, rinvpow, run_inv;

  template<typename T>
  static ull __mod(T x){
    x %= mod;
    return x < 0 ? mod + x : x;
  }
  static ull calc_mod(ull a){
    ull au = a >> 61, ad = a & mask61;
    ull b = au + ad;
    return b >= mod ? b - mod : b;
  }
  static ull add(ull a, ull b){
    ull c = a + b;
    return c >= mod ? c - mod : c;
  }
  static ull sub(ull a, ull b){
    return a < b ? a + mod - b : a - b;
  }
  static ull mul(ull a, ull b){
    ull au = a >> 31, ad = a & mask31;
    ull bu = b >> 31, bd = b & mask31;
    ull c = ad * bu + au * bd;
    ull cu = c >> 30, cd = c & mask30;
    return calc_mod(au * bu * 2 + cu + (cd << 31) + ad * bd);
  }
  // a^b
  static ull pow(ull a, ull b){
    ull c = 1;
    while(b){
      if(b & 1) c = mul(c, a);
      a = mul(a, a);
      b >>= 1;
    }
    return c;
  }
  
  static void init(int max_n, bool use_run = false){
    assert(0 < max_n || rpow.empty()); // 2回呼び出されると壊れる
    r = random_number(2, mod - 1);
    ull g;
    std::tie(g, rinv) = inv_ext_gcd(r, mod);
    rpow.resize(max_n + 1);
    rinvpow.resize(max_n + 1);
    run_inv.resize(max_n + 1);
    rpow[0] = rinvpow[0] = run_inv[0] = 1;
    for(int i = 1; i <= max_n; i++){
      rpow[i] = mul(rpow[i - 1], r);
      rinvpow[i] = mul(rinvpow[i - 1], rinv);
      if(use_run) run_inv[i] = pow(sub(rpow[i], 1), mod - 2); // 1 / (rpow[i] - 1)
    }
  }
  static ull hash(const std::string &s){
    assert(s.size() <= rpow.size());
    ull res = 0;
    for(int i = 0; i < s.size(); i++) res = add(res, mul(rpow[i], s[i]));
    return res;
  }
  // 要素aが負の場合やaがmod以上の場合, aとa + modが同一視される
  template<typename T>
  static ull hash(const std::vector<T> &s){
    assert(s.size() <= rpow.size());
    ull res = 0;
    for(int i = 0; i < s.size(); i++) res = add(res, mul(rpow[i], __mod(s[i])));
    return res;
  }
  // table[i] = [0, i)のハッシュテーブル
  template<typename T>
  static std::vector<ull> hash_table(const std::vector<T> &s){
    assert(s.size() <= rpow.size());
    std::vector<ull> res(s.size() + 1);
    res[0] = 0;
    for(int i = 0; i < s.size(); i++) res[i + 1] = add(res[i], mul(rpow[i], __mod(s[i])));
    return res;
  }
  // xがk文字目にあった時のハッシュ
  template<typename T>
  static ull get(int k, T x){
    assert(k < rpow.size());
    return mul(rpow[k], __mod(x));
  }

  // 長さalenのハッシュがaの文字の末尾にハッシュがbの文字を結合する
  // O(1)またはO(log alen)
  static ull concat(ull a, ull alen, ull b){
    if(alen < rpow.size()) return add(a, mul(rpow[alen], b));
    return add(a, mul(pow(r, alen), b));
  }
  // 長さalenのハッシュがaの文字をk回繰り返す
  // a + a * rpow[alen] + a * rpow[2 * alen] ...
  // = a * (rpow[k * alen] - 1) / (rpow[alen] - 1)
  static ull run(ull a, ull alen, ull k){
    if(alen == 0) return 0;
    if(k <= 1) return k ? a : 0;
    __uint128_t K = k * alen;
    if(K < (int)rpow.size()) return mul(a, mul(sub(rpow[K], run_inv[0]), run_inv[alen]));
    // rpow[alen]
    unsigned long long r_alen = (alen < rpow.size() ? rpow[alen] : pow(r, alen));
    // rpow[k * alen]
    unsigned long long r_K = pow(r_alen, k);
    // 1 / (rpow[alen] - 1)
    r_alen = pow(sub(r_alen, run_inv[0]), mod - 2);
    return mul(a, mul(sub(r_K, run_inv[0]), r_alen));
  }
};
unsigned long long rolling_hash61::r = 0;
unsigned long long rolling_hash61::rinv = 0;
std::vector<unsigned long long> rolling_hash61::rpow;
std::vector<unsigned long long> rolling_hash61::rinvpow;
std::vector<unsigned long long> rolling_hash61::run_inv;

*/

struct rolling_hash_string{
  using ull = unsigned long long;
  using rh = montgomery_rolling_hash;
private:
  int n;
  std::vector<ull> sum;
public:
  template<typename T>
  rolling_hash_string(const std::vector<T> &_v): n(_v.size()), sum(rh::hash_table(_v)){}
  int size() const{
    return n;
  }
  ull get(int k) const{
    assert(k < n);
    return rh::mul(rh::sub(sum[k + 1], sum[k]), rh::rinvpow[k]);
  }
  template<typename T>
  void push_back(T x){
    ull h = rh::add(sum.back(), rh::mul(rh::__mod(x), rh::rpow[n]));
    sum.push_back(h);
    n++;
  }
  void pop_back(){
    assert(n);
    sum.pop_back();
    n--;
  }
  // 全体のハッシュ
  ull hash_all() const{
    return sum.back();
  }
  // [l, r)のハッシュ
  // l == rなら0
  ull hash_range(int l, int r) const{
    assert(0 <= l && l <= r && r <= n);
    if(!l) return sum[r];
    return rh::mul(rh::sub(sum[r], sum[l]), rh::rinvpow[l]);
  }
  // 結合後のサイズがmontgomeryテーブルのサイズを超えない
  void concat(const rolling_hash_string &b){
    int m = b.size();
    assert(n + m < rh::rpow.size());
    unsigned long long x = rh::rpow[n], y = sum.back();
    for(int i = 1; i <= m; i++) sum.push_back(rh::add(y, rh::mul(b.sum[i], x)));
    n += m;
  }
  // aの[l1...]とbの[l2...]のlcp
  static int lcp(const rolling_hash_string &a, const rolling_hash_string &b, int l1, int l2){
    int L = 0, R = std::min(a.size() - l1, b.size() - l2) + 1;
    while(R - L > 1){
      int mid = (L + R) >> 1;
      if(a.hash_range(l1, l1 + mid) == b.hash_range(l2, l2 + mid)) L = mid;
      else R = mid;
    }
    return L;
  }
  // これの[l1...]とbの[l2...]のlcp
  int lcp(const rolling_hash_string &b, int l1, int l2) const{
    return lcp(*this, b, l1, l2);
  }
};

struct dynamic_rolling_hash_string{
  using ull = unsigned long long;
  using rh = montgomery_rolling_hash;
private:
  int n, h;
  std::vector<ull> sum, point;
  void __init(std::vector<ull> &v){
    sum.resize(1, 0);
    sum.insert(sum.begin() + 1, v.begin(), v.end());
    for(int i = 1; i <= n; i++){
      int nxt = i + (i & (-i));
      if(nxt <= n) sum[nxt] = rh::add(sum[nxt], sum[i]);
    }
    if(n) h = 31 - __builtin_clz(n);
  }
  void __update(int k, ull x){
    for(int i = k + 1; i <= n; i += (i & (-i))){
      sum[i] = rh::add(sum[i], x);
    }
  }
  ull __query(int r) const{
    ull res = 0;
    for(int k = r; k > 0; k -= (k & (-k))){
      res = rh::mr.add(res, sum[k]);
    }
    return rh::mr.fix(res);
  }
  ull __query(int l, int r) const{
    ull a = __query(r), b = __query(l);
    return rh::sub(a, b);
  }
  ull __query_all() const{
    return __query(n);
  }
public:
  template<typename T>
  dynamic_rolling_hash_string(const std::vector<T> &_v): n(_v.size()){
    std::vector<ull> table(n);
    point.resize(n);
    for(int i = 0; i < n; i++){
      table[i] = rh::get(i, _v[i]);
      point[i] = rh::__mod(_v[i]);
    }
    __init(table);
  }
  int size() const{
    return n;
  }
  template<typename T>
  void set(int k, T x){
    ull _x = rh::__mod(x), _y = get(k);
    point[k] = _x;
    __update(k, rh::mul(rh::sub(_x, _y), rh::rpow[k]));
  }
  ull get(int k) const{
    assert(0 <= k && k < n);
    return point[k];
  }
  // 全体のハッシュ
  ull hash_all() const{
    return __query_all();
  }
  // [0, r)のハッシュ
  ull hash_range(int r) const{
    return __query(r);
  }
  // [l, r)のハッシュ
  // l == rなら0
  ull hash_range(int l, int r) const{
    assert(0 <= l && l <= r && r <= n);
    if(l == 0) return __query(r);
    return rh::mul(__query(l, r), rh::rinvpow[l]);
  }
  static int lcp(const dynamic_rolling_hash_string &a, const dynamic_rolling_hash_string &b, int l1, int l2){
    int v = 1 << a.h, H = a.h;
    ull s = rh::sub(0, a.__query(l1));
    ull hash_b = b.__query(l2);
    while(H--){
      int len = v - l1;
      if(a.n < v || l2 + len > b.size()){
        v -= 1 << H;
        continue;
      }
      ull hash_tmp = rh::add(s, a.sum[v]);
      if(v <= l1){
        s = hash_tmp;
        v += 1 << H;
        continue;
      }
      bool f = rh::mul(rh::sub(b.__query(l2 + len), hash_b), rh::rinvpow[l2]) != rh::mul(hash_tmp, rh::rinvpow[l1]);
      if(f){
        v -= 1 << H;
      }else{
        s = hash_tmp;
        v += 1 << H;
      }
    }
    if(v == a.n + 1) return a.n - l1;
    int len = v - l1;
    bool f = l2 + len > b.size() || b.hash_range(l2, l2 + len) != rh::mul(rh::add(s, a.sum[v]), rh::rinvpow[l1]);
    return f ? v - 1 - l1 : v - l1;
  }
  // len(a) == len(b)かつ0からのlcp O(loglen(a))
  static int lcp_same_size(const dynamic_rolling_hash_string &a, const dynamic_rolling_hash_string &b){
    assert(a.size() == b.size());
    int v = 1 << a.h, H = a.h;
    while(H--){
      if(a.n < v) v -= 1 << H;
      else if(a.sum[v] != b.sum[v]) v -= 1 << H;
      else v += 1 << H;
    }
    if(v == a.n + 1) return a.n;
    return a.sum[v] != b.sum[v] ? v - 1 : v;
  }
  static int lcp(const dynamic_rolling_hash_string &a, const rolling_hash_string &b, int l1, int l2){
    int v = 1 << a.h, H = a.h;
    ull s = rh::sub(0, a.__query(l1));
    while(H--){
      int len = v - l1;
      if(a.n < v || l2 + len > b.size()){
        v -= 1 << H;
        continue;
      }
      ull hash_tmp = rh::add(s, a.sum[v]);
      if(v <= l1){
        s = hash_tmp;
        v += 1 << H;
        continue;
      }
      bool f = b.hash_range(l2, l2 + len) != rh::mul(hash_tmp, rh::rinvpow[l1]);
      if(f){
        v -= 1 << H;
      }else{
        s = hash_tmp;
        v += 1 << H;
      }
    }
    if(v == a.n + 1) return a.n - l1;
    int len = v - l1;
    bool f = l2 + len > b.size() || b.hash_range(l2, l2 + len) != rh::mul(rh::add(s, a.sum[v]), rh::rinvpow[l1]);
    return f ? v - 1 - l1 : v - l1;
  }
  static int lcp(const rolling_hash_string &a, const dynamic_rolling_hash_string &b, int l1, int l2){
    return lcp(b, a, l2, l1);
  }
};

struct dynamic_rolling_hash_string_persistent{
  using ull = unsigned long long;
  using rh = montgomery_rolling_hash;
private:
  struct node{
    node *l, *r;
    int sz;
    ull val;
    node(ull x): l(nullptr), r(nullptr), sz(1), val(x){}
    node(node *l, node *r): l(l), r(r), sz(l->sz + r->sz), val(rh::add(l->val, rh::mul(r->val, rh::rpow[l->sz]))){}
  };
  template<typename T>
  static node *__build(const std::vector<T> &v, int l, int r){
    if(l + 1 == r) return new node(v[l]);
    int mid = (l + r) / 2;
    return new node(__build(v, l, mid), __build(v, mid, r));
  }
  static node *__set(node *v, int k, ull x, int l, int r){
    if(r - l == 1) return new node(x);
    int mid = (l + r) / 2;
    if(k < mid){
      return new node(__set(v->l, k, x, l, mid), v->r);
    }else{
      return new node(v->l, __set(v->r, k, x, mid, r));
    }
  }
  static ull __get(node *v, int k, int l, int r){
    if(r - l == 1) return v->val;
    int mid = (l + r) / 2;
    if(k < mid) return __get(v->l, k, l, mid);
    else return __get(v->r, k, mid, r);
  }
  static ull __hash_range(node *v, int b){
    int l = 0, r = v->sz;
    b = std::min(b, r);
    int lsz = 0;
    ull res = 0;
    while(v){
      if(b <= l) return res;
      if(b >= r) return rh::add(res, rh::mul(v->val, rh::rpow[lsz]));
      int mid = (l + r) / 2;
      if(b <= mid){
        r = mid;
        v = v->l;
      }else{
        res = rh::add(res, rh::mul(v->l->val, rh::rpow[lsz]));
        lsz += mid - l;
        l = mid;
        v = v->r;
      }
    }
    return res;
  }
  static ull __hash_range(node *v, int a, int b, int l, int r){
    if(!v || b <= l || r <= a) return 0;
    if(a <= l && r <= b) return v->val;
    int mid = (l + r) / 2;
    return rh::add(__hash_range(v->l, a, b, l, mid), rh::mul(__hash_range(v->r, a, b, mid, r), rh::rpow[mid - l]));
  }
public:
  dynamic_rolling_hash_string_persistent(){}
  template<typename T>
  static node *build(const std::vector<T> &v){
    int n = v.size();
    return __build(v, 0, n);
  }
  static int size(node *v){
    return v ? v->sz : 0;
  }
  static node *set(node *v, int k, ull x){
    assert(k < size(v));
    return __set(v, k, x, 0, v->sz);
  }
  static ull get(node *v, int k){
    assert(k < size(v));
    return __get(v, k, 0, v->sz);
  }
  static ull hash_range(node *v, int r){
    assert(0 <= r && r <= size(v));
    if(r == 0) return 0;
    return __hash_range(v, r);
  }
  static ull hash_range(node *v, int l, int r){
    if(l >= r) return 0;
    assert(0 <= l && r <= size(v));
    return __hash_range(v, l, r, 0, v->sz);
  }
  static int __lcp(node *v, const rolling_hash_string &b, int l1, int l2){
    if(!v || l2 >= b.size()) return 0;
    if(l1 == 0 && l2 + v->sz <= b.size() && v->val == b.hash_range(l2, l2 + v->sz)) return v->sz;
    int szl = v->l ? v->l->sz : 0;
    if(szl <= l1) return __lcp(v->r, b, l1 - szl, l2);
    int left_cross = szl - l1, L = __lcp(v->l, b, l1, l2);
    if(L != left_cross) return L;
    l2 += left_cross;
    return L + __lcp(v->r, b, 0, l2);
  }
  static int __lcp(node *v, const dynamic_rolling_hash_string &b, int l1, int l2, int &oksz, ull &ok){
    if(!v || l2 >= b.size()) return 0;
    if(l1 == 0 && l2 + v->sz <= b.size()){
      ull rh = rh::mul(rh::sub(b.hash_range(l2 + v->sz), ok), rh::rinvpow[l2]);
      if(v->val == rh){
        ok = rh::add(ok, rh::mul(v->val, rh::rpow[oksz]));
        oksz += v->sz;
        return v->sz;
      }
    }
    int szl = v->l ? v->l->sz : 0;
    if(szl <= l1) return __lcp(v->r, b, l1 - szl, l2, oksz, ok);
    int left_cross = szl - l1, L = __lcp(v->l, b, l1, l2, oksz, ok);
    if(L != left_cross) return L;
    l2 += left_cross;
    return L + __lcp(v->r, b, 0, l2, oksz, ok);
  }
  static int __lcp(node *v, node *u, int l1, int l2, int &oksz, ull &ok){
    if(!v || l2 >= size(u)) return 0;
    if(l1 == 0 && l2 + v->sz <= size(u)){
      ull rh = rh::mul(rh::sub(hash_range(u, l2 + v->sz), ok), rh::rinvpow[l2]);
      if(v->val == rh){
        ok = rh::add(ok, rh::mul(v->val, rh::rpow[oksz]));
        oksz += v->sz;
        return v->sz;
      }
    }
    int szl = v->l ? v->l->sz : 0;
    if(szl <= l1) return __lcp(v->r, u, l1 - szl, l2, oksz, ok);
    int left_cross = szl - l1, L = __lcp(v->l, u, l1, l2, oksz, ok);
    if(L != left_cross) return L;
    l2 += left_cross;
    return L + __lcp(v->r, u, 0, l2, oksz, ok);
  }
  static int lcp(node *v, const rolling_hash_string &b, int l1, int l2){
    if(!v) return 0;
    return __lcp(v, b, l1, l2);
  }
  static int lcp(const rolling_hash_string &b, node *v, int l1, int l2){
    if(!v) return 0;
    return __lcp(v, b, l2, l1);
  }
  static int lcp(node *v, const dynamic_rolling_hash_string &b, int l1, int l2){
    int oksz = l2;
    ull ok = b.hash_range(l2);
    return __lcp(v, b, l1, l2, oksz, ok);
  }
  static int lcp(const dynamic_rolling_hash_string &a, node *v, int l1, int l2){
    return lcp(v, a, l2, l1);
  }
  static int lcp(node *v, node *u, int l1, int l2){
    int oksz = l2;
    ull ok = hash_range(u, 0, l2);
    return __lcp(v, u, l1, l2, oksz, ok);
  }
};
#endif