#ifndef _WAVELET_MATRIX2D_H_
#define _WAVELET_MATRIX2D_H_
#include "wavelet_matrix.hpp"

template<typename Idx, typename Val>
struct wavelet_matrix2d{
private:
  int n, h;
  std::vector<Idx> X, Y;
  std::vector<Val> Z;
  using __bitvector = bitvector_memory;
  std::vector<__bitvector> bv;
  std::vector<wavelet_matrix_arbitrary_radix<8>> wm;
  inline int lb_x(Idx c){return std::lower_bound(X.begin(), X.end(), c) - X.begin();}
  inline int lb_y(Idx c){return std::lower_bound(Y.begin(), Y.end(), c) - Y.begin();}
  inline int lb_z(Val c){return std::lower_bound(Z.begin(), Z.end(), c) - Z.begin();}
  void build(std::vector<std::pair<int, int>> v){
    std::vector<std::pair<int, int>> tmp(v.size());
    wm.resize(h + 1);
    bv.resize(h + 1);
    std::queue<std::tuple<int, int, int>> q; // {h, l, r}
    q.push({h, 0, v.size()});
    std::vector<bool> bits(n);
    std::vector<int> y(n);
    while(!q.empty()){
      auto [b, l, r] = q.front();
      q.pop();
      if(b){
        int lcnt = 0, rcnt = 0;
        for(int i = l; i < r; i++){
          bool dir = (v[i].first >> (b - 1)) & 1; 
          bits[i] = dir;
          y[i] = v[i].second;
          if(dir) tmp[rcnt++] = v[i];
          else v[l + (lcnt++)] = v[i];
        }
        int mid = l + lcnt;
        for(int i = 0; i < rcnt; i++) v[mid + i] = tmp[i];
        if(mid > l) q.push({b - 1, l, mid});
        if(r > mid) q.push({b - 1, mid, r});
        if(r == n){
          wm[b] = wavelet_matrix_arbitrary_radix<8>(y, Y.size());
          bv[b] = __bitvector(bits);
        }
      }else{
        for(int i = l; i < r; i++) y[i] = v[i].second;
        if(r == n) wm[b] = wavelet_matrix_arbitrary_radix<8>(y, Y.size());
      }
    }
  }
  int __range_freq(int d, int lx, int rx, int Lx, int Rx, int lz, int rz, int Lz, int Rz, int ly, int ry){
    if(lx == rx || Lx == Rx || rz <= Lz || Rz <= lz) return 0;
    else if(lz <= Lz && Rz <= rz) return wm[d].range_freq(lx, rx, ly, ry);
    int L0 = bv[d].rank0(Lx), R0 = bv[d].rank0(Rx) - L0;
    int l0 = bv[d].rank0(lx) - L0, r0 = bv[d].rank0(rx) - L0;
    return __range_freq(d - 1, Lx + l0, Lx + r0, Lx, Lx + R0, lz, rz, Lz, (Lz + Rz) / 2, ly, ry) +
            __range_freq(d - 1, lx + (R0 - l0), rx + (R0 - r0), Lx + R0, Rx, lz, rz, (Lz + Rz) / 2, Rz, ly, ry);
  }
  int __range_kth_smallest(int lx, int rx, int ly, int ry, int k){
    if(lx >= rx || ly >= ry) return -1;
    int Lx = 0, Rx = n, ret = 0;
    for(int i = h; i > 0; i--){
      int L0 = bv[i].rank0(Lx), R0 = bv[i].rank0(Rx) - L0;
      int l0 = bv[i].rank0(lx) - L0, r0 = bv[i].rank0(rx) - L0;
      int left_cnt = wm[i - 1].range_freq(Lx + l0, Lx + r0, ly, ry);
      if(left_cnt <= k){
        k -= left_cnt;
        ret += 1 << (i - 1);
        lx += R0 - l0, rx += R0 - r0, Lx += R0;
      }else{
        lx = Lx + l0, rx = Lx + r0, Rx = Lx + R0;
      }
    }
    if(wm[0].range_freq(lx, rx, ly, ry) <= k) return -1;
    return ret;
  }
  int __range_kth_largest(int lx, int rx, int ly, int ry, int k){
    if(lx >= rx || ly >= ry) return -1;
    int Lx = 0, Rx = n, ret = 0;
    for(int i = h; i > 0; i--){
      int L0 = bv[i].rank0(Lx), R0 = bv[i].rank0(Rx) - L0;
      int l0 = bv[i].rank0(lx) - L0, r0 = bv[i].rank0(rx) - L0;
      int right_cnt = wm[i - 1].range_freq(lx + (R0 - l0), rx + (R0 - r0), ly, ry);
      if(right_cnt <= k){
        k -= right_cnt;
        lx = Lx + l0, rx = Lx + r0, Rx = Lx + R0;
      }else{
        ret += 1 << (i - 1);
        lx += R0 - l0, rx += R0 - r0, Lx += R0;
      }
    }
    if(wm[0].range_freq(lx, rx, ly, ry) <= k) return -1;
    return ret;
  }
public:
  using point = std::tuple<Idx, Idx, Val>;
  wavelet_matrix2d(std::vector<point> _v): n(_v.size()){
    if(n == 0) return;
    std::vector<std::pair<int, int>> v(n); // {z, y}
    std::sort(_v.begin(), _v.end(), [](point &a, point &b){
      return std::get<0>(a) < std::get<0>(b);
    });
    for(int i = 0; i < n; i++){
      auto [x, y, z] = _v[i];
      X.push_back(x);
      Y.push_back(y);
      Z.push_back(z);
    }
    std::sort(Z.begin(), Z.end());
    Z.erase(std::unique(Z.begin(), Z.end()), Z.end());
    std::sort(Y.begin(), Y.end());
    Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
    for(int i = 0; i < n; i++){
      auto [x, y, z] = _v[i];
      v[i] = {lb_z(z), lb_y(y)};
    }
    h = 0;
    while((1 << h) < Z.size()) h++;
    build(v);
  }
  // [lx, rx) × [ly, ry)の[lz, rz)の要素数
  int range_freq(Idx lx, Idx rx, Idx ly, Idx ry, Val lz, Val rz){
    int lxi = lb_x(lx), rxi = lb_x(rx);
    int lyi = lb_y(ly), ryi = lb_y(ry);
    int lzi = lb_z(lz), rzi = lb_z(rz);
    if(lxi >= rxi || lyi >= ryi || lzi >= rzi) return 0;
    return __range_freq(h, lxi, rxi, 0, n, lzi, rzi, 0, 1 << h, lyi, ryi);
  }
  // [lx, rx) × [ly, ry)のk番目に小さい要素, ない場合は-1
  Val range_kth_smallest(Idx lx, Idx rx, Idx ly, Idx ry, int k){
    int lxi = lb_x(lx), rxi = lb_x(rx);
    int lyi = lb_y(ly), ryi = lb_y(ry);
    if(lxi >= rxi || lyi >= ryi) return -1;
    int ret = __range_kth_smallest(lxi, rxi, lyi, ryi, k);
    return ret == -1 ? -1 : Z[ret];
  }
  // [lx, rx) × [ly, ry)のk番目に大きい要素, ない場合は-1
  Val range_kth_largest(Idx lx, Idx rx, Idx ly, Idx ry, int k){
    int lxi = lb_x(lx), rxi = lb_x(rx);
    int lyi = lb_y(ly), ryi = lb_y(ry);
    if(lxi >= rxi || lyi >= ryi) return -1;
    int ret = __range_kth_largest(lxi, rxi, lyi, ryi, k);
    return ret == -1 ? -1 : Z[ret];
  }
};
#endif