#ifndef _PSEUDO_TREE_H_
#define _PSEUDO_TREE_H_
#include <vector>
#include <cassert>
#include <algorithm>

// 0-indexedのセグメントツリーを模した木
template<typename Idx = int>
struct pseudo_segment_tree{
  static constexpr int bitlen = sizeof(Idx) * 8;
  Idx N, M;
  pseudo_segment_tree(){}
  pseudo_segment_tree(Idx n): N(n){
    M = 1;
    while(M < N) M <<= 1;
  }
  // aの深さ
  int depth(Idx a){
    if(bitlen <= 32) return 31 - __builtin_clz(a + 1);
    return 63 - __builtin_clzll(a + 1);
  }
  // aが表す区間の幅
  Idx width(Idx a){
    return M >> depth(a);
  }
  // aが葉か
  bool is_leaf(Idx a){
    return M - 1 <= a;
  }
  // a, bの最短距離
  Idx dist(Idx a, Idx b){
    return depth(a) + depth(b) - 2 * depth(lca(a, b));
  }
  // aのk個親, 深さを超える場合は-1
  Idx la(Idx a, int k){
    if(depth(a) < k) return -1;
    return ((a + 1) >> k) - 1;
  }
  // lca
  Idx lca(Idx a, Idx b){
    a++, b++;
    int da = depth(a), db = depth(b);
    if(da > db) std::swap(a, b), std::swap(da, db);
    b >>= (db - da);
    if(a == b) return a - 1;
    int msb_diff = (bitlen <= 32 ? 31 - __builtin_clz(a ^ b) : 63 - __builtin_clzll(a ^ b)) + 1;
    return (a >> msb_diff) - 1;
  }
  // aが対応する区間
  std::pair<Idx, Idx> index_to_range(Idx a){
    assert(0 <= a && a < 2 * M - 1);
    int dep = depth(a);
    Idx offset = (a + 1) - ((Idx)1 << dep), wid = M >> dep;
    return std::make_pair(offset * wid, (offset + 1) * wid);
  }
  // 区間[l, r)に対応するノード番号(左が先)
  std::vector<Idx> range_to_index(Idx l, Idx r){
    l = std::max(l, 0), r = std::min(r, N);
    assert(l <= r);
    l += M, r += M;
    std::vector<Idx> left, right;
    while(l < r){
      if(l & 1) left.push_back((l++) - 1);
      if(r & 1) right.push_back((--r) - 1);
      l >>= 1;
      r >>= 1;
    }
    std::reverse(right.begin(), right.end());
    left.insert(left.end(), right.begin(), right.end());
    return left;
  }
  // 葉a( < N) から根まで辿るときのノード番号(底が先)
  std::vector<Idx> leaf_to_root(Idx a){
    assert(0 <= a && a < N);
    a += M - 1;
    std::vector<Idx> ret{a};
    while(a){
      a = (a - 1) >> 1;
      ret.push_back(a);
    }
    return ret;
  }
};


// 0-indexedのk分木を模した木
// 頂点iから ki + 1, ki + 2....ki + kに辺が伸びている(nを超える場合はなし)
// (= 頂点iから (i - 1) / kに辺が伸びている(0からはなし))
template<typename Idx, int k>
struct pseudo_k_ary_tree{
  static constexpr int bitlen = sizeof(Idx) * 8;
  Idx N, M;
  std::vector<Idx> Lelem; // 各深さの最左ノード
  std::vector<Idx> kpow;
  pseudo_k_ary_tree(){}
  pseudo_k_ary_tree(Idx n): N(n){
    assert(n);
    M = 1;
    Lelem.push_back(0);
    while(M < N){
      Lelem.push_back(M);
      // Mは最大でNK程度になり, N, kが大きいとMがオーバーフローする可能性がある
      assert((std::numeric_limits<Idx>::max() - 1) / k >= M);
      M = (M * k + 1);
    }
    Idx p = 1;
    for(int i = 0; i < Lelem.size(); i++){
      kpow.push_back(p);
      p *= k;
    }
  }
  int height(){
    return Lelem.size();
  }
  // aの深さ
  int depth(Idx a){
    int ret = 0;
    while(a){
      a = (a - 1) / k;
      ret++;
    }
    return ret;
  }
  // {aの深さ, aと同じ深さのノードでaより小さいものの数}
  std::pair<int, Idx> index_sibling(Idx a){
    int d = depth(a);
    return {d, a - Lelem[d]};
  }
  // 深さが最も深いノードの数
  Idx num_deepest(){
    return N - Lelem.back();
  }
  // 葉の数
  Idx num_leaf(){
    Idx nd = num_deepest();
    Idx ALLLEAF = M - Lelem.back();
    return nd + (ALLLEAF - nd) / k;
  }
  // aの部分木に含まれる最も深いノードの数
  Idx num_subdeepest(Idx a){
    auto [d, si] = index_sibling(a);
    int hdiff = (int)Lelem.size() - d;
    // 完全k分木ならk ^ (h - 1 - d)個の葉がある
    return std::max(Idx(0), num_deepest() - si * kpow[hdiff - 1]);
  }
  // aの部分木に含まれる葉の数
  Idx num_subleaf(Idx a){
    auto [d, si] = index_sibling(a);
    int hdiff = (int)Lelem.size() - d;
    // 完全k分木ならk ^ (h - 1 - d)個の葉がある
    Idx subdeep = std::max(Idx(0), num_deepest() - si * kpow[hdiff - 1]);
    return subdeep + (kpow[hdiff - 1] - subdeep) / k;
  }
  // aの部分木のサイズ
  Idx num_subtree(Idx a){
    auto [d, si] = index_sibling(a);
    int hdiff = (int)Lelem.size() - d;
    // 完全k分木ならk ^ (h - 1 - d)個の葉がある
    Idx subdeep = std::max(Idx(0), num_deepest() - si * kpow[hdiff - 1]);
    return Lelem[hdiff - 1] + subdeep;
  }
  // aが葉か
  bool is_leaf(Idx a){
    return Lelem.back() <= a;
  }
  Idx dist(Idx a, Idx b){
    return dist2(a, b).second;
  }
  // a, bの{lca, 最短距離}
  std::pair<Idx, Idx> dist2(Idx a, Idx b){
    Idx d = 0;
    while(a != b){
      if(a < b) std::swap(a, b);
      a = (a - 1) / k;
      d++;
    }
    return {a, d};
  }
  // 親, ない場合は-1
  Idx parent(Idx a){
    return a ? (a - 1) / k : -1;
  }
  // aのt個親, 深さを超える場合は-1
  Idx la(Idx a, int t){
    for(int i = 0; i < t; i++){
      if(!a) return -1;
      a = (a - 1) / k;
    }
    return a;
  }
  // lca
  Idx lca(Idx a, Idx b){
    return dist2(a, b).first;
  }
  // {ノードの深さ, その部分木の深さh-1のノードがいくつ欠けているか}でノードを分類すると, その種類数は高々3h
  // {個数, ノードの深さ, 深さh-1のノードがいくつ欠けているか}を返す
  std::vector<std::tuple<Idx, Idx, Idx>> depth_frequency_decompose(){
    std::vector<std::tuple<Idx, Idx, Idx>> ret;
    Idx x = N - 1, nd = 1;
    int h = (int)Lelem.size();
    for(int d = h - 1; d >= 0; d--){
      Idx L = x - Lelem[d], R = kpow[d] - 1 - L;
      if(L) ret.push_back({L, d, 0});
      if(R) ret.push_back({R, d, kpow[h - 1 - d]});
      ret.push_back({1, d, kpow[h - 1 - d] - nd});
      if(d){
        x--;
        nd += kpow[h - 1 - d] * (x % k);
        x /= k;
      }
    }
    return ret;
  }
  // aの部分木の頻度テーブル(ans[i] := aの部分木に含まれaとの距離がiのノード数)
  std::vector<Idx> depth_frequency(Idx a){
    if(a >= N) return {};
    std::vector<Idx> ret;
    int t = 0;
    while(a < M){
      if(N <= a){
        ret.push_back(std::max(Idx(0), kpow[t++] - (a - N + 1)));
        return ret;
      }else{
        ret.push_back(kpow[t++]);
      }
      a = a * k + k;
    }
    return ret;
  }
  // aの部分木の頂点でaとの距離がdの頂点の数
  Idx count_dist_subtree(Idx a, int d){
    if(d < 0 || a >= N) return 0;
    int da = depth(a);
    return __count_dist_subtree(a, da, d);
  }
  Idx __count_dist_subtree(Idx a, int da, int d){
    if(d < 0 || a >= N) return 0;
    int h = Lelem.size();
    if(da + d >= h) return 0;
    if(da + d < h - 1){
      return kpow[d];
    }else{
      Idx ldeep = (a - Lelem[da]) * kpow[d];
      return std::min(kpow[d], std::max(Idx(0), num_deepest() - ldeep));
    }
  }
  // aとの距離がdの頂点の数
  Idx count_dist(Idx a, int d){
    if(d < 0 || a >= N) return 0;
    Idx ans = 0;
    int da = depth(a);
    while(d >= 0){
      ans += __count_dist_subtree(a, da, d);
      if(!a) return ans;
      ans -= __count_dist_subtree(a, da, d - 2);
      d--, da--;
      a = (a - 1) / k;
    }
    return ans;
  }
};

#endif