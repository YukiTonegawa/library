/*
直線クラス. 線分も兼ねている.

vec: 方向ベクトル
norm: ノルム
abs: 原点からの距離
proj: 正射影
refl: 線対称な点
各関数名は以下の規則に従っている.

接頭辞(1文字, 動作)
i: 交差判定
c: 交点
d: 距離
接尾辞(2文字, 引数)
l: 直線
s: 線分
p: 点
*/

namespace D2 {
/*
  depends on
  point
*/
struct line {
    point a, b;
    line(){}
    line(point a, point b) : a(a), b(b) {}

    point vec() const { return b-a; }
    ld abs() const { return vec().abs(); }
    ld norm() const { return vec().norm(); }
    point proj(const point &p) const { return a+vec().proj(p-a); }
    point refl(const point &p) const { return proj(p)*2-p; }
};
std::ostream &operator<<(std::ostream &os, const line &l) { os << l.a << " - " << l.b; return os; }

bool ill(const line &l, const line &m) { return std::abs(l.vec().det(m.vec())) > eps; }
bool ils(const line &l, const line &s) { return ccw(l.a,l.b,s.a)*ccw(l.a,l.b,s.b)<=0; }
bool isl(const line &s, const line &l) { return ils(l,s); }
bool iss(const line &s, const line &t) { return ils(s,t) and ils(t,s); }
point cll(const line &l, const line &m) { return l.a+l.vec()*(m.vec().det(m.a-l.a)/m.vec().det(l.vec())); }
ld dlp(const line &l, const point &p) { return (l.proj(p)-p).abs(); }
ld dpl(const point &p, const line &l) { return dlp(l,p); }
ld dll(const line &l, const line &m) { return ill(l,m) ? 0.0 : dlp(l,m.a); }
ld dps(const point &p, const line &s) { return ccw(s.a,s.b,s.proj(p)) ? std::min((s.a-p).abs(), (s.b-p).abs()) : (s.proj(p)-p).abs(); }
ld dsp(const line &s, const point &p) { return dps(p,s); }
ld dss(const line &s, const line &m) { return iss(s,m)? 0.0 : std::min(std::min(dps(m.a,s),dps(m.b,s)), std::min(dps(s.a,m),dps(s.b,m))); }

}
