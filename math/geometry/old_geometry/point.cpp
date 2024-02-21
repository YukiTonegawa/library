/*
complexではなくpointクラスを作った.

norm: ノルム
abs: 原点からの距離
arg: 角度(−π
−
π
からπ
π
)
dot: 内積
cross,det: 外積(行列式)
proj: 直線(0,0),(x,y)への正射影
ccwはa->b->cの順で進むとき
1ならccw(反時計回り), -1ならcw(時計回り), 2なら折り返し, -2なら直進, 0なら線分上に折り返し
*/


namespace D2 {
const ld eps = 1e-7;
bool eq(const ld &a, const ld &b){ return std::abs(a-b)<eps;}
struct point {
    ld x, y;
    point():x(0),y(0) {}
    point(ld x, ld y):x(x),y(y) {}

    point operator+=(const point &p) { x+=p.x; y+=p.y; return *this; }
    point operator-=(const point &p) { x-=p.x; y-=p.y; return *this; }
    point operator*=(const point &p) { ld tx = x; x=x*p.x-y*p.y; y=y*p.x+tx*p.y; return *this; }
    point operator*=(const ld &d) { x*=d; y*=d; return *this; }
    point operator/=(const ld &d) { x/=d; y/=d; return *this; }
    point operator+(const point &p) const { return point(*this)+=p; }
    point operator-(const point &p) const { return point(*this)-=p; }
    point operator*(const ld &d) const { return point(*this)*=d; }
    point operator*(const point &p) const { return point(*this)*=p; }
    point operator/(const ld &d) const { return point(*this)/=d; }
    point operator-() const { return point(-x,-y); }
    bool operator<(const point &p) const { return std::abs(x- p.x)>eps ? x < p.x-eps : y < p.y-eps; }
    bool operator==(const point &p) const { return std::abs(x-p.x)<eps and std::abs(y-p.y)<eps; }
    bool operator!=(const point &p) const { return !(point(*this)==p); }

    ld norm() const { return x*x+y*y; }
    ld abs() const { return sqrt(norm()); }
    ld arg() const { return atan2(y, x); }

    ld dot(const point &p) const { return x*p.x+y*p.y; }
    ld cross(const point &p) const { return x*p.y-y*p.x; }
    ld det(const point &p) const { return x*p.y-y*p.x; }
    point proj(const point &p) const { ld k = dot(p)/norm(); return point(x*k,y*k); }
    ld pi(){
      ld a = arg();
      if(a < 0) return (2*M_PI + a);
      return a;
    }
};
point operator*(const ld &d, const point &p) { point q = p*d; return q; }
std::ostream &operator<<(std::ostream &os, const point &p) { os << "(" << p.x << "," << p.y << ")"; return os; }

point polar(const ld &r, const ld &th) { return point(r*std::cos(th),r*std::sin(th)); }
ld abs(const point &p) { return p.abs(); }
ld norm(const point &p) { return p.norm(); }
ld dot(const point &a, const point &b) { return a.dot(b); }
ld cross(const point &a, const point &b) { return a.cross(b); }
ld det(const point &a, const point &b) { return a.det(b); }
ld arg(const point &a) { return a.arg(); }
/*
 *           1
 *  ------------------
 *    -2 |a  0  b| 2
 *  ------------------
 *          -1
 */
int ccw(point a, point b, point c)
{
    b -= a, c -= a;
    if (b.cross(c) > eps)        return +1;
    if (b.cross(c) < -eps)       return -1;
    if (b.dot(c) < -eps)         return -2; // c -- a -- b
    if (b.norm() < c.norm()-eps) return +2; // a -- b -- c
    return 0;
}
}
