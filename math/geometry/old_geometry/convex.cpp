/*
求まる凸包は反時計回りの順になっている.
*/

namespace D2 {
std::vector<point> convex_hull(std::vector<point> ps) {
    int n = ps.size(), k = 0;
    if(n < 3) return ps;
    sort(begin(ps),end(ps));
    std::vector<point> ret(2*n);
    for (int i = 0; i < n; ret[k++] = ps[i++])
        while (k >= 2 and ccw(ret[k-2], ret[k-1], ps[i]) < 1) k--;
    for (int i = n-2, t = k+1; i >= 0; ret[k++] = ps[i--])
        while (k >= t and ccw(ret[k-2], ret[k-1], ps[i]) < 1) k--;
    ret.resize(k-1);
    return ret;
}


ld max_x=-1e9, max_y = -1e9, min_x = 1e9, min_y = 1e9;
// 点Pが凸包の境界線上および内部にあるか
int convex_contains(vector<point> &P, point p) {
  if(P.size()==0||P.size()==1) return 0;
  int n = P.size();
  point g = (P[0] + P[n/3] + P[2*n/3]) / 3.0; // inner-point
  int a = 0, b = n;

  //凸包が線の場合
  point x = P[1] - P[0];
  point y = P[n-1] - P[0];
  if(abs(x.arg()-y.arg()) < eps){
    point pl = P[0] - p;
    point pr = p - P[0];
    if((abs(pl.arg()-x.arg())<eps || abs(pr.arg()-x.arg()) < eps)
      &&((max_x-p.x)*(p.x-min_x)>0||(max_y-p.y)*(p.y-min_y)>0)) {
        return 1;
    }
    else return 0;
  }
  while (b-a>1) { // invariant: c is in fan g-P[a]-P[b]
    int c = (a + b) / 2;
    if (cross(P[a]-g, P[c]-g) > 0){ // angle < 180 deg
      if (cross(P[a]-g, p-g)>0&&cross(P[c]-g, p-g)<0) {
        b = c;
      }else{
        a = c;
      }
    }else{
      if (cross(P[a]-g, p-g) < 0 && cross(P[c]-g, p-g) > 0){
        a = c;
      }else{
        b = c;
      }
    }
  }
  b %= n;
  if (cross(P[a]-p, P[b]-p)<0) return 0;
  if (cross(P[a]-p, P[b]-p)>0) return 2;
  return 1;
}
}
