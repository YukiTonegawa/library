using namespace D2;
ld theta(point a, point b){
  return acos(dot(a, b)/sqrtl(dot(a, a)*dot(b, b)));
}
//3点から半径取得
ld get_radius(point a, point b, point c){
  return sqrtl(dot((b-c), (b-c)))/sinl(theta(b-a,c-a))/2;
}
bool equal(ld a, ld b){
  if(a < b + 1e-4 && a + 1e-4 > b) return true;
  return false;
}
bool equal(point a, point b){
  if(equal(a.x, b.x)&&equal(a.y, b.y)) return true;
  return false;
}
