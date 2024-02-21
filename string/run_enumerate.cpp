/*
// Reference:
// D. Gusfield,
// Algorithms on Strings, Trees, and Sequences: Computer Science and
// Computational Biology
template <class T> vector<int> z_algorithm(const vector<T>& s) {
  int n = int(s.size());
  if (n == 0) return {};
  vector<int> z(n);
  z[0] = 0;
  for (int i = 1, j = 0; i < n; i++) {
    int& k = z[i];
    k = (j + z[j] <= i) ? 0 : min(j + z[j] - i, z[i - j]);
    while (i + k < n && s[k] == s[i + k]) k++;
    if (j + z[j] < i + z[i]) j = i;
  }
  z[0] = n;
  return z;
}

vector<int> z_algorithm(const string& s) {
  int n = int(s.size());
  vector<int> s2(n);
  for (int i = 0; i < n; i++) {
    s2[i] = s[i];
  }
  return z_algorithm(s2);
}

typedef array<int, 3> ary;
void Main_Lorentz(const string& s, const string& rev, vector<ary>& ret, int l=0, int r=-1){
  int n = (int)s.size();
  if(r==-1) r = (int)s.size();
  if(r-l<2) return;
  int mid = (l + r)/2;

  string L = rev.substr(n-mid, mid - l);
  string R = s.substr(mid, r - mid);
  vector<int> z1 = z_algorithm(L);
  vector<int> z2 = z_algorithm(R + "#" + s.substr(l, r - l));
  int len = r - l;
  for(int i=0;i<mid-l;i++){
    int x = i;
    int y = z1[i];
    if(i==0) x = mid - l, y = 0;
    int z = z2[len+1-x];
    int lidx = mid - x - y, ridx = mid + z;
    if(x>0&&z>0&&(y+z)>=x) ret.push_back({x, lidx, ridx});//周期, [l, r)
  }
  z1 = z_algorithm(R);
  z2 = z_algorithm(L + "#" + rev.substr(n-r, r - l));
  for(int i=0;i<r-mid;i++){
    int x = i;
    int y = z1[i];
    if(i==0) x = r-mid, y = 0;
    int z = z2[len+1-x];
    int ridx = n - (n - mid - x - y), lidx = n - (n - mid + z);
    if(x>0&&z>0&&(y+z)>=x) ret.push_back({x, lidx, ridx});//周期, [l, r)
  }
  Main_Lorentz(s, rev, ret, l, mid);
  Main_Lorentz(s, rev, ret, mid, r);
}

vector<ary> Main_Lorentz(const string& s){
  int n = (int)s.size();
  //(t, l, r)は周期tを持つ
  //同じ(l, r)での周期を最小にする
  //同じ(l, 周期)でのrを最大化する
  //同じ(r, 周期)でのlを最小化する
  vector<ary> tmp, ret;
  string t = s;
  reverse(all(t));
  Main_Lorentz(s, t, tmp);
  sort(all(tmp), [&](const ary& x, const ary& y){
    if(x[0]!=y[0]) return x[0] < y[0];
    if(x[1]!=y[1]) return x[1] < y[1];
                   return x[2] > y[2];
  });
  int l = 1000000000, r = -1000000000;
  set<array<int, 2>> st;
  for(int i=0;i<(int)tmp.size();i++){
    ary v = tmp[i];
    if(i&&v[0]!=tmp[i-1][0]) l = 1000000000, r = -1000000000;
    if(v[1]>=l&&v[2]<=r) continue;
    if(n<v[2]) continue;
    l = min(l, v[1]);
    r = max(r, v[2]);
    if(st.find({v[1], v[2]})!=st.end()) continue;
    ret.push_back(v);
    st.insert({v[1], v[2]});
  }
  return ret;
}


template<typename T, T INF>
void Main_Lorentz(const vector<T>& s, const vector<T>& rev, vector<ary>& ret, int l=0, int r=-1){
  int n = (int)s.size();
  if(r==-1) r = (int)s.size();
  if(r-l<2) return;
  int mid = (l + r)/2;

  vector<T> L(mid-l), R(r-mid), R2(r-mid+1+r-l), L2(mid-l+1+r-l);
  for(int i=0;i<mid-l;i++) L[i] = L2[i] = rev[n-mid+i];
  for(int i=0;i<r-mid;i++) R[i] = R2[i] = s[mid+i];
  L2[mid-l] = R2[r-mid] = INF;//insert #
  for(int i=0;i<r-l;i++){
    R2[r-mid+1+i] = s[l+i];
    L2[mid-l+1+i] = rev[n-r+i];
  }

  vector<int> z1 = z_algorithm<T>(L);
  vector<int> z2 = z_algorithm<T>(R2);
  int len = r - l;
  for(int i=0;i<mid-l;i++){
    int x = i;
    int y = z1[i];
    if(i==0) x = mid - l, y = 0;
    int z = z2[len+1-x];
    int lidx = mid - x - y, ridx = mid + z;
    if(x>0&&z>0&&(y+z)>=x) ret.push_back({x, lidx, ridx});//周期, [l, r)
  }
  z1 = z_algorithm<T>(R);
  z2 = z_algorithm<T>(L2);

  for(int i=0;i<r-mid;i++){
    int x = i;
    int y = z1[i];
    if(i==0) x = r-mid, y = 0;
    int z = z2[len+1-x];
    int ridx = n - (n - mid - x - y), lidx = n - (n - mid + z);
    if(x>0&&z>0&&(y+z)>=x) ret.push_back({x, lidx, ridx});//周期, [l, r)
  }
  Main_Lorentz<T, INF>(s, rev, ret, l, mid);
  Main_Lorentz<T, INF>(s, rev, ret, mid, r);
}

template<typename T, T INF>
vector<ary> Main_Lorentz(const vector<T>& s){
  int n = (int)s.size();
  //(t, l, r)は周期tを持つ
  //同じ(l, r)での周期を最小にする
  //同じ(l, 周期)でのrを最大化する
  //同じ(r, 周期)でのlを最小化する
  vector<ary> tmp, ret;
  vector<T> t = s;
  reverse(all(t));
  Main_Lorentz<T, INF>(s, t, tmp);
  sort(all(tmp), [&](const ary& x, const ary& y){
    if(x[0]!=y[0]) return x[0] < y[0];
    if(x[1]!=y[1]) return x[1] < y[1];
                   return x[2] > y[2];
  });
  int l = 1000000000, r = -1000000000;
  set<array<int, 2>> st;
  for(int i=0;i<(int)tmp.size();i++){
    ary v = tmp[i];
    if(i&&v[0]!=tmp[i-1][0]) l = 1000000000, r = -1000000000;
    if(v[1]>=l&&v[2]<=r) continue;
    if(n<v[2]) continue;
    l = min(l, v[1]);
    r = max(r, v[2]);
    if(st.find({v[1], v[2]})!=st.end()) continue;
    ret.push_back(v);
    st.insert({v[1], v[2]});
  }
  return ret;
}
*/
