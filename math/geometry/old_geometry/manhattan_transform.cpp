template<typename T>
pair<T, T> manhattan_transform(pair<T, T> v){
  return {v.first+v.second, v.first-v.second};
}
template<typename T>
vector<T> manhattan_transform(vector<T> v){
  int k = v.size();
  vector<T> res(1<<(k-1), v[0]);
  for(int i=0;i<(1<<(k-1));i++){
    for(int j=0;j<(k-1);j++){
      if(!bit_kth(i, j)) res[i] += v[j+1];
      else res[i] -= v[j+1];
    }
  }
  return res;
}
