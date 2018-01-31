#include <iostream>
#include <iterator>
#include <vector>
#include <numeric>
#include <algorithm>
using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<
typedef pair<float,int> p_ind;
typedef vector<p_ind> sorted_p;

vector<float> ecdf(unsigned int p_size) {
  vector<float> v(p_size);
  iota(all(v),1);
  for_each(all(v),[=](float& x){x /= p_size;});
  return v;
}
bool greater1(float x) {return x>1;}

void bh_fdr(const vector<float>& p, vector<float>& q) {
  int T = p.size();
  sorted_p pvals(T);
  for (int i = 0; i < T;i++) {
    pvals[i].first = p[i];
    pvals[i].second=i;
  }
  sort(all(pvals));
  vector<float> ecdffactor;
  ecdffactor = ecdf(T);
  vector<float> pvals_corrected_raw(T);
  transform(all(pvals), begin(ecdffactor), begin(pvals_corrected_raw),
  	    [](p_ind pv, float x){return pv.first/x;});
  //pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
  //first reverse pvals_corrected_raw
  reverse(all(pvals_corrected_raw));
  vector<float> pvals_corrected(T);
  pvals_corrected[0]=pvals_corrected_raw[0];
  for (int i=1;i<T;i++) {
    pvals_corrected[i]=min(pvals_corrected_raw[i], pvals_corrected[i-1]);
  }
  //reverse pvals_corrected
  reverse(all(pvals_corrected));
  //replace all values greater than 1 by 1
  replace_if(all(pvals_corrected),greater1,1);
  vector<pair<int,float>> qvals(T);
  for (int i=0;i<T;i++) {
    qvals[i] = make_pair(pvals[i].second, pvals_corrected[i]);
  }
  sort(all(qvals));
  for (int i=0;i<T;i++) {
    q[i] = qvals[i].second;
  }
}
