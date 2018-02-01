#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

#define all(t) begin(t), end(t)
#define EPS1 0.001
#define EPS2 1.0e-8

float probks(float alam) {
//Kolmogorov-Smirnov probability function.
  float fac=2.0,sum=0.0,termbf=0.0;
  float a2 = -2.0*alam*alam;
  for (int j=1;j<=100;j++) {
    float term=fac*exp(a2*j*j);
    sum += term;
    if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
    fac = -fac;//Alternating signs in sum.
    termbf=fabs(term);
  }
  return 1.0;
  //Get here only by failing to converge.
}

pair<float, float> kstwo(const vector<float>& data1, const vector<float>& data2) {
/*
Given an array data1[1..n1] , and an array data2[1..n2] , 
this routine returns the K–
S statistic d , and the significance level prob for the null hypothesis 
that the data sets are drawn from the same distribution.
Small values of prob show that the cumulative distribution
function of data1 is significantly different from that of data2 .
The arrays data1 and data2 are modified by being sorted into ascending order.
*/
  //sort(all(data1));
  //presorted to avoid repeated sorting
  //sort(all(data2));
  //presorted
  float en1=data1.size();
  float en2=data2.size();
  float d=0.0;
  unsigned long j1=0,j2=0;
  float d1,d2,dt,en,fn1=0.0,fn2=0.0;
  //data1 en data2 are 0-indexed
  while (j1 < data1.size() && j2 < data2.size()) {
    d1=data1[j1];
    d2=data2[j2];
    if (d1 <= d2) {j1++;fn1 = j1/en1;}
    //if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=j1++/en1;
    if (d2 <= d1) {j2++;fn2 = j2/en2;}
    dt=fabs(fn2-fn1);
    if (dt > d) d = dt;//the max
  }
  en=sqrt(en1*en2/(en1+en2));
  float prob=probks((en+0.12+0.11/en)*d);
  return make_pair(d, prob);
}
