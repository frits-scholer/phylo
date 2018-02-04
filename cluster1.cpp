#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

#include "glpk.h"

using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<

int main() {
  ifstream is("constraints1.in");
  if (!is) {
    cerr << "Could not open file\n";
    return 1;
  }


  int N;
  is >> N;
  glp_prob *mip = glp_create_prob();
  glp_set_prob_name(mip, "cluster");
  glp_set_obj_dir(mip, GLP_MAX);
  
  glp_add_cols(mip, N);
  //set object
  for (int i=1;i<=N;i++) {
    string s;
    int l;
    is >> s >> l;
    glp_set_col_name(mip, i, s.c_str());
    glp_set_col_kind(mip, i, GLP_BV);
    glp_set_obj_coef(mip, i, l);
  }
  /*
  vector<int> x(N);
  vector<int> ia(N+1), ja(N+1);
  vector<double> ar[N+1];
  double Z;
  */
  is.close();

  glp_delete_prob(mip);
}
