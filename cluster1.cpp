#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>

#include "glpk.h"

using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<
const unsigned int MAX_COEFF_NR = 500000;
const float FDR = 0.05;
const int C = 1000;

int main() {
  ifstream is("constraints2.in");
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
  int ia[1+MAX_COEFF_NR], ja[1+MAX_COEFF_NR];
  double ar[1+MAX_COEFF_NR];
  glp_add_rows(mip, N);
  long indx{1};
  for (int i=1;i<=N;i++) {
    glp_set_row_bnds(mip, i, GLP_DB, 0.0, 1.0);
    int cur_node;
    is >> cur_node;
    ia[indx]=i;ja[indx]=cur_node;ar[indx]=1;indx++;
    //cout << cur_node << ": ";
    char ch;
    is.get(ch);//:
    string s;
    getline(is,s);
    istringstream iss(s);
    while (true) {
      iss >> cur_node;
      if (!iss) break;
      ia[indx]=i;ja[indx]=cur_node;ar[indx]=1;indx++;
      //cout << cur_node << " ";
    }
    //cout << endl;
  }
  int row_n = N + 1;
  while (true) {
    int cur_node1, cur_node2;
    float qij;
    is >> cur_node1;
    if (!is) break;
    is >> cur_node2 >> qij;
    glp_add_rows(mip, 1);//add constraint row
    glp_set_row_bnds(mip, row_n, GLP_UP, 0.0, 2*C + FDR - qij);
    ia[indx]=row_n;ja[indx]=cur_node1;ar[indx]=C;indx++;
    ia[indx]=row_n;ja[indx]=cur_node2;ar[indx]=C;indx++;
    row_n++;
    //cout << cur_node1 sp cur_node2 sp qij << endl;
  }
  is.close();
  glp_load_matrix(mip, indx-1, ia, ja, ar);
  glp_iocp parm;
  glp_init_iocp(&parm);
  parm.presolve = GLP_ON;
  int err = glp_intopt(mip, &parm);
  cout << err << endl;
  cout << "Object value: " << glp_mip_obj_val(mip) << endl;
  for (int i=1;i<=N;i++) {
    if (glp_mip_col_val(mip, i) == 1 ) cout << glp_get_col_name(mip, i) << endl;
  }
  glp_delete_prob(mip);
}
