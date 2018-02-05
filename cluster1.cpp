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
  int ia[1+N*N], ja[1+N*N];
  double ar[1+N*N];
  glp_add_rows(mip, N);
  long indx{1};
  for (int i=1;i<=N;i++) {
    glp_set_row_bnds(mip, i, GLP_DB, 0.0, 1.0);
    int cur_node;
    is >> cur_node;
    ia[indx]=i;ja[indx]=cur_node;ar[indx]=1;indx++;
    cout << cur_node << ": ";
    char ch;
    is.get(ch);//:
    string s;
    getline(is,s);
    istringstream iss(s);
    while (true) {
      iss >> cur_node;
      if (!iss) break;
      ia[indx]=i;ja[indx]=cur_node;ar[indx]=1;indx++;
      cout << cur_node << " ";
    }
    cout << endl;
  }
  is.close();
  glp_load_matrix(mip, indx-1, ia, ja, ar);
  glp_iocp parm;
  glp_init_iocp(&parm);
  parm.presolve = GLP_ON;
  int err = glp_intopt(mip, &parm);
  cout << "Object value: " << glp_mip_obj_val(mip) << endl;
  for (int i=1;i<=N;i++) {
    cout << glp_get_col_name(mip, i) << " = " << glp_mip_col_val(mip, i) << endl;
  }
  glp_delete_prob(mip);
}
