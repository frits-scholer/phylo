#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <algorithm>
#include <string>
#include <cctype>
#include <map>
#include <utility>
#include <iterator>
#include <list>
#include "kstwo.hpp"
#include "bh_fdr.hpp"
#include "glpk.h"
using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<
#define tb << "\t" <<

const float epsilon = 0.00001;
const unsigned int MAX_COEFF_NR = 10000000;
const float C = 1000;

enum Token_value {
  NAME,END,
  LP='(',RP=')',NUMBER=':',COMMA=','
};

Token_value curr_tok=END;
string string_value;
float number_value;
int ml;

map<string,int> node_indx;

struct node;
typedef vector<float> distribution;
typedef list<node*> nodelist;
typedef vector<node*> nodevector;
typedef pair<float, node*> node_mean;

struct node {
  node *parent;
  nodelist children;
  int nrleaves;
  float distance;
  string info;
  bool selected;
  bool isleaf;
  distribution D;//needed for KS
  node(): parent(nullptr), nrleaves(0), selected(false), isleaf(false) {}
};

vector<pair<node*,node*>> zero_nodes;
nodevector merge_nodes;

bool ladder_sort(const node *a, node *b) {
  return (a->nrleaves < b->nrleaves) ||
    (a->nrleaves == b->nrleaves && a->distance < b->distance);
}

Token_value get_token(istream& is) {
  char ch;
  do {
    if (!is.get(ch)) return curr_tok = END;
  }  while (ch !=';' && isspace(ch));
  switch(ch) {
  case '(':case ')':case ',':
    return curr_tok = Token_value(ch);
  case ';':
    return curr_tok = END;
  case ':':    
    is >> number_value;
    return curr_tok = NUMBER;
  case '\'':
    string_value.clear();
    while (is.get(ch)) {if (ch == '\'') break;string_value.push_back(ch);}
    return curr_tok = NAME;
  default:
    is.putback(ch);
    string_value.clear();
    if (curr_tok == LP || curr_tok == COMMA) {
      getline(is, string_value,':');
      is.putback(':');
      return curr_tok = NAME;}
    cerr << "ERROR\n";
    return curr_tok = END;
  }
}

node* build_tree(vector<node*>& leaves) {
  cerr << "filename?\n";
  string fname;
  cin >> fname;
  ifstream is(fname.c_str());
  if (!is) {
    cerr << "Could not open file\n";
    return nullptr;
  }
  node *root = new node;
  node *cur_node = root;
  root->parent = root;
  int internal_count = 0;
  stack<node *> A;
  node *nptr;
  while (is) {
    get_token(is);
    if (curr_tok == END) break;
    switch(curr_tok) {
    case LP:
      A.push(cur_node);
      cur_node->info = to_string(internal_count);
      internal_count++;
      nptr = new node;
      nptr->parent = cur_node;
      cur_node->children.push_back(nptr);
      cur_node = nptr;
      break;
    case RP:
      cur_node = A.top();
      A.pop();
      break;
    case NAME:
      cur_node->info = string_value;
      cur_node->isleaf = true;
      leaves.push_back(cur_node);
      break;
    case NUMBER:
      cur_node->distance = number_value;
      break;
    case COMMA:
      cur_node = A.top();
      nptr = new node;
      nptr->parent = cur_node;
      cur_node->children.push_back(nptr);
      cur_node = nptr;
      break;
    case END:
      cout << "error";
    }
  }
  return root;
}
  
int nr_of_children(node *root) {
  return root->children.size();
}

bool is_root(node *root) {
  return root == root->parent;
}

void rzb_nodes(nodelist& lst) {
  for (auto nptr : lst) {
    if (nptr->isleaf) continue;
    rzb_nodes(nptr->children);
    if (nptr->distance > epsilon) continue;
    node *from = nptr;
    node *to = nptr->parent;
    while (!is_root(to) && to->distance <= epsilon) to = to->parent;
    zero_nodes.push_back(make_pair(from, to));
  }
}
void rzn(nodelist& lst) {
  for (auto nptr : lst) {
    if (nptr->isleaf) continue;
    rzn(nptr->children);
    int nch = nr_of_children(nptr);
    if (nch != 1) continue;
    merge_nodes.push_back(nptr);
  }
}

void append_stream(ostream& os, node* root) {
  if (root->isleaf) {os << root->info << ':' << root->distance;return;}
  os << '(';
  for (auto it = begin(root->children);it != end(root->children);it++) {
    if (it != begin(root->children)) os << ',';
    append_stream(os, *it);
  }
  os << ')';
  if (is_root(root)) os << ";";
  else os << ':' << root->distance;
}

void write_newick(node *root) {
  cerr << "filename of output newick file?\n";
  string fname;
  cin >> fname;
  ofstream os(fname.c_str());
  if (!os) {
    cerr << "Could not open file\n";
    return;
  }
  append_stream(os, root);
  os.close();
}

int nrLeaves(node *root) {
  if (root->isleaf) {root->nrleaves = 1;return 1;}
  int N=0;
  for (auto it = begin(root->children);it != end(root->children);it++) {
    N +=  nrLeaves(*it);
  }
  root->nrleaves = N;
  return N;
}

void select_clades(node *root) {
  if (root->isleaf) return;
  for (auto it = begin(root->children);it != end(root->children);it++) {
    select_clades(*it);
  }    
  if (root->nrleaves >= ml ) root->selected = true;
  else root->selected = false;
  return;
}

float ancestor_distance(node* z, node* w) {//w is descendant of z
  float dist = 0;
  while (w != z) {
    dist += w->distance;
    w = w->parent;
  }
  return dist;
}

node* common_ancestor(node* x, node* y) {
  node* ax = x->parent;
  while (true) {
    node* ay = y->parent;
    while (true) {
      if (ax == ay) return ax;
      if (is_root(ay)) break;
      ay = ay->parent;
    }
    ax = ax->parent;
  }
}

void process_distance(node* x, node* y) {
  node* z = common_ancestor(x, y);
  float ild = ancestor_distance(z, x) + ancestor_distance(z, y);
  while (true) {
    z->D.push_back(ild);//generate D even if internal node is not selected
    if (is_root(z)) return;
    z = z->parent;
  }
}

float calc_distance(node* x, node* y) {
  node* z = common_ancestor(x, y);
  return  ancestor_distance(z, x) + ancestor_distance(z, y);
}

void calc_mean(node *root, vector<node_mean>& v) {
  for (auto it = begin(root->children);it != end(root->children);it++) {
    calc_mean(*it, v);
  }
  if (!(root->isleaf) && root->selected) {
    float S = accumulate(all(root->D),0.0);
    v.push_back(node_mean(S/(root->D).size(), root));
    //auto a = v.back();
    //cerr << a.first sp a.second->info << endl;
  }
}

void printLeaves(node *root) {
  for (auto it = begin(root->children);it != end(root->children);it++) {
    printLeaves(*it);
  }
  if (root->isleaf) cout << root->info << '\n';
}

void printAncestors(node *root) {//prints all nontrivial ancestors
  node* z = root;
  while (true) {
    z = z->parent;
    if (is_root(z)) break;
    cout << z->info << " ";
  }
}

void printNodes(node *root) {
  for (auto it = begin(root->children);it != end(root->children);it++) {
    printNodes(*it);
  }
  if (!root->isleaf) {
    cout << root->info sp root->children.size() sp ':';
    printAncestors(root);
    cout << endl;
  }
}

void ladderize(node *root) {
  for (auto it = begin(root->children);it != end(root->children);it++) {
    ladderize(*it);
  }
  root->children.sort(ladder_sort);
}

bool is_ancestor(node* x, node* y) {//x is ancestor of y
  node* z = y->parent;
  do {
    if (x == z) return true;
    z = z->parent;
  } while (!is_root(z));
  return x==z;
} 

void preSort(node *root) {
  for (auto it = begin(root->children);it != end(root->children);it++) {
    preSort(*it);
  }
  if (!(root->isleaf)) sort(all(root->D));
}

void IndexAncestors(node *root, int *ia, int *ja, double *ar, int i, long& indx) {
  node* z = root;
  do {
    z = z->parent;
    if (node_indx[z->info]) {
      ia[indx]=i;
      ja[indx] = node_indx[z->info];
      ar[indx]=1;
      indx++;
    }
  } while (!is_root(z));
}

void show_event(string s, clock_t& tm) {
  tm = clock()-tm;
  cerr <<  "\t" << s << " " << (double) tm/CLOCKS_PER_SEC << " s "<< endl;
}

void cleanup_zero_nodes() {
  for (auto zn : zero_nodes) {
    auto it = begin(zn.second->children);//position to move to
    for (auto nptr : zn.first->children) nptr->parent = zn.second;
    zn.second->children.splice(it, zn.first->children);
    zn.first->parent->children.remove(zn.first);
  }
}

void cleanup_merge_nodes() {
  for (auto mn : merge_nodes) {
    auto it = begin(mn->parent->children);//position to move to
    for (auto nptr : mn->children) {
      nptr->parent = mn->parent;
      nptr->distance += mn->distance;
    }
    mn->parent->children.splice(it, mn->children);
    mn->parent->children.remove(mn);
  }
}

int main() {
  nodevector leaves;
  node *root = build_tree(leaves);
  if (!root) return 1;
  cerr << "Do you want to remove zero length branches?(y/n)\n";
  char zl;
  cin >> zl;
  float gamma;
  cerr << "Minimum nr of leaves of considered clades?\n";
  cin >> ml;
  cerr << "Gamma factor?\n";
  cin >> gamma;
  float FDR;
  cerr << "FDR (False Discovery Rate)?\n";
  cin >> FDR;
  //end of input
  clock_t tm=clock();
  if (zl == 'y') {
    for (auto chld : root->children) {
      rzb_nodes(chld->children);
    }
    cleanup_zero_nodes();
    rzn(root->children);
    cleanup_merge_nodes();
  }
  nrLeaves(root);
  ladderize(root);
  write_newick(root);//give output nwk filename
  cout << "Minimum nr of leaves:" sp ml sp "Gamma: " sp gamma sp "FDR: " sp FDR << endl;
  //start timer
  //printNodes(root);
  select_clades(root);
  root->selected = true;
  //set up interleaf distances
  for (auto il=begin(leaves)+1;il != end(leaves);il++) {
    for (auto jl = begin(leaves);jl != il;jl++) {
      process_distance(*il,*jl);
    }
  }
  vector<node_mean> means;
  calc_mean(root, means);
  sort(all(means));
  //for (auto a : means) cerr << a.first sp a.second->info << endl;
  size_t msize = means.size();
  float gm;
  if (msize & 1) gm = means[msize/2].first;
  else gm = (means[msize/2-1].first + means[msize/2].first)/2.0;
  cout << "size of means = " sp msize sp "grand median = " << gm << endl;
  auto it = lower_bound(all(means),node_mean(gm,nullptr));
  msize = distance(it, end(means));
  float mad;
  if (msize & 1) mad = (it + msize/2)->first - gm;
  else mad = ((it + msize/2-1)->first + (it + msize/2)->first)/2.0 - gm;
  cout << "mad = " << mad << endl;
  float upperbound = gm + gamma * mad;
  it = upper_bound(all(means), node_mean(upperbound+0.00000001, nullptr));
  means.erase(it, end(means));
  cout << "size of filtered means = " << means.size() << endl;
  //presort all data
  preSort(root);
  nodevector sel_nodes;
  for_each(means.rbegin(),means.rend(),[&](node_mean m){sel_nodes.push_back(m.second);});
  vector<bool> sel;
  vector<float> p;
  int N = sel_nodes.size();  
  for (int i = 1;i < N;i++) {
    for (int j = 0;j < i;j++) {
      sel.push_back(false);
      //cout << sel_nodes[i]->info sp sel_nodes[j]->info << '\t';
      if (is_ancestor(sel_nodes[i], sel_nodes[j]) ||
	  is_ancestor(sel_nodes[j], sel_nodes[i])) {
	auto ks = kstwo(sel_nodes[i]->D, sel_nodes[j]->D);
	p.push_back(ks.second);
	sel.back()=true;
	}
      else {
	node* ca = sel_nodes[i]->parent;
	if (ca != sel_nodes[j]->parent) continue; 
	auto ksi = kstwo(sel_nodes[i]->D, ca->D);
	auto ksj = kstwo(sel_nodes[j]->D, ca->D);
	p.push_back(max(ksi.second, ksj.second));
	sel.back()=true;
      }
    }
  }
  //cerr << p.size() << endl;
  vector<float> q(p.size());
  bh_fdr(p,q);
  //auxiliary output
  /*
  for (auto n: sel_nodes) {
    cout << n->info << ": ";
    printAncestors(n);
    cout << endl;
  }
  for (auto n: sel_nodes) {
    cout << n->info << ": ";
    printLeaves(n);
    cout << endl;
  }
  ostream_iterator<float> os(cout, " ");
  for (auto n: sel_nodes) {
    cout << n->info << ": ";
    copy(all(n->D),os);
    cout << endl;
  }
  for (unsigned int i = 1;i < N;i++) {
    for (unsigned int j = 0;j < i;j++) {
      cout << sel_nodes[i]->info sp sel_nodes[j]->info << '\t'
	   << calc_distance(sel_nodes[i], sel_nodes[j]) << endl;
    }
  }
  */
  cerr << "Do you want to find clusters with glpk(y/n)?\n";
  char cl;
  cin >> cl;
  if (cl == 'y') {
  //start clustering
  glp_prob *mip = glp_create_prob();
  glp_set_prob_name(mip, "cluster");
  glp_set_obj_dir(mip, GLP_MAX);
  glp_add_cols(mip, N);
  long indx{1};
  //set object
  for (int i=1;i<=N;i++) {
    node_indx[sel_nodes[i-1]->info] = indx;indx++;
    int lv = sel_nodes[i-1]->nrleaves;
    glp_set_col_name(mip, i, sel_nodes[i-1]->info.c_str());
    glp_set_col_kind(mip, i, GLP_BV);
    glp_set_obj_coef(mip, i, lv);
  }
  int *ia = new int [1+MAX_COEFF_NR];
  int *ja = new int [1+MAX_COEFF_NR];
  double *ar = new double [1+MAX_COEFF_NR];
  glp_add_rows(mip, N);
  indx = 1;
  for (int i=1;i<=N;i++) {
    glp_set_row_bnds(mip, i, GLP_DB, 0.0, 1.0);
    ia[indx]=i;ja[indx]=i;ar[indx]=1;indx++;
    if (sel_nodes[i-1] != root) IndexAncestors(sel_nodes[i-1], ia, ja, ar, i, indx);
    }
  auto itq = begin(q);
  auto isel = begin(sel);
  int row_n = N + 1;
  for (int i = 1;i < N;i++) {
    for (int j = 0;j < i;j++) {
      if (*isel) {
	//cout << sel_nodes[i]->info sp sel_nodes[j]->info << '\t' << *itq sp 2*C + FDR - *itq << endl;
	glp_add_rows(mip, 1);//add constraint row
	glp_set_row_bnds(mip, row_n, GLP_UP, 0.0, 2*C + FDR - *itq);
	ia[indx]=row_n;ja[indx]=node_indx[sel_nodes[i]->info];ar[indx]=C;indx++;
	ia[indx]=row_n;ja[indx]=node_indx[sel_nodes[j]->info];ar[indx]=C;indx++;
	row_n++;
	itq++;
      }
      isel++;
    }
  }
  glp_load_matrix(mip, indx-1, ia, ja, ar);
  glp_simplex(mip, NULL);
  glp_iocp parm;
  glp_init_iocp(&parm);
  //parm.presolve = GLP_ON;
  int err = glp_intopt(mip, &parm);
  cout << "Return value (should be 0): " << err << endl;
  cout << "Nr of clustered leaves: " << glp_mip_obj_val(mip) << endl;
  for (int i=1;i<=N;i++) {
    if (glp_mip_col_val(mip, i) == 1 ) {
      cout << glp_get_col_name(mip, i) sp ":" << endl;
      printLeaves(sel_nodes[i-1]);
      cout << endl;
    }
  }
  glp_delete_prob(mip);
  delete[] ia;
  delete[] ja;
  delete[] ar;
  }
  show_event("total time", tm);

}
