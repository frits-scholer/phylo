#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stack>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <cctype>
#include <map>
#include <utility>
#include <iterator>
#include <list>
#include <iomanip>
#include "kstwo.hpp"
#include "bh_fdr.hpp"
using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<
#define tb << "\t" <<
#define tb2 << "\t\t" <<
#define tb3 << "\t\t\t" <<

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
    cerr << "internal error in get_token\n";
    return curr_tok = END;
  }
}

node* build_tree(const string& tree_name, nodevector& leaves) {
  ifstream is(tree_name);//open nwk file
  if (!is) {
    cerr << "Error: could not open tree file\n";
    exit(EXIT_FAILURE);
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
      cerr << "internal error in build_tree";
    }
  }
  if (!root) {
    cerr << "Error: could not construct input tree\n";
    exit(EXIT_FAILURE);
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
  bool nfirst{false};
  for_each(all(root->children),[&](node *nd){
      if (nfirst) os << ','; else nfirst = true;
      append_stream(os, nd);});
  os << ')';
  if (is_root(root)) os << ";";
  else os << ':' << root->distance;
}

void write_newick(const string& tname, node* root) {
  ofstream os(tname);
  if (!os) {
    cerr << "Error: could not open output file\n";
    exit(EXIT_FAILURE);
  }
  append_stream(os, root);
}

int nrLeaves(node *root) {
  if (root->isleaf) {root->nrleaves = 1;return 1;}
  int N=0;
  for_each(all(root->children),[&](node *nd){N +=  nrLeaves(nd);});
  root->nrleaves = N;
  return N;
}

void select_clades(node *root, int cs) {
  if (root->isleaf) return;
  for_each(all(root->children),[&](node *nd){select_clades(nd,cs);});
  if (root->nrleaves >= cs ) root->selected = true;
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
    z->D.push_back(ild);
    if (is_root(z)) return;
    z = z->parent;
  }
}

float calc_distance(node* x, node* y) {
  node* z = common_ancestor(x, y);
  return  ancestor_distance(z, x) + ancestor_distance(z, y);
}

void calc_mean(node *root, vector<node_mean>& v) {
  for_each(all(root->children),[&](node *nd){calc_mean(nd,v);});
  if (!(root->isleaf) && root->selected) {
    float S = accumulate(all(root->D),0.0);
    v.push_back(node_mean(S/(root->D).size(), root));
  }
}

void printLeaves(node *root, node *cl, map<node*,node*>& tc) {
  for_each(all(root->children),[&](node *nd){printLeaves(nd, cl, tc);});
  if (root->isleaf) tc[root] = cl;
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
  if (!(root->isleaf)) {sort(all(root->D));root->selected=false;}
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

void write_ancestors(ostream&os, nodevector& sel_nodes) {
  for_each(all(sel_nodes),[&](node *nd){
      os << 'N' << nd->info << ": ";
      node* z = nd;
      while (!is_root(z)) {
	z = z->parent;
	if (z->selected) os << z->info << " ";
      }
      os << endl;
    });
  os << endl;
}

void write_qvals(ostream&os, nodevector& sel_nodes, vector<float>& q, vector<bool> sel) {
  int N = sel_nodes.size();
  auto it = begin(sel);
  auto itq = begin(q);
  for (int i = 1;i < N;i++) {
    for (int j = 0;j < i;j++) {
      if (*it) {
	os << 'I' << sel_nodes[i]->info << 'J' << sel_nodes[j]->info << 'Q' << *itq << endl;
	itq++;
      }
      it++;
    }
  }
  os << endl;
}

void write_distributions(ostream&os, nodevector& sel_nodes) {
    ostream_iterator<float> osf(os, " ");
    for_each(all(sel_nodes),[&](node* nd){
	os << 'N' << nd->info << ": ";
	copy(all(nd->D),osf);
	os << endl;
      });
    os << endl;
}

void write_internode_distances(ostream&os, nodevector& sel_nodes) {
  int N = sel_nodes.size();
  for (int i = 1;i < N;i++) {
    for (int j = 0;j < i;j++) {
	os << 'I' << sel_nodes[i]->info << 'J' << sel_nodes[j]->info << ": "
	   << calc_distance(sel_nodes[i], sel_nodes[j]) << endl;
    }
  }
  os << endl;
}

int main(int argc, char* argv[]) {
  istream* is;//A Bjarne Stroustrup trick to allow redirection or an input file
  switch(argc) {
  case 1:
    is = &cin;
    break;
  case 2:
    is = new ifstream(argv[1]);
    if (!(*is)) {
      cerr << "Error: could not open input file\n";
      exit(EXIT_FAILURE);
    }
    break;
  default:
    cerr << "too many arguments\n";
    exit(EXIT_FAILURE);
  }
  string tree_name;
  (*is) >> tree_name;
  nodevector leaves;
  node *root = build_tree(tree_name, leaves);
  clock_t tm=clock();
  char zl;
  (*is) >> zl;//input collapse zerolength
  //Start looping to complete the input
  vector<int> cs;vector<float> fdr;vector<float> gamma;
  while (*is) {
    int vcs;
    (*is) >> vcs;
    if (!(*is)) break;
    float vfdr;
    (*is) >> vfdr;
    float vgamma;
    (*is) >> vgamma;
    cs.push_back(vcs);fdr.push_back(vfdr);gamma.push_back(vgamma);
  }
  //end of input
  if (zl == 'y') {
    rzb_nodes(root->children);
    cleanup_zero_nodes();
    rzn(root->children);
    cleanup_merge_nodes();
  }
  nrLeaves(root);
  ladderize(root);
  string tname = "reordered_" + tree_name;
  write_newick(tname, root);
  //end of output reordered tree
  //set up interleaf distances
  for (auto il=begin(leaves)+1;il != end(leaves);il++) {
    for (auto jl = begin(leaves);jl != il;jl++) {
      process_distance(*il,*jl);
    }
  }
  tree_name.erase(tree_name.find('.'));
  //From here on iterate over parameter sets
  ofstream os("tmp.out");
  for (unsigned int t = 0;t < cs.size();t++) {
    select_clades(root, cs[t]);
    vector<node_mean> means;
    calc_mean(root, means);
    sort(all(means));
    size_t msize = means.size();
    float gm;
    if (msize & 1) gm = means[msize/2].first;
    else gm = (means[msize/2-1].first + means[msize/2].first)/2.0;
    auto it = lower_bound(all(means),node_mean(gm,nullptr));
    msize = distance(it, end(means));
    float mad;
    if (msize & 1) mad = (it + msize/2)->first - gm;
    else mad = ((it + msize/2-1)->first + (it + msize/2)->first)/2.0 - gm;
    float upperbound = gm + gamma[t] * mad;
    it = upper_bound(all(means), node_mean(upperbound+0.00000001, nullptr));
    means.erase(it, end(means));
    //presort all data
    preSort(root);//turn selected off
    nodevector sel_nodes;
    for_each(means.rbegin(),means.rend(),[&](node_mean m){
	sel_nodes.push_back(m.second);
	m.second->selected = true;});//set selected true if selected
    vector<bool> sel;
    vector<float> p;
    int N = sel_nodes.size();  
    for (int i = 1;i < N;i++) {
      for (int j = 0;j < i;j++) {
	sel.push_back(false);
	//cout << sel_nodes[i]->info sp sel_nodes[j]->info << '\t';
	/*
	if (is_ancestor(sel_nodes[i], sel_nodes[j]) ||
	    is_ancestor(sel_nodes[j], sel_nodes[i])) {
	  auto ks = kstwo(sel_nodes[i]->D, sel_nodes[j]->D);
	  p.push_back(ks.second);
	  sel.back()=true;
	}
	else */{
	  node* ca = common_ancestor(sel_nodes[i],sel_nodes[j]);
	  //if (ca != sel_nodes[j]->parent) continue; 
	  auto ksi = kstwo(sel_nodes[i]->D, ca->D);
	  auto ksj = kstwo(sel_nodes[j]->D, ca->D);
	  p.push_back(max(ksi.second, ksj.second));
	  sel.back()=true;
	}
      }
    }
    vector<float> q(p.size());
    bh_fdr(p,q);
    write_ancestors(os, sel_nodes);
    write_qvals(os, sel_nodes, q, sel);
    write_distributions(os, sel_nodes);
    write_internode_distances(os, sel_nodes);
  }
  show_event("total time", tm);
}
