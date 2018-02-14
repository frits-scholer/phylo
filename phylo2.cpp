#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <algorithm>
#include <string>
#include <cctype>
#include <map>
#include "kstwo.hpp"
#include "bh_fdr.hpp"
using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<

enum Token_value {
  NAME,END,
  LP='(',RP=')',COLON=':',COMMA=','
};

Token_value curr_tok=END;
string string_value;
float number_value;
int ml;

typedef vector<float> distribution;
map<string,int> node_indx;

struct node {
  node* ltree;
  node* rtree;
  node* backlink;
  float distance;
  string info;
  bool selected;
  bool isleaf;
  vector<float> D;//needed for KS
  node(): ltree(nullptr), rtree(nullptr),backlink(nullptr),isleaf(false){}
};

typedef pair<float, node*> node_mean;
typedef vector<node*> nodelist;

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
    return curr_tok = COLON;
  case '\'':
    string_value.clear();
    while (is.get(ch)) {if (ch == '\'') break;string_value.push_back(ch);}
    return curr_tok = NAME;
  default:
    is.putback(ch);
    string_value.clear();
    if (curr_tok == LP || curr_tok == COMMA) {
      getline(is, string_value,':');is.putback(':');return curr_tok = NAME;}
    cerr << " ERROR\n";
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
  node *root = new node();
  node *cur_node = root;
  root->backlink = root;
  int internal_count = 0;
  stack<node *> A;
  while (is) {
    get_token(is);
    if (curr_tok == END) break;
    switch(curr_tok) {
    case LP:
      A.push(cur_node);
      cur_node->info = to_string(internal_count);
      internal_count++;
      cur_node->ltree = new node();
      cur_node->ltree->backlink = cur_node;
      cur_node = cur_node->ltree;
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
    case COLON:
      cur_node->distance = number_value;
      break;
    case COMMA:
      cur_node = A.top();
      cur_node->rtree = new node();
      cur_node->rtree->backlink = cur_node;
      cur_node = cur_node->rtree;
      break;
    case END:
      cout << "error";
    }
  }
  return root;
}

int select_clades(node *root) {//returns nr of leaves
  if (root->isleaf) return 1;
  int nr_leaves = select_clades(root->ltree) + select_clades(root->rtree);
  if (nr_leaves >= ml ) {
    //cout << root->info << '[' << nr_leaves << ']' << endl;
    root->selected = true;
  }
  else root->selected = false;
  return nr_leaves;
}

float ancestor_distance(node* z, node* w) {//w is descendant of z
  float dist = 0;
  while (w != z) {
    dist += w->distance;
    w = w->backlink;
  }
  return dist;
}

bool search(node* root, node* child) {
  if (root->ltree == child) return true;
  if (root->rtree == child) return true;
  if (root->ltree) return search(root->ltree, child);
  if (root->rtree) return search(root->rtree, child);
  return false;
}

node* common_ancestor(node* x, node* y) {
  node* z = x->backlink;
  while (!search(z,y)) {
    if (z==z->backlink) break;
    else z = z->backlink;
  }
  return z;
}

void process_distance(node* x, node* y) {
  node* z = common_ancestor(x, y);
  float ild = ancestor_distance(z, x) + ancestor_distance(z, y);
  while (true) {
    z->D.push_back(ild);//generate D even if internal node is not selected
    if (z == z->backlink) return;
    z = z->backlink;
  }
}

void calc_mean(node *root, vector<node_mean>& v) {
  if (root->ltree) calc_mean(root->ltree, v);
  if (!(root->isleaf) && root->selected) {
    float S = accumulate(all(root->D),0.0);
    v.push_back(node_mean(S/(root->D).size(), root));
  }
  if (root->rtree) calc_mean(root->rtree, v);
}

void printLeaves(node *root, int& lv) {
  if (root->ltree) printLeaves(root->ltree, lv);
  if (root->rtree) printLeaves(root->rtree, lv);
  if (root->isleaf) {lv++;/*cout << root->info << ' ';*/}
}

bool is_ancestor(node* x, node* y) {//x is ancestor of y
  node* z = y->backlink;
  do {
    if (x == z) return true;
    z = z->backlink;
  } while (z != z->backlink);
  return x==z;
} 

void preSort(node *root) {
  if (root->ltree) preSort(root->ltree);
  if (root->rtree) preSort(root->rtree);
  if (!(root->isleaf)) sort(all(root->D));
}

void printIndexAncestors(node *root) {
  node* z = root;
  while (true) {
    z = z->backlink;
    if (z == z->backlink) break;
    if (node_indx[z->info]) cout << z->info << " ";
  }
}

void show_event(string s, clock_t& tm) {
  tm = clock()-tm;
  cerr <<  "\t" << s << " " << (double) tm/CLOCKS_PER_SEC << " s "<< endl;
}

int main() {
  vector<node*> leaves;
  node* root = build_tree(leaves);
  if (!root) return 1;
  int gamma;
  cerr << "Minimum nr of leaves of considered clades?\n";
  cin >> ml;
  cerr << "Gamma factor?\n";
  cin >> gamma;
  //cout << "Minimum nr of leaves:" sp ml sp "Gamma: " sp gamma << endl;
  //start timer

  clock_t tm=clock();
  cout << leaves.size() << endl;

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
  size_t msize = means.size();
  float gm;
  if (msize & 1) gm = means[msize/2].first;
  else gm = (means[msize/2-1].first + means[msize/2].first)/2.0;
  //cout << "size of means = " sp msize sp "grand median = " << gm << endl;
  auto it = lower_bound(all(means),node_mean(gm,nullptr));
  msize = distance(it, end(means));
  float mad;
  if (msize & 1) mad = (it + msize/2)->first - gm;
  else mad = ((it + msize/2-1)->first + (it + msize/2)->first)/2.0 - gm;
  //cout << "mad = " << mad << endl;
  float upperbound = gm + gamma * mad;
  it = upper_bound(all(means), node_mean(upperbound+0.00000001, nullptr));
  means.erase(it, end(means));
  //cout << "size of filtered means = " << means.size() << endl;
  cout << means.size() << endl;
  //presort all data
  preSort(root);
  nodelist sel_nodes;
  for_each(means.rbegin(),means.rend(),[&](node_mean m){sel_nodes.push_back(m.second);});

  int indx{1};
  for (auto n: sel_nodes) {
    cout << 'n' << n->info << " ";
    node_indx[n->info] = indx;indx++;
    int lv{0};
    printLeaves(n, lv);
    cout << lv << endl;
  }
  for (unsigned int i = 0;i < sel_nodes.size();i++) {
    cout << sel_nodes[i]->info << ": ";
    printIndexAncestors(sel_nodes[i]);
    cout << endl;
  }

  vector<bool> sel;
  vector<float> p;
  for (unsigned int i = 1;i < sel_nodes.size();i++) {
    for (unsigned int j = 0;j < i;j++) {
      sel.push_back(false);
      //cout << sel_nodes[i]->info sp sel_nodes[j]->info << '\t';
      if (is_ancestor(sel_nodes[i], sel_nodes[j]) || is_ancestor(sel_nodes[j], sel_nodes[i])) {
	auto ks = kstwo(sel_nodes[i]->D, sel_nodes[j]->D);
	p.push_back(ks.second);
	sel.back()=true;
	}
      else {
	node* ca = sel_nodes[i]->backlink;
	if (ca != sel_nodes[j]->backlink) continue; 
	//node* ca = common_ancestor(sel_nodes[i], sel_nodes[j]);
	auto ksi = kstwo(sel_nodes[i]->D, ca->D);
	auto ksj = kstwo(sel_nodes[j]->D, ca->D);
	p.push_back(max(ksi.second, ksj.second));
	sel.back()=true;
      }
    }
  }

  vector<float> q(p.size());
  bh_fdr(p,q);
  auto isel = begin(sel);
  auto itq = begin(q);
  for (unsigned int i = 1;i < sel_nodes.size();i++) {
    for (unsigned int j = 0;j < i;j++) {
      if (*isel) {
	if (*itq > 0.9) cout << sel_nodes[i]->info sp sel_nodes[j]->info << '\t' << *itq << endl;
	itq++;
      }
      isel++;
    }
  }

  show_event("total time", tm);
}
