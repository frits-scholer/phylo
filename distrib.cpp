#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <algorithm>
#include <string>
#include <cctype>
using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<

enum Token_value {
  NORMAL,NAME,NUMBER,END,
  LP='(',RP=')',COLON=':',COMMA=','
};

Token_value curr_tok=NORMAL;
string string_value;
double number_value;
int ml;

typedef vector<double> distribution;

struct node {
  node* ltree;
  node* rtree;
  node* backlink;
  double distance;
  string info;
  bool selected;
  bool isleaf;
  double sum_ild;//current sum of interleaf distances
  int count_ild;
  vector<double> D;//needed for KS
  node(): ltree(nullptr), rtree(nullptr),backlink(nullptr),isleaf(false),
	  sum_ild(0.0), count_ild(0){}
};

typedef pair<double, node*> node_mean;

Token_value get_token(istream& is) {
  char ch;
  do {
    if (!is.get(ch)) return curr_tok = END;
  }  while (ch !=';' && isspace(ch));
  switch(ch) {
  case '(':case ')':case ':':case ',':
    return curr_tok = Token_value(ch);
  case ';':
    return curr_tok = END;      
  case '0':case '1':case '2':case '3':case '4':case '5':case '6':
    case '7':case '8':case '9':
    is.putback(ch);
    is >> number_value;
    return curr_tok = NUMBER;
  case '\'':
    string_value.clear();
    while (is.get(ch)) {if (ch == '\'') break;string_value.push_back(ch);}
    return curr_tok = NAME;
  default:
  return curr_tok = NORMAL;
  }
}

node* build_tree(vector<node*>& leaves) {
  cerr << "filename?\n";
  string fname;
  cin >> fname;
  ifstream is(fname.c_str());
  //ifstream is("toy_tree_small.nwk");
  //ifstream is("test.nwk");
  //ifstream is("7330_leaves.nwk");
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
      break;
    case NUMBER:
      cur_node->distance = number_value;
      break;
    case COMMA:
      cur_node = A.top();
      cur_node->rtree = new node();
      cur_node->rtree->backlink = cur_node;
      cur_node = cur_node->rtree;
      break;
    case NORMAL:case END:
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

double ancestor_distance(node* z, node* w) {//w is descendant of z
  double dist = 0;
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

node* common_ancestor(node* x, node* y) {//y occurs before x
  node* z = x->backlink;
  while (!search(z,y)) {if (z==z->backlink) break;else z = z->backlink;}
  return z;
}

double distance(node* x, node* y) {//y occurs before x
  node* z = common_ancestor(x, y);
  return ancestor_distance(z, x) + ancestor_distance(z, y);
}

void process_distance(node* x, node* y) {//y occurs before x
  node* z = common_ancestor(x, y);
  double ild = ancestor_distance(z, x) + ancestor_distance(z, y);
  while (true) {//skip small clades
    if (z->selected) break;
    if (z == z->backlink) break;//root is always selected
    z = z->backlink;
  }
  //from here on up all clades are selected
  while (true) {//
    z->sum_ild += ild;
    z->count_ild++;
    if (z == z->backlink) break;//root is always selected
    z = z->backlink;
  }
}

void inOrder(node *root, vector<node_mean>& v) {
  if (root->ltree) inOrder(root->ltree, v);
  if (!(root->isleaf) && root->selected)
      v.push_back(node_mean(root->sum_ild/root->count_ild, root));
  if (root->rtree) inOrder(root->rtree, v);
}

void clear(node *root) {
  if (root->ltree) clear(root->ltree);
  if (root->isleaf || root->selected)
    root->selected = false;
  if (root->rtree) clear(root->rtree);
}

void preOrder(node *root) {
  if (root->ltree) preOrder(root->ltree);
  if (root->rtree) preOrder(root->rtree);
  if (root->isleaf) root->selected = true;
}

void process_distance2(node* x, node* y) {//y occurs before x
  node* z = common_ancestor(x, y);
  double ild = ancestor_distance(z, x) + ancestor_distance(z, y);
  while (true) {//skip unselected clades
    if (z->selected) break;
    if (z == z->backlink) break;//root is always selected
    z = z->backlink;
  }
  //from here on up all clades are selected
  while (true) {//
    if (z->selected) (z->D).push_back(ild);
    if (z == z->backlink) break;//root is always selected
    z = z->backlink;
  }
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
  cout << "Minimum nr of leaves:" sp ml sp "Gamma: " sp gamma << endl;
  select_clades(root);
  root->selected = true;
  //set up interleaf distances
  for (auto il=begin(leaves)+1;il != end(leaves);il++) {
    for (auto jl = begin(leaves);jl != il;jl++) {
      process_distance(*il,*jl);
    }
  }
  vector<node_mean> means;
  inOrder(root, means);
  sort(all(means));
  size_t msize = means.size();
  double gm;
  if (msize & 1) gm = means[msize/2].first;
  else gm = (means[msize/2-1].first + means[msize/2].first)/2.0;
  cout << "size of means = " sp msize sp "grand median = " << gm << endl;
  auto it = lower_bound(all(means),node_mean(gm,nullptr));
  msize = distance(it, end(means));
  double mad;
  if (msize & 1) mad = (it + msize/2)->first - gm;
  else mad = ((it + msize/2-1)->first + (it + msize/2)->first)/2.0 - gm;
  cout << "mad = " << mad << endl;
  double upperbound = gm + gamma * mad;
  it = upper_bound(all(means), node_mean(upperbound+0.000001, nullptr));
  means.erase(it, end(means));
  cout << "size of filtered means = " << means.size() << endl;
  clear(root);//turn off the whole tree
  for (auto m: means) {
    m.second->selected = true;//turn on selected nodes
    preOrder(m.second);//turn on selected leaves
  }
  for (auto il=begin(leaves)+1;il != end(leaves);il++) {
    if (!((*il)->selected)) continue;
    for (auto jl = begin(leaves);jl != il;jl++) {
      if (!((*jl)->selected)) continue;
      process_distance2(*il,*jl);//add the interleaf distance to the selected nodes
    }
  }
  for (auto m: means) {
    cout << m.second->info << ": ";
    for (double d: m.second->D) cout << d << " ";
    cout << endl;
  }
  cout << "finished!" << endl;
}
