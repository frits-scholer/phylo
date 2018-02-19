#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <algorithm>
#include <string>
#include <cctype>
#include <map>
#include <utility>
using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<
const unsigned int MAX_COEFF_NR = 10000000;
const float C = 1000;

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

multimap<float, pair<node*,string>, greater<float>> ordered_nodes;

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

float calc_distance(node* x, node* y) {
  node* z = common_ancestor(x, y);
  return  ancestor_distance(z, x) + ancestor_distance(z, y);
}

void printLeaves(node* realroot, node *root) {
  if (root->ltree) printLeaves(realroot, root->ltree);
  if (root->rtree) printLeaves(realroot, root->rtree);
  if (root->isleaf) {
    float cur_dist = ancestor_distance(realroot, root);
    cout << root->info << '\t' << cur_dist << '\n';
  }
}

void nrLeaves(node *root, int& lv) {
  if (root->ltree) nrLeaves(root->ltree, lv);
  if (root->rtree) nrLeaves(root->rtree, lv);
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

void printAncestors(node *root) {
  node* z = root;
  while (true) {
    z = z->backlink;
    if (z == z->backlink) break;
    cout << z->info << " ";
  }
}

float printNodes(node* realroot, node *root) {
  float lmin{0.0}, rmin{0.0}, node_min;
  if (root->isleaf) node_min = ancestor_distance(realroot, root);
  //calc minimum root to leaf
  else {
    lmin = printNodes(realroot, root->ltree);
    rmin = printNodes(realroot, root->rtree);
    node_min = min(lmin, rmin);
    //cout << root->info << ": " << node_min << endl;
    ordered_nodes.insert(make_pair(node_min, make_pair(root,root->info)));
  }
  cout << root->info << ": " << node_min << endl;
  return node_min;
}

void show_event(string s, clock_t& tm) {
  tm = clock()-tm;
  cerr <<  "\t" << s << " " << (double) tm/CLOCKS_PER_SEC << " s "<< endl;
}

int main() {
  vector<node*> leaves;
  clock_t tm=clock();
  node* root = build_tree(leaves);
  if (!root) return 1;
  multimap<float, string, greater<float>> ordered_leaves;
  //cout << "leaf order:\n";
  for (auto leaf : leaves) {
    //cout << leaf->info << '\t' << ancestor_distance(root, leaf) << endl;
    ordered_leaves.insert(make_pair(ancestor_distance(root, leaf), leaf->info));
  }
  cout << "\nLeaves in branch order:\n";
  for (auto it = ordered_leaves.begin();it != ordered_leaves.end();++it)
    cout << it->second << '\t' << it->first << endl;
  /*
  */

  cout << "\nNodes in visiting order:\n";
  printNodes(root, root);
  cout << "\nNodes in branch order:\n";
  for (auto it = ordered_nodes.begin();it != ordered_nodes.end();++it)
    cout << it->second.second << '\t' << it->first << endl;
  //cout << "\nLeaves In order:\n";
  //printLeaves(root, root);
  show_event("total time", tm);
}
