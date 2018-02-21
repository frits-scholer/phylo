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
const float epsilon = 0.00001;

enum Token_value {
  NAME,END,
  LP='(',RP=')',NUMBER=':',COMMA=','
};

Token_value curr_tok=END;
string string_value;
float number_value;
int ml;
struct node;

typedef vector<float> distribution;
typedef vector<node*> nodelist;
map<string,int> node_indx;


struct node {
  nodelist subtree;
  node* backlink;
  float distance;
  string info;
  bool selected;
  bool isleaf;
  distribution D;//needed for KS
  node(): backlink(nullptr),isleaf(false){}
};

typedef pair<float, node*> node_mean;

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

node* build_tree(nodelist& leaves) {
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
      cur_node->subtree.push_back(new node());
      cur_node->subtree.back()->backlink = cur_node;
      cur_node = cur_node->subtree.back();
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
      cur_node->subtree.push_back(new node());
      cur_node->subtree.back()->backlink = cur_node;
      cur_node = cur_node->subtree.back();
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
  if (find(all(root->subtree), child) != end(root->subtree)) return true;
  for (auto stree: root->subtree)
    if (stree) return search(stree, child);
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
  for (auto stree: root->subtree)
    if (stree) printLeaves(realroot, stree);
  if (root->isleaf) {
    float cur_dist = ancestor_distance(realroot, root);
    cout << root->info << '\t' << cur_dist << '\n';
  }
}

void nrLeaves(node *root, int& lv) {
  for (auto stree: root->subtree)
    if (stree) nrLeaves(stree, lv);
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

void printNodes(node *root) {
  for (auto stree: root->subtree)
    if (stree) printNodes(stree);
  cout << root->info << ": " << root->distance << endl;
}

void collapse(node *root) {
  for (auto stree: root->subtree)
    collapse(stree);
  if (root->distance > epsilon) return;
  
  node * r = root->backlink;
  if (r != r->backlink) {//take out subtree if parent is not the real root
    r->subtree.erase(find(all(r->subtree),root));
    //r->backlink->subtree.push_back(root);
  }
  //else r->subtree.push_back(root);
  root->distance = r->distance;
  cout << root->info sp r->info sp r->subtree.size() << endl;//if r not real root make root a sibling of r

}

void show_event(string s, clock_t& tm) {
  tm = clock()-tm;
  cerr <<  "\t" << s << " " << (double) tm/CLOCKS_PER_SEC << " s "<< endl;
}

int main() {
  nodelist leaves;
  clock_t tm=clock();
  node* root = build_tree(leaves);
  if (!root) return 1;
  collapse(root);
  printNodes(root);
  show_event("total time", tm);
}
