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

struct node {
  node * ltree;
  node * rtree;
  node* backlink;
  double distance;
  string info;
  bool isleaf;
  node(): ltree(nullptr), rtree(nullptr),backlink(nullptr),isleaf(false) {}
};

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
  /*
  cout << "filename? ";
  string fname;
  cin >> fname;
  ifstream is(fname.c_str());
  */
  ifstream is("toy_tree_small.nwk");
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


void inOrder(node *root) {
    if (root->ltree) inOrder(root->ltree);
    /*if (root->isleaf)*/ cout << root->info << " ";
    cout << ": " << root->distance << " <- " << root->backlink->info << endl;
    if (root->rtree) inOrder(root->rtree);
}

int main() {
  vector<node*> leaves;
  node* root = build_tree(leaves);
  //root->info = "root";
  if (!root) return 1;
  //for_each(all(leaves),[](const node* l){cout << l->info sp l->distance << endl;});
  inOrder(root);
  cout << endl;

  for (auto il=begin(leaves)+1;il != end(leaves);il++) {
    cout << (*il)->info << ':';
    for (auto jl = begin(leaves);jl != il;jl++) {
      cout << (*jl)->info << (char)LP << distance(*il,*jl) << (char)RP << ' ';
    }
    cout << endl;
  }

}
