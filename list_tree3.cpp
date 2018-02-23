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
using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<
#define tb << "\t" <<

const float epsilon = 0.00001;

enum Token_value {
  NAME,END,
  LP='(',RP=')',NUMBER=':',COMMA=','
};

Token_value curr_tok=END;
string string_value;
float number_value;

struct node;

typedef vector<float> distribution;
typedef vector<node*> nodelist;

struct node {
  node *parent;
  node *child;
  node *sibling;
  float distance;
  string info;
  bool selected;
  bool isleaf;
  distribution D;//needed for KS
  node(): parent(nullptr), child(nullptr), sibling(nullptr), selected(false),isleaf(false) {}
};

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
      cur_node->child = new node;
      cur_node->child->parent = cur_node;
      cur_node = cur_node->child;
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
      nptr->sibling = cur_node->child;
      cur_node->child = nptr;
      cur_node->child->parent = cur_node;
      cur_node = cur_node->child;
      break;
    case END:
      cout << "error";
    }
  }
  return root;
}
  
void printNodes(node *root) {
  cout << "node:" sp root->info
       tb "ancestor:" sp root->parent->info
       tb "distance:" sp root->distance << endl;
  node *nptr = root->child;
  while (nptr) {
    printNodes(nptr);
    nptr = nptr->sibling;
  }
}

bool is_root(node *root) {
  return root == root->parent;
}

void show_event(string s, clock_t& tm) {
  tm = clock()-tm;
  cerr <<  "\t" << s << " " << (double) tm/CLOCKS_PER_SEC << " s "<< endl;
}

int main() {
  nodelist leaves;
  clock_t tm=clock();
  node *root = build_tree(leaves);
  if (!root) return 1;
  printNodes(root);
  show_event("total time", tm);
}
