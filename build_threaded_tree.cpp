#include <iostream>
#include <fstream>
#include <stack>
#include <algorithm>
#include <string>
#include <cctype>

using namespace std;

enum Token_value {
  NORMAL,NAME,NUMBER,END,
  LP='(',RP=')',COLON=':',COMMA=','
};

Token_value curr_tok=NORMAL;
string string_value;
double number_value;

struct node {
  node* ltree;
  node* rtree;
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

void inOrder(node *root) {
    if (root->ltree) inOrder(root->ltree);
    if (root->isleaf) cout << root->info << " ";
    cout << ": " << root->distance << '-' << root->backlink << endl;
    if (root->rtree) inOrder(root->rtree);
}

int main() {
  /*
  cout << "filename? ";
  string fname;
  cin >> fname;
  ifstream is(fname.c_str());
  */
  ifstream is("toy_tree_small.nwk");
  if (!is) {
    cerr << "could not open file\n";
    return 1;
  }
  int level(0);
  node *root = new node();
  root->backlink = root;
  root->info = "root";
  node *cur_node = root;
  stack<node *> A;
  while (is) {
    get_token(is);
    if (curr_tok == END) break;
    switch(curr_tok) {
    case LP:
      level++;
      A.push(cur_node);
      cur_node->ltree = new node();
      cur_node->ltree->backlink = cur_node;
      cur_node = cur_node->ltree;
      break;
    case RP:
      level--;
      cur_node = A.top();
      A.pop();
      break;
    case NAME:
      cur_node->info = string_value;
      cur_node->isleaf = true;
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
  cout << cur_node->info << endl;
  inOrder(root);
  cout << endl;
}
