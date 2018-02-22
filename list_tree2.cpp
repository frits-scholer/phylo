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
  nodelist subtree;
  node *backlink;
  float distance;
  string info;
  bool selected;
  bool isleaf;
  distribution D;//needed for KS
  node(): backlink(nullptr), selected(false),isleaf(false) {}
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
      cur_node->subtree.push_back(new node);
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
      cur_node->subtree.push_back(new node);
      cur_node->subtree.back()->backlink = cur_node;
      cur_node = cur_node->subtree.back();
      break;
    case END:
      cout << "error";
    }
  }
  return root;
}
  
void printNodes(node *root) {
  if (!root->isleaf) cout << "node:" sp root->info
       tb "ancestor:" sp root->backlink->info
       tb "distance:" sp root->distance
       tb "subtrees:" sp root->subtree.size() << endl;
  else cout << root->info
       tb "ancestor:" sp root->backlink->info
	 tb "distance:" sp root->distance << endl;
  for (auto stree: root->subtree)
    if (stree) printNodes(stree);
}

bool is_root(node *root) {
  return root == root->backlink;
}

void show_event(string s, clock_t& tm) {
  tm = clock()-tm;
  cerr <<  "\t" << s << " " << (double) tm/CLOCKS_PER_SEC << " s "<< endl;
}

void remove_zero_branches(node *n) {
  for (auto child : n->subtree) remove_zero_branches(child);//postorder
  if (n->isleaf) {
    if (n->distance <= epsilon) n->selected = true;
    return;
  }
  //n is an internal node from here on
  if (n->distance <= epsilon) n->selected = true;
  else {
    if (n->subtree.size() <= 1) n->selected = true;
  }
  for (auto it = begin(n->subtree); it != end(n->subtree);it++) {
    auto child = *it;
    if (child->selected) {
      //copy(all(child->subtree), back_inserter(n->subtree));
      cerr << child->info << endl;
      for (auto grandchild : child->subtree) n->subtree.push_back(grandchild);
      delete child;
      //delete child from subtree
      n->subtree.erase(it);
    }
    else {
      for ( auto jt = begin(child->subtree); jt != end(child->subtree);) {
	auto grandchild = *jt;
	if (grandchild->selected) {
	  n->subtree.push_back(grandchild);
	  //delete grandchild from child subtree
	  child->subtree.erase(jt);
	}
	else jt++;
      }
      it++;
    }
  }  
}

int main() {
  nodelist leaves;
  clock_t tm=clock();
  node *root = build_tree(leaves);
  if (!root) return 1;
  printNodes(root);
  cout << "------------------------------------------------------------------------\n";
  remove_zero_branches(root);
  printNodes(root);
  show_event("total time", tm);
}
