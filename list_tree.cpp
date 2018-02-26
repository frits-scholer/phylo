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

bool less_dist(const node* a, const node* b) {
  return a->distance < b->distance;
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
  
int nr_of_children(node *root) {
  int n = 0;
  node *nptr = root->child;
  while (nptr) {
    n++;
    nptr = nptr->sibling;
  }
  return n;
}
  
void printNodes(node *root) {
  int n = nr_of_children(root);
  cout << "node:" sp root->info
    tb "ancestor:" sp root->parent->info
    tb "distance:" sp root->distance
    tb "children:" sp n << endl;
  node *nptr = root->child;
  while (nptr) {
    printNodes(nptr);
    nptr = nptr->sibling;
  }
}

bool is_root(node *root) {
  return root == root->parent;
}

void rzb(node *root) {
  node *nptr = root->child;
  while (nptr) {
    node *sptr = nptr->sibling;//this might be invalidated
    rzb(nptr);
    nptr = sptr;
  }
  if (root->isleaf) {
    if (is_root(root->parent) || root->distance > epsilon) return;
    node *pptr = root->parent;
    node *cptr = pptr->child;
    if (cptr == root) {pptr->child = root->sibling;}//root is not child anymore
    else {
      node *dptr;
      do {
	dptr = cptr;
	cptr = cptr->sibling;
      } while (cptr != root);
      dptr->sibling = root->sibling;
    }
    cptr = pptr->parent->child;
    pptr->parent->child = root;
    root->parent = pptr->parent;//root is grandchild now
    root->sibling = cptr;
    root->distance = pptr->distance<=epsilon?0:pptr->distance;
  }
}

void rzb_nodes(node *root) {
  node *nptr = root->child;
  while (nptr) {
    node *sptr = nptr->sibling;//this might be invalidated
    rzb_nodes(nptr);
    nptr = sptr;
  }
  //a non-leaf
  if (!root->isleaf) {
    if (is_root(root) || root->distance > epsilon) return;
    root->distance = 0;
    //children become grandchildren
    node *pptr = root->parent;
    node *cptr = pptr->child;
    node *dptr = root->child;
    if (!dptr) return;
    pptr->child = root->child;//first child becomes first grandchild
    node *eptr;
    while (dptr) {
      eptr = dptr;
      dptr->parent = pptr;//siblings of first child become grandchildren
      dptr = dptr->sibling;
    }
    eptr->sibling = cptr;//new siblings become siblings of old siblings
    root->child = nullptr;
  }
}

void sort_by_distance(node *root) {
  if (root->isleaf) return;
  //root is not a leaf
  node *cptr = root->child;
  while (cptr) {
    sort_by_distance(cptr);
    cptr = cptr->sibling;
  }
  cptr = root->child;
  list<node*> a_list;
  while (cptr) {
    a_list.push_back(cptr);
    cptr = cptr->sibling;
  }
  a_list.sort(less_dist);

  auto it = begin(a_list);
  root->child = *it;
  while (it != end(a_list)) {
    auto jt = it;
    it++;
    if (it != end(a_list)) (*jt)->sibling = *it;
    else (*jt)->sibling = nullptr;
  }

}

void rzn(node *root) {
  if (root->isleaf) return;
  //root is not a leaf
  node *cptr = root->child;
  while (cptr) {
    node *dptr = cptr->sibling;
    rzn(cptr);
    cptr = dptr;
  }
  //at this point zero and one nodes of children are removed
  int nch = nr_of_children(root);
  if (nch > 1) return;//nothing more to do
  node *pptr = root->parent;
  node *dptr;
  cptr = pptr->child;
  if (nch == 0) {
	if (cptr == root) {//reset child
	  //cerr << "C";
	  pptr->child = root->sibling;
	  delete root;
	  return;
	}
	while (cptr != root) {
	  dptr = cptr;
	  cptr = cptr->sibling;
	}
	dptr->sibling = root->sibling;
	delete root;
	return;
  }
  //case root has single child
  dptr = root->child;
  pptr->child = dptr;//child becomes grandchild
  dptr->sibling = cptr;//new child becomes sibling of old siblings
  dptr->parent = pptr;
  dptr->distance = dptr->distance + root->distance;
  while (cptr != root) {
    dptr = cptr;
    cptr = cptr->sibling;
  }
  dptr->sibling = root->sibling;
  delete root;
}

void show_event(string s, clock_t& tm) {
  tm = clock()-tm;
  cerr <<  "\t" << s << " " << (double) tm/CLOCKS_PER_SEC << " s "<< endl;
}

void append_stream(ostream& os, node* root) {
  if (root->isleaf) {os << root->info << ':' << root->distance;return;}
  os << '(';
  node *nptr = root->child;
  while (nptr) {
    append_stream(os, nptr);
    nptr = nptr->sibling;
    if (nptr) os << ',';
    else os << "):" << root->distance;
  }
  if (is_root(root)) os << ')';
}

void write_newick(node *root) {
  cerr << "filename?\n";
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

int main() {
  nodelist leaves;
  clock_t tm=clock();
  node *root = build_tree(leaves);
  if (!root) return 1;
  rzb(root);
  rzb_nodes(root);
  rzn(root);
  sort_by_distance(root);
  printNodes(root);
  write_newick(root);
  show_event("total time", tm);
}
