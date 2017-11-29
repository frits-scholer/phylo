#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <iterator>
#include <cctype>

using namespace std;
#define all(t) begin(t), end(t)
#define sp << ' ' <<

enum Token_value {
  NORMAL,NAME,NUMBER,END,
  LP='(',RP=')',COLON=':',COMMA=','
};

Token_value curr_tok=NORMAL;
string string_value;
double number_value;
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
  return NORMAL;
  }
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
  cout << level << endl;
  while (is) {
    get_token(is);
    if (curr_tok == END) break;
    switch(curr_tok) {
    case LP:
      level++;
      cout << "level " << level << endl;
      break;
    case RP:
      level--;
      cout << "level" << level;
      break;
    case NAME:
      cout << string_value << ' ';
      break;
    case COLON:
      cout << ':';
      break;
    case NUMBER:
      cout << number_value << endl;
      break;
    case COMMA:
      cout << ',';
      break;
    case NORMAL:case END:
      cout << "error";
    }
  }
  cout << endl;
}
