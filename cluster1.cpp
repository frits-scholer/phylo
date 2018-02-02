#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

#define all(t) begin(t), end(t)
#define sp << " " <<

int main() {
  ifstream is("constraints1.in");
  if (!is) {
    cerr << "Could not open file\n";
    return 1;
  }
  is.close();
}
