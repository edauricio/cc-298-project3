#include <vector>
#include <iostream>
#include <string>

using namespace std;

enum class EnumType { a = 1 , b };

int main() {
  EnumType test = EnumType::a;
  cout << EnumType::a << endl;

  return 0;
}