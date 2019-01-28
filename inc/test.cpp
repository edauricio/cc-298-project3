#include <iostream>
#include <fstream>
#include "json.hpp"

using json = nlohmann::json;

int main() {

  json j;
  std::ifstream in("numerics");
  j = json::parse(in);

  std::string test = j["object"]["currency"];
  std::cout << test << std::endl;

}