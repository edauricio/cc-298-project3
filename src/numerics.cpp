#include <fstream>
#include "Numerics.h"
#include "json.hpp"

using json = nlohmann::json;

void Method::parseScheme(std::string file) {
  json Data;
  std::ifstream inFile(file);
  Data = json::parse(inFile);

  Name = Data["Method"]["Scheme"];

  if (Data["Method"]["TimeOrder"] == "first")
    TimeOrder = 1;
  else if (Data["Method"]["TimeOrder"] == "second")
    TimeOrder = 2;

  if (Data["Method"]["SpaceOrder"] == "first")
    SpaceOrder = 1;
  else if (Data["Method"]["SpaceOrder"] == "second")
    SpaceOrder = 2;

  inFile.close();
}