#ifndef NUMERICS_H
#define NUMERICS_H

#include <string>

class Method {
  public:
    Method() = default;
    explicit Method(std::string sFile) : File(sFile) { parseScheme(File); };
    void parse(std::string sFile) { File = sFile; parseScheme(File); };
    std::string &name() { return Name; };
    std::string &expimp() { return ExpImp; };
    std::string &file() { return File; };
    int &timeOrder() { return TimeOrder; };
    int &spaceOrder() { return SpaceOrder; };
  private:
    std::string Name, ExpImp, File;
    int TimeOrder = 1, SpaceOrder = 1;
    void parseScheme(std::string);
};

#endif