#ifndef LINEARSYSSOLVERS_H
#define LINEARSYSSOLVERS_H

#include <vector>

std::vector<double> Trid(const std::vector<std::vector<double> > &, const std::vector<double> &);

std::vector<double> Penta(const std::vector<std::vector<double> > &, const std::vector<double> &);

std::vector<double> Gauss(const std::vector<std::vector<double> > &, const std::vector<double> &);

std::vector<std::vector<double> > Block(const std::vector<std::vector<std::vector<double> > >&, const std::vector<std::vector<std::vector<double> > >&, const std::vector<std::vector<std::vector<double> > >&, const std::vector<std::vector<double> >&);

#endif