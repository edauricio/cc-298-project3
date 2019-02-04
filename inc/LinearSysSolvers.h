#ifndef LINEARSYSSOLVERS_H
#define LINEARSYSSOLVERS_H

#include <vector>

template <typename T>
std::vector<T> Trid(const std::vector<std::vector<T> > &, const std::vector<T> &);

template <typename T>
std::vector<T> Penta(const std::vector<std::vector<T> > &, const std::vector<T> &);

template <typename T>
std::vector<T> Gauss(const std::vector<std::vector<T> > &, const std::vector<T> &);

template <typename T>
std::vector<std::vector<T> > Block(const std::vector<std::vector<std::vector<T> > >&, const std::vector<std::vector<std::vector<T> > >&, const std::vector<std::vector<std::vector<T> > >&, const std::vector<std::vector<T> >&);

#endif