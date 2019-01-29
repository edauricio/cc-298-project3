#ifndef OPERATORS_H
#define OPERATORS_H

#include <vector>
#include <iostream>

using namespace std;

/*Operator << on an object of vector class:
Returns the column j of matrix A.
Usage:
A<<j
*/
template <typename T>
vector<T> operator<< (const vector<vector<T> > &m, const size_t &j) {
	if (j >= m[0].size()) {
		cout << "Operator Matrix<<col called on a column beyond matrix maximum. Exiting. . ." << endl;
		exit(-1);
	} else {
		vector<T> temp(m.size());
		for (auto i = 0; i != m.size(); ++i)
			temp[i] = m[i][j];
		return temp;
	}
}

/*Operator >> on an object of vector<vector<vector> > class:
Returns the column j of matrix of vectors A[i][j][k], i.e., returns a vector<vector> of size [i][k].
Usage:
A>>j
*/
template <typename T>
vector<vector<T> > operator>> (const vector<vector<vector<T> > > &mv, const size_t &j) {
	if (j >= mv[0].size()) {
		cout << "Operation Matrix>>col called on a column beyond matrix maximum. Exiting. . . " << endl;
		cout << mv[0].size() << " " << j << endl;
		exit(-1);
	} else {
		vector<vector<T> > temp(mv.size());
		for (auto &v : temp)
			v.resize(mv[0][0].size());
		for (auto i = 0; i != temp.size(); ++i)
			for (auto k = 0; k != temp[i].size(); ++k)
				temp[i][k] = mv[i][j][k];
		return temp;
	}
}

template <typename T>
vector<T> operator+ (const vector<T> &v1, const vector <T> &v2) {
	if (v1.size() == v2.size()) {
		vector<T> temp(v1.size());
		for (auto i = 0; i != v1.size(); ++i)
			temp[i] = v1[i]+v2[i];
		return temp;
	}	else {
		cout << "Operation Vector+Vector called on vectors of different size. Exiting. . ." << endl;
		exit(-1);
	}
}

template <typename T>
vector<vector<T> > operator+ (const vector<vector<T> > &m1, const vector<vector<T> > &m2) {
	if (m1.size() != m2.size() || m1[0].size() != m2[0].size()) {
		cout << "Operation Matrix+Matrix called on matrices of different size. Exiting. . ." << endl;
		exit(-1);
	} else {
		vector<vector<T> > temp(m1.size());
		for (auto i = 0; i != temp.size(); ++i)
			temp[i].resize(m1[0].size());
		for (auto i = 0; i != temp.size(); ++i)
			for (auto j = 0; j != temp[i].size(); ++j)
				temp[i][j] = m1[i][j] + m2[i][j];
		return temp;
	}
}

template <typename T>
vector<T> operator- (const vector<T> &v1, const vector<T> &v2) {
	if (v1.size() != v2.size()) {
		cout << "Operation Vector+Vector called on vectors of different size. Exiting. . ." << endl;
		exit(-1);
	} else {
		vector<T> temp(v1.size());
		for (auto i = 0; i != v1.size(); ++i)
			temp[i] = v1[i]-v2[i];
		return temp;
	}
}

template <typename T>
vector<vector<T> > operator- (const vector<vector<T> > &m1, const vector<vector<T> > &m2) {
	if (m1.size() != m2.size() || m1[0].size() != m2[0].size()) {
		cout << "Operation Matrix+Matrix called on matrices of different size. Exiting. . ." << endl;
		exit(-1);
	} else {
		vector<vector<T> > temp(m1.size());
		for (auto i = 0; i != temp.size(); ++i)
			temp[i].resize(m1[0].size());
		for (auto i = 0; i != temp.size(); ++i)
			for (auto j = 0; j != temp[i].size(); ++j)
				temp[i][j] = m1[i][j] - m2[i][j];
		return temp;
	}
}

template <typename T>
vector<T> operator* (const double &s, const vector<T> &v) {
	vector<T> temp(v.size());
	for (auto i = 0; i != temp.size(); ++i)
		temp[i] = s*v[i];
	return temp;
}

template <typename T>
vector<T> operator* (const vector<T> &v, const double &s) {
  vector<T> temp(v.size());
  for (auto i = 0; i != temp.size(); ++i)
    temp[i] = s*v[i];
  return temp;
}

template <typename T>
vector<vector<T> > operator* (const double &s, const vector<vector<T> > &m) {
	vector<vector<T> > temp(m.size());
	for (auto &v : temp)
		v.resize(m[0].size());
	for (auto i = 0; i != m.size(); ++i)
		for (auto j = 0; j != m[i].size(); ++j)
			temp[i][j] = s*m[i][j];
	return temp;
}

template <typename T>
vector<vector<T> > operator* (const vector<vector<T> > &m, const double &s) {
  vector<vector<T> > temp(m.size());
  for (auto &v : temp)
    v.resize(m[0].size());
  for (auto i = 0; i != m.size(); ++i)
    for (auto j = 0; j != m[i].size(); ++j)
      temp[i][j] = s*m[i][j];
  return temp;
}

template <typename T>
vector<T> operator* (const vector<T> &v1, const vector <T> &v2) {
	if (v1.size() == v2.size()) {
		vector<T> temp(v1.size());
		for (auto i = 0; i != v1.size(); ++i)
			temp[i] = v1[i]*v2[i];
		return temp;
	} else {
		cout << "Operation Vector*Vector called on vectors of different size. Exiting. . ." << endl;
		exit(-1);
	}
}

template <typename T>
vector<vector<T> > operator* (const vector<vector<T> > &m1, const vector<vector<T> > &m2) {
	if (m1[0].size() != m2.size()) {
		cout << "Operation Matrix*Matrix called on matrices of inconsistent sizes. Exiting. . ." << endl;
		exit(-1);
	} else {
		vector<vector<T> > temp(m1.size());
		for (auto i = 0; i != m1.size(); ++i)
			temp[i].resize(m2[0].size());
		for (auto i = 0; i != temp.size(); ++i)
			for (auto j = 0; j != temp[0].size(); ++j)
				for (auto k = 0; k != m2.size(); ++k)
					temp[i][j] = temp[i][j] + m1[i][k]*m2[k][j];
		return temp;
	}
}

template <typename T>
vector<T> operator/ (const vector<T> &v, const double &s) {
	vector<T> temp(v.size());
	for (auto i = 0; i != v.size(); ++i)
		temp[i] = v[i]/s;
	return temp;
}

template <typename T>
vector<vector<T> > operator/ (const vector<vector<T> > &m, const double &s) {
	vector<vector<T> > temp(m.size());
	for (auto &v : temp)
		v.resize(m[0].size());
	for (auto i = 0; i != temp.size(); ++i)
		for (auto j = 0; j != temp[i].size(); ++j)
			temp[i][j] = m[i][j]/s;
	return temp;
}

template <typename T>
vector<T> operator* (const vector<T> &v, const vector<vector<T> > &m) {
	if (v.size() != m.size()) {
		cout << "Operation Vector*Matrix called with inconsistent sizes of operands. Exiting. . ." << endl;
		exit(-1);
	} else {
		vector<T> temp(m[0].size());
		for (auto j = 0; j != temp.size(); ++j) {
			for (auto i = 0; i != v.size(); ++i) {
				temp[j] = temp[j] + v[i]*m[i][j];
			}
		}
		return temp;
	}
}

template <typename T>
vector<T> operator* (const vector<vector<T> > &m, const vector<T> &v) {
	if (m[0].size() != v.size()) {
		cout << "Operation Vector*Matrix called with inconsistent sizes of operands. Exiting. . ." << endl;
		exit(-1);
	} else {
		vector<T> temp(m.size());
		for (auto i = 0; i != temp.size(); ++i) {
			for (auto j = 0; j != m.size(); ++j) {
				temp[i] = temp[i] + m[i][j]*v[j];
			}
		}
		return temp;
	}
}

#endif