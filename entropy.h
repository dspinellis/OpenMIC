#ifndef ENTROPY_H
#define ENTROPY_H

#include <vector>
#include <cmath>	// log2

#include "Partition.h"

// Return the Shannon entropy of a probability vector P
// See http://www.scholarpedia.org/article/Entropy#Shannon_entropy
template <typename T>
double
H(const T &p)
{
	double sum = 0;

	for (typename T::const_iterator i = p.begin(); i != p.end(); i++)
		sum += *i * log2(*i);
	return -sum;
}

// Return the Shannon entropy of the specified partition of a set of npoints
double
H(const Partition &part)
{
	vector <double> p;
	int npoints = 0;
	for (Partition::const_iterator i = part.begin(); i != part.end(); i++) {
		p.push_back(i->size());
		npoints += i->size();
	}
	// Convert cardinalities to probability weights
	for (vector <double>::iterator i = p.begin(); i != p.end(); i++)
		*i /= npoints;
	return H(p);
}

#endif /* ENTROPY_H */
