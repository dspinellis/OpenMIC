#include <vector>
#include <cassert>
#include <cmath>	// log2

#include "debug.h"
#include "Partition.h"
#include "entropy.h"

static const double EPSILON = 1e-8;

// Return the Shannon entropy of a probability vector P
// See http://www.scholarpedia.org/article/Entropy#Shannon_entropy
template <typename T>
double
H(const T &p)
{
	double sum = 0;

	for (typename T::const_iterator i = p.begin(); i != p.end(); i++)
		if (*i != 0)		// XXX To handle emtpy grid areas. Not specified in the paper
			sum += *i * log2(*i);
	/*
	 * Entropy is equal to -sum
    	 * (c) The entropy of a partition is nonnegative and equal to zero if and only if one of the
	 * elements Ai of the partition has measure 1 (and all other elements have measure zero).
	 * (d) The entropy of a partition into n sets is highest for the measure which assigns equal
	 * values 1n to these sets. The entropy then equals log2n .
	 */
	if (DP())
		cout << "H=" << -sum << endl;
	assert(-sum >= 0);
	assert(-sum <= log2(p.size()) + EPSILON);
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
