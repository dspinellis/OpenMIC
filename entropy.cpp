/*-
 * Copyright 2011-2012 Diomidis Spinellis
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

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
		if (*i != 0)
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

template double H(const vector <double> &p);
