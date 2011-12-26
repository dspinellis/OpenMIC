#ifndef ADJUSTABLEPARTITION_H
#define ADJUSTABLEPARTITION_H

#include <set>
#include <vector>
#include <cassert>

#include "Partition.h"
#include "entropy.h"
#include "debug.h"

/*
 * A partition whose points can be moved and added among its
 * defining clumps.
 */
class AdjustablePartition {
public:
	typedef set <Partition::const_iterator> Points;		// Points storage type
private:
	const Partition &clumps;	// All possible horizontal partitions
	const Partition &q;		// Vertical partition
	double hq;			// Entropy of q
	Points points;			// The ordinals of clumps over which we partition
	vector <int> cached_partition_points;
	static const int CACHE_EMPTY = -1;
	// Return the orginal of a points iterator
	inline int ordinal(Points::value_type i) const { return i - clumps.begin(); }
public:
	// Create a new adjustable partition of size 2
	AdjustablePartition(const Partition &c, const Partition &qin, double hqin, Points::value_type start, Points::value_type end) :
		clumps(c), q(qin), hq(hqin), cached_partition_points(2, CACHE_EMPTY) {
			assert(start > clumps.begin());
			assert(start <= clumps.end());
			assert(end > clumps.begin());
			assert(end <= clumps.end());
			points.insert(start);
			points.insert(end);
	}

	// Return a new AdjustablePartition with an added point
	AdjustablePartition add_point(Points::value_type p) const {
		assert(p > clumps.begin());
		assert(p <= clumps.end());
		AdjustablePartition r(*this);
		if (points.find(p) == points.end())
			return (r);
		r.points.insert(p);
		// Invalidate caches
		r.cached_partition_points.insert(r.cached_partition_points.begin() + ordinal[p], CACHE_EMPTY);
		if (cached_partition_points.begin() + ordinal(p) + 1 != cached_partition_points.end())
			r.cached_partition_points[ordinal(p) + 1] = CACHE_EMPTY;
		return (r);
	}

	// Return the number of points in the specified horizontal partition p
	inline int partition_points(Points::value_type p) {
		if (cached_partition_points[ordinal(p)] != CACHE_EMPTY)
			return cached_partition_points[ordinal(p)];
		int n = 0;
		Points::value_type start(points.find(p));
		if (start != points.begin())
			;//start --;
		for (Points::value_type i = start; i != p; i++)
			n += i->size();
		cached_partition_points[ordinal(p)] = n;
		cout << var(p) << var(n) << endl;
		return n;
	}

	// Return H(p)
	double hp() {
		// Create probability vector for H(p)
		vector <double> pp;
		int npoints = 0;
		for (int i = 0; i < points.size(); i++) {
			int n = partition_points(i);
			pp.push_back(n);
			npoints += n;
		}
		// Convert cardinalities to probability weights
		for (vector <double>::iterator i = pp.begin(); i != pp.end(); i++)
			*i /= npoints;
		return H(pp);
	}
};

#endif /* ADJUSTBLEPARTITION_H */
