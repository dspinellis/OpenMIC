#ifndef EXTENSIBLEPARTITION_H
#define EXTENSIBLEPARTITION_H

#include <set>
#include <vector>
#include <cassert>
#include <iterator>

#include "Partition.h"
#include "entropy.h"
#include "debug.h"

// An output iterator used for counting generated elements
class CounterOutputIterator :
	public iterator<output_iterator_tag, void, void, void, void>
{
public:
	explicit CounterOutputIterator(int &c) : counter(c) {}

	CounterOutputIterator & operator*() { return *this; }
	template <typename T> void operator=(T const& rhs) { }
	CounterOutputIterator & operator=(const CounterOutputIterator &rhs) {
		counter = rhs.counter;
		return *this;
	}
	CounterOutputIterator operator++() {
		CounterOutputIterator new_obj(*this);
		counter++;
		return new_obj;
	}
	CounterOutputIterator & operator++(int) { counter++; return *this; }
private:
	int &counter;
};

/*
 * A partial partition over the X axis clumps whose endpoint can be extended.
 * The clump numbers specified are 1-based.
 */
class ExtensiblePartition {
public:
	typedef vector <int> Points;		// Points storage type
private:
	const Partition &clumps;	// All possible partitions on the horizontal axis
	const Partition &q;		// Partition along the vertical axis
	double hq;			// Entropy of q
	Points points;			// The ordinals of clumps over which we partition
	int last_column_points;		// Number of points in last column
	vector <int> last_column_partitioned_points;	// Number of points in last column partitioned by rows
public:
	// Create a new adjustable partition of size 2
	ExtensiblePartition(const Partition &c, const Partition &qin, double hqin, int partition, int end) :
		clumps(c), q(qin), hq(hqin), points(3) {
			assert(partition > 0);
			assert(partition < clumps.size());
			assert(end >= 0);
			assert(end <= clumps.size());
			assert(partition <= end);
			points[0] = 0;
			points[1] = partition;
			points[2] = end;
	}

	// Return a new ExtensiblePartition with an added point at its end
	ExtensiblePartition add_point(int p) const {
		assert(p >= points.back());
		ExtensiblePartition r(*this);
		if (p == points.back())
			return (r);
		r.points.push_back(p);
		return (r);
	}

	// Return the number of points in the specified horizontal axis partition ending in p
	inline int number_of_horizontal_partition_points(int p) {
		assert(p > 0);
		assert(p <= points.size());
		int n = 0;

		Partition::const_iterator start(clumps.begin());
		advance(start, points[p - 1]);

		Partition::const_iterator end(clumps.begin());
		advance(end, points[p]);

		for (Partition::const_iterator i = start; i != end; i++)
			n += i->size();
		if (DP())
			cout << var(p) << var(n) << endl;
		return n;
	}

	// Return the points in the specified horizontal axis partition p
	set <const Point *> horizontal_partition_points(int p) {
		assert(p > 0);
		assert(p <= points.size());

		set <const Point *> result;

		Partition::const_iterator start(clumps.begin());
		advance(start, points[p - 1]);

		Partition::const_iterator end(clumps.begin());
		advance(end, points[p]);

		/*
		 * Form the union of all sets.
		 * This is efficient, because
		 * insert range is linear in N if the range is already sorted by value_comp().
		 */
		for (Partition::const_iterator i = start; i != end; i++)
			result.insert(i->begin(), i->end());
		return result;
	}

	// Return H(p)
	double hp() {
		// Create probability vector for H(p)
		// TODO: Cache the following two as members and adjust them in add_point
		vector <double> pp;
		int npoints = 0;
		for (int i = 1; i < points.size(); i++) {
			int n = number_of_horizontal_partition_points(i);
			pp.push_back(n);
			npoints += n;
		}
		// Convert cardinalities to probability weights
		for (vector <double>::iterator i = pp.begin(); i != pp.end(); i++)
			*i /= npoints;
		return H(pp);
	}

	// Return H(p, q)
	double hpq() {
		// Create probability vector for H(p, q)
		// TODO: Cache the following two as members and adjust them in add_point
		vector <double> ppq;
		int npoints = 0;
		for (Partition::const_iterator i = q.begin(); i != q.end(); i++) {
			for (int j = 1; j < points.size(); j++) {
				set <const Point *> hpoints(horizontal_partition_points(j));
				int n = 0;
				CounterOutputIterator count_elements(n);
				set_intersection(hpoints.begin(), hpoints.end(), i->begin(), i->end(), count_elements);
				ppq.push_back(n);
				npoints += n;
				if (DP()) {
					cout << "x=" << j << " y=" << distance(q.begin(), i);
					cout << "\nPoints along the horizontal axis:\n";
					copy(i->begin(), i->end(), ostream_iterator<const Point *>(cout, " "));
					cout << "\nPoints along the vertical axis:\n";
					copy(hpoints.begin(), hpoints.end(), ostream_iterator<const Point *>(cout, " "));
					cout << "Intersection:\n";
					set_intersection(hpoints.begin(), hpoints.end(), i->begin(), i->end(), ostream_iterator<const Point *>(cout, " "));
				}
			}
		}
		if (DP()) {
			copy(ppq.begin(), ppq.end(), ostream_iterator<double>(cout, "\t"));
			cout << endl << var(npoints) << endl;
		}
		// Convert cardinalities to probability weights
		for (vector <double>::iterator i = ppq.begin(); i != ppq.end(); i++)
			*i /= npoints;
		return H(ppq);
	}
};

#endif
