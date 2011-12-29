#include <vector>
#include <set>
#include <limits>
#include <cassert>
#include <algorithm>	// sort
#include <iostream>	// cout
#include <fstream>	// ifstream
#include <cstdio>	// perror
#include <cstdlib>	// exit
#include <cmath>	// log2
#include <iterator>	// ostream_iterator

#include "Point.h"
#include "Partition.h"
#include "ExtensiblePartition.h"

#include "entropy.h"
#include "debug.h"

typedef vector< vector<double> > matrix;

/*
 * Missing from STL.  See Effective STL item 36 and
 * http://stackoverflow.com/questions/1448817/why-there-is-no-stdcopy-if-algorithm
 */
template <class InputIterator, class OutputIterator, class Predicate>
OutputIterator copy_if(InputIterator begin, InputIterator end,
                       OutputIterator result, Predicate pred) {
	while (begin != end) {
		if (pred(*begin))
			*result++ = *begin;
		    ++begin;
	}
	return result;
}

using namespace std;

// Read a vector from the specified file
void
read_vector(const char *name, vector <double> &v)
{
	ifstream vfile(name);
	if (!vfile.is_open()) {
		perror(name);
		exit(1);
	}
	for (;;) {
		double d;
		vfile >> d;
		if (!vfile.good())
			break;
		v.push_back(d);
	}
	cout << name << endl;
	for (vector <double>::const_iterator i = v.begin(); i != v.end(); i++)
		cout << *i << endl;
}

template <typename T>
void
show_vector(const vector <T> &v)
{
	copy(v.begin(), v.end(), ostream_iterator<T>(cout, "\t"));
	cout << endl;
}

template <typename T>
void
show_matrix(const vector <vector <T> > &m)
{
	for (typename vector <vector <T> >::const_iterator i = m.begin(); i != m.end(); i++)
		show_vector(*i);
}

/*
 * Algorithm 3
 * "Returns a map Q: D -> {1, ...., y}  such that Q(p) is the row assignment of the point p and there
 * is approximately the same number of points in each row"
 */
Partition
equipartition_y_axis(const vector <Point> &points, int y)
{
	assert(y > 1);

	// Create vector of pointers to points sorted by y
	vector <const Point *> data(points.size());
	for (int i = 0; i < points.size(); i++)
		data[i] = &(points[i]);
	sort(data.begin(), data.end(), less_y());

	int npoints = data.size();
	int i = 0;			// Input position in data
	int desired_row_size = npoints / y;
	Partition q(1);
	int current_row = 0;		// Output position in q
	int currently_assigned = 0;	// Equivalent of #
	if (DP())
		cout << "\nn=" << npoints << endl;
	do {
		if (DP())
			cout << "current_row=" << current_row << " i=" << i << " currently_assigned=" << currently_assigned << " desired_row_size=" << desired_row_size << endl;
		// Line 6: Cardinality of S is all that is needed; exploit ordering by y
		int same_points = 1;	// Number of points with same y (|S|)
		for (int j = i + 1; j < npoints && data[j]->y == data[i]->y; j++)
			same_points++;
		if (DP())
			cout << var(currently_assigned) << var(same_points) << var(desired_row_size) <<
			    " LHS=" << abs(currently_assigned + same_points - desired_row_size) <<
			    " RHS=" << abs(currently_assigned - desired_row_size) << endl;
		if (currently_assigned == 0 ||
		    // Distance from target to handle tie breaks
		    abs(currently_assigned + same_points - desired_row_size) <= abs(currently_assigned - desired_row_size)) {
			if (DP())
				cout << "Assign points to row " << current_row << endl;
			for (int j = 0; j < same_points; j++) {
				q[current_row].insert(data[i + j]);
				if (DP())
					cout << "Assign point " << i + j << " to row " << current_row << endl;
			}
		    	i += same_points;
			currently_assigned += same_points;
			if (DP())
				cout << "i=" << i << " currently_assigned=" << currently_assigned << " current_row=" << current_row << endl;
			if (y - current_row)
				desired_row_size = (npoints - i + currently_assigned) / (y - current_row);
			else
				desired_row_size = numeric_limits<int>::max();
		} else {
			current_row++;
			if (DP())
				cout << "Advance current_row to " << current_row << endl;
			q.push_back(Partition::value_type());
			currently_assigned = 0;
		}
	} while(i < npoints);
	return q;
}

/*
 * Partition data by "drawing x-axis partition lines only between runs of consecutive points that fall
 * in the same row of the y-axis partition Q."
 * "Return the minimal partition that separates every pair of points that lie in distinct clumps."
 * Not listed in pseudocode.
 */
Partition
get_clumps_partition(const vector <Point> &points, const Partition &q)
{
	// Create vector of pointers to points sorted by x
	vector <const Point *> data(points.size());
	for (int i = 0; i < points.size(); i++)
		data[i] = &(points[i]);
	sort(data.begin(), data.end(), less_x());

	// Create a map from a point ordinal to its y partition
	vector <const Partition::value_type *> ypartition_map(points.size());
	for (Partition::const_iterator i = q.begin(); i != q.end(); i++)
		for (Partition::value_type::const_iterator j = i->begin(); j != i->end(); j++)
			ypartition_map[*j - &*(points.begin())] = &*i;

	Partition clumps;
	Partition::value_type const *current_y_partition = NULL;
	for (int i = 0; i < data.size(); i++) {
		if (DP())
			cout << "Look at point " << i << ": " << *data[i] << endl;
		// Indirect through data to get correct point ordinals
		if (ypartition_map[data[i] - &*points.begin()] != current_y_partition) {
			clumps.push_back(Partition::value_type());	// Start a new partition
			current_y_partition = ypartition_map[data[i] - &*points.begin()];
		}
		clumps.back().insert(data[i]);
	}
	return clumps;
}

/*
 * Repartition by merging true clumps into superclumps in a way that aims to have each
 * superclump contain approximately the same number of points returning at most max_clumps
 * clumps.
 * npoints is the total number of points.
 * Not listed in pseudocode.
 */
Partition
get_superclumps_partition(const Partition &input_partitions, int npoints, int max_clumps)
{
	assert(!input_partitions.empty());
	assert(npoints > 0);
	assert(max_clumps > 0);

	if (input_partitions.size() <= max_clumps)
		return (input_partitions);

	int points_per_partition = npoints / max_clumps;
	Partition q(1);
	int output_partition = 0;	// Output position in q
	int currently_assigned = 0;	// Points assigned in this iteration
	int total_assigned = 0;		// Points assigned over all iterations
	Partition::const_iterator i = input_partitions.begin();
	do {
		if (DP()) {
			cout << var(abs(currently_assigned + (int)i->size() - points_per_partition)) << endl;
			cout << var(abs(currently_assigned - points_per_partition)) << endl;
		}
		if (currently_assigned == 0 ||
		    // Distance from target to handle tie breaks
		    abs(currently_assigned + (int)i->size() - points_per_partition) <= abs(currently_assigned - points_per_partition)) {
			q[output_partition].insert(i->begin(), i->end());
			currently_assigned += i->size();
			total_assigned += i->size();
		    	i++;
			if (max_clumps - output_partition)
				points_per_partition = (npoints - total_assigned + currently_assigned) / (max_clumps - output_partition);
			else
				points_per_partition = numeric_limits<int>::max();
			if (DP()) {
				cout << var(npoints) << var(total_assigned) << var(currently_assigned) << endl;
				cout << var(output_partition) << var(points_per_partition) << endl;
			}
		} else {
			output_partition++;
			q.push_back(Partition::value_type());
			currently_assigned = 0;
		}
	} while(i != input_partitions.end());
	return q;
}


// Return the point ordinals corresponding to each clump
vector <int>
get_clump_point_ordinals(const Partition &clumps)
{
	vector <int> result;
	result.reserve(clumps.size());
	result.push_back(0);
	for (Partition::const_iterator i = clumps.begin(); i != clumps.end(); i++)
		result.push_back(result.back() + i->size());
	assert(result.size() == clumps.size() + 1);
	return result;
}

/*
 * Algorithm 2
 * "Returns a list of scores (I_2 ... I_x) such that each I_l is the maximum value of I(P;Q) over all
 * partitions P of size l."
 * Max_clumps (\^k in the text) is the maximum number of clumps to use.
 */
vector <double>
optimize_x_axis(const vector <Point> &points, const Partition &q, int x, int max_clumps)
{
	assert(x > 1);

	Partition clumps(get_superclumps_partition(get_clumps_partition(points, q), points.size(), max_clumps));

	vector <int> c(get_clump_point_ordinals(clumps));

	int k = clumps.size();		// Compared to Algorithm 2 this is k + 1

	matrix I(k, vector <double> (x + 1));
	vector < vector <ExtensiblePartition> > P(k, vector <ExtensiblePartition> (x + 1));

	double hq = H(q);

	ExtensiblePartition::set_clumps(&clumps);
	ExtensiblePartition::set_q(&q);

	// Find the optimal partitions of size 2 for elements of 2 to k clumps
	for (int t = 2; t < k; t++) {
		// Find the best partition point s
		int maxs = 0;
		double maxh = -numeric_limits<double>::max();
		vector <ExtensiblePartition> cand(t + 1);		// Candidate partitions
		for (int s = 1; s <= t; s++) {
			cand[s] = ExtensiblePartition(s, t);
			double hdiff = cand[s].hp() - cand[s].hpq();
			if (hdiff > maxh) {
				maxs = s;
				maxh = hdiff;
			}
			if (DP())
				cout << var(t) << var(s) << var(hdiff) << var(maxh) << var(maxs) << endl;
		}
		assert(maxs != 0);
		P[t][2] = cand[maxs];
		I[t][2] = hq + maxh;
	}

	// Inductively build the rest of the table of optimal partitions

	// Build up for larger and larger partitions
	for (int l = 3; l <= x; l++)
		// Try various clump points
		for (int t = 2; t < k; t++) {
			int maxs = 0;
			double maxf = -numeric_limits<double>::max();
			vector <ExtensiblePartition> cand(t + 1);		// Candidate partitions
			// Find the best partition to use from the previously found partitions
			for (int s = 2; s <= t; s++) {
				cand[s] = P[s][l - 1].add_point(t);
				if (cand[s].number_of_columns() < l)
					continue;
				double sum = 0;
				double column_points = cand[s].number_of_horizontal_partition_points(l);
				if (column_points == 0)
					continue;
				for (int i = 1; i <= q.size(); i++) {
					double cell_points = cand[s].number_of_cell_points(i, l);
					if (cell_points == 0)
						continue;
					sum += cell_points / c[t] * log2(cell_points / column_points);
				}
				double f = (double)c[s] / (double) c[t] * (I[s][l - 1] - hq) + sum;
				if (DP())
					cout << var(l) << var(t) << var(s) << var(f) << endl;
				if (f > maxf) {
					maxs = s;
					maxf = f;
				}
			}
			if (maxs == 0) {
				P[t][l] = P[t][l - 1];
				I[t][l] = I[t][l - 1];
			} else {
				P[t][l] = cand[maxs];
				I[t][l] = hq + P[t][l].hp() - P[t][l].hpq();
			}
		}
	return I[k - 1];
}

/*
 * Algorithm 4
 * "Returns a set of mutual information scores (0, 0, I_{2,y} ... I_{x,y}) such that I_{i,j} is
 * heuristically close to the highest achievable mutual information score using i rows and j columns."
 * Max_clumps (\^k in the text) is the maximum number of clumps to use.
 */
vector <double>
max_mi(vector <Point> &data, int x, int y, int max_clumps)
{
	assert(x > 1);
	assert(y > 1);
	assert(max_clumps > 1);

	Partition q(equipartition_y_axis(data, y));
	return optimize_x_axis(data, q, x, max_clumps);
}

// Functor for flipping x, y
struct flip : public unary_function<const Point &, Point> {
	Point operator()(const Point &p) { return Point(p.y, p.x); }
};

/*
 * Algorithm 5
 * Return \forall where x * y < b the matrix with the highest information content
 * The clump factor (c in the text) "determines by what factor clumps may outnumber columns
 * when OptimizeXAxis is called. When trying to partition the x-axis into x columns, the
 * algorithm will start with at most cx clumps."
 */
matrix
characteristic_matrix(vector <Point> &data, int b, int clump_factor)
{
	assert(clump_factor > 0);
	assert(b > 3);

	// data2 (D\bot) is (y1, x1), (y2, x2) ...
	vector <Point> data2;
	transform(data.begin(), data.end(), back_inserter(data2), flip());

	// Calculate the information content matrix (lines 2-6)
	matrix mi(2, vector<double>(b / 2, 0));
	matrix mi2(2, vector<double>(b / 2, 0));
	for (int y = 2; y <= b / 2; y++) {
		int x = b / y;
		if (1 || DP()) {
			cout << "x=" << x << " y=" << y << " b=" << b << endl;
			vector <double> mmi(max_mi(data, x, y, clump_factor * x));
			cout << "max_mi" << endl;
			show_vector(mmi);
			mi.push_back(mmi);
		} else
			mi.push_back(max_mi(data, x, y, clump_factor * x));
		mi2.push_back(max_mi(data2, x, y, clump_factor * x));
	}

	// Fill-in the characteristic matrix (lines 7-10)
	matrix cm(b / 2 + 1, vector<double>(b / 2 + 1, 0));
	for (int x = 2; x <= b / 2; x++)
		for (int y = 2; y <= b / 2; y++) {
			if (x * y > b)
				continue;
			double ixy = max(mi[y][x], mi2[x][y]);
			cm[y][x] = ixy / min(log2(x), log2(y));
		}
	return cm;
}

void test_equipartition();
void test_get_clumps_partition();
void test_get_superclumps_partition();
void test_H();
void test_ExtensiblePartition();
void test_CounterOutputIterator();
void test_get_clump_point_ordinals();

// Return the maximal information coefficient
double
mic(const matrix &cm)
{
	double result = -numeric_limits<double>::max();

	for (matrix::const_iterator i = cm.begin(); i != cm.end(); i++)
		for (matrix::value_type::const_iterator j = i->begin(); j != i->end(); j++)
			if (*j > result)
				result = *j;
	return result;
}

// Return the maximum asymmetry score
double
mas(const matrix &cm)
{
	double result = -numeric_limits<double>::max();

	for (int i = 0; i < cm.size(); i++)
		for (int j = 0; j < cm[i].size(); j++) {
			double v = fabs(cm[i][j] - cm[j][i]);
			if (v > result)
				result = v;
		}
	return result;
}

// Return the maximum edge value
double
mev(const matrix &cm)
{
	double result = -numeric_limits<double>::max();

	for (int i = 0; i < cm.size(); i++)
		if (cm[i][2] > result)			// Or maybe 1 XXX?
			result = cm[i][2];
	for (int i = 0; i < cm[2].size(); i++)
		if (cm[2][i] > result)
			result = cm[2][i];		// Or maybe 1 XXX?
	return result;
}

// Return the complexity
double
mcn(const matrix &cm, double mic, double epsilon)
{
	double result = numeric_limits<double>::max();

	for (int i = 0; i < cm.size(); i++)
		for (int j = 0; j < cm[i].size(); j++) {
			if (cm[i][j] < (1 - epsilon) * mic)
				continue;
			double v = log2((i + 1) * (j + 1));
			if (v < result)
				result = v;
		}
	return result;
}


int
main(int argc, char *argv[])
{
	vector <Point> points;
	double grid_exponent = 0.6;
	int clumping = 15;

#ifdef TEST
	test_equipartition();
	test_get_clumps_partition();
	test_get_superclumps_partition();
	test_H();
	test_ExtensiblePartition();
	test_CounterOutputIterator();
	test_get_clump_point_ordinals();
	cout << "All tests finished" << endl;
#endif

	// Read space-separated points
	ifstream pfile(argv[1]);
	if (!pfile.is_open()) {
		perror(argv[1]);
		exit(1);
	}
	for (;;) {
		Point p;
		pfile >> p.x >> p.y;
		if (!pfile.good())
			break;
		points.push_back(p);
	}
	// Print the points read
	cout << argv[1] << endl;
	for (vector <Point>::const_iterator i = points.begin(); i != points.end(); i++)
		cout << i->x << ' ' << i->y << endl;

	int b = pow(points.size(), grid_exponent);
	if (b < 4) {
		cerr << "not enough points" << endl;
		exit(1);
	}

	matrix cm(characteristic_matrix(points, b, clumping));

	cout << "Characteristic matrix:" << endl;
	show_matrix(cm);
	cout << endl;

	// Report results
	cout << "X var,Y var,MIC (strength),MAS (non-monotonicity),"
		"MEV (functionality),MCN (complexity)" << endl;
	double m;
	cout << "x, y, " <<
		(m = mic(cm)) << ',' <<
		mas(cm) << ',' <<
		mev(cm) << ',' <<
		mcn(cm, m, 0) << endl;
	return 0;
}

#ifdef TEST
/*
 * Return a partition of points as indicated by their corresponding ordinals. 
 */
static Partition
point_to_ptr(const vector <Point> &points, const vector <int> ordinals)
{
	// Calculate number of partitions
	set <int> ordinal_set(ordinals.begin(), ordinals.end());
	Partition result(ordinal_set.size());

	int prev = -1;
	int n = 0;
	for (vector <int>::const_iterator i = ordinals.begin(); i != ordinals.end(); i++)
		result[*i].insert(&points[n++]);
	return result;
}

void
test_equipartition()
{
	//           1       2       3       4       5       6       7       8       9       10      11      12      13        14
	Point p[] = {{1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}, {6, 6}, {7, 7}, {7, 7}, {7, 7}, {8, 8}, {9, 9}, {10, 10}, {11, 11}};

	{	// 2 elements into 2 rows
		vector <Point> test(p, p + 2);
		Partition got(equipartition_y_axis(test, 2));
		Partition expect(point_to_ptr(test, {0, 1}));
		if (DP()) {
			show_vector(test);
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}

	{	// 3 elements into 3 rows
		vector <Point> test(p, p + 3);
		Partition got(equipartition_y_axis(test, 3));
		Partition expect(point_to_ptr(test, {0, 1, 2}));
		if (DP()) {
			show_vector(test);
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 6 elements into 3 rows
		vector <Point> test(p, p + 6);
		Partition got(equipartition_y_axis(test, 3));
		Partition expect(point_to_ptr(test, {0, 0, 1, 1, 2, 2, }));
		if (DP()) {
			show_vector(test);
			cout << expect;
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 3 elements into 2 rows
		vector <Point> test(p, p + 3);
		Partition got(equipartition_y_axis(test, 2));
		Partition expect(point_to_ptr(test, {0, 1, 1}));
		if (DP()) {
			show_vector(test);
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 8 elements into 3 rows
		vector <Point> test(p, p + 8);
		Partition got(equipartition_y_axis(test, 3));
		Partition expect(point_to_ptr(test, {0, 0, 1, 1, 1, 2, 2, 2, }));
		if (DP()) {
			show_vector(test);
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 9 elements into 3 rows with tie
		vector <Point> test(p, p + 9);
		Partition got(equipartition_y_axis(test, 3));
		Partition expect(point_to_ptr(test, {0, 0, 0, 1, 1, 1, 1, 2, 2, }));
		if (DP()) {
			show_vector(test);
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 10 elements into 5 rows with two ties
		vector <Point> test(p, p + 10);
		Partition got(equipartition_y_axis(test, 5));
		Partition expect(point_to_ptr(test, {0, 0, 1, 1, 2, 2, 2, 3, 3, 3, }));
		if (DP()) {
			show_vector(test);
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 2 elements into 2 rows unsorted
		vector <Point> test(p, p + 2);
		swap(*test.begin(), *(test.begin() + 1));
		Partition expect(point_to_ptr(test, {1, 0}));
		Partition got(equipartition_y_axis(test, 2));
		if (DP()) {
			show_vector(test);
			cout << expect;
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 22 elements into 2 rows with 20 ties
		vector <Point> test(22, Point(10,10));
		test[0] = Point(1,1);
		test[1] = Point(2,2);
		Partition got(equipartition_y_axis(test, 2));
		vector <int> expect_ordinals(22, 1);
		expect_ordinals[0] = expect_ordinals[1] = 0;
		Partition expect(point_to_ptr(test, expect_ordinals));
		if (DP()) {
			show_vector(test);
			cout << got;
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
}

void
test_get_clumps_partition()
{
	/*
	 * 3         x
	 * 2       x
	 * 1   x x
	 * 0 x         x
	 *   0 1 2 3 4 5
	 *
	 * Consider the above points.
	 * Their Y axis equipartition would be {{(0,0), (5,0)}, {(1, 1), (2, 1)}, {(3,2), (4,3)}}
	 * Partition ordinals:                    0      0        1       1         2      2
	 * The corresponding clumps would be {{(0,0)},  {(1, 1), (2, 1)}, {(3,2), (4, 3)}, {(5,0)}}
	 * Partition ordinals:                    0      1        1         2      2         3
	 */

	Point p[] = {{0, 0}, {1, 1}, {3, 2}, {2, 1}, {5, 0}, {4, 3}};

	{	// Nonconsecutive and consecutive points
		vector <Point> test(p, p + 6);

		// Six points into three bins
		Partition got_y(equipartition_y_axis(test, 3));
		Partition expect_y(point_to_ptr(test, {0, 1, 2, 1, 0, 2}));

		Partition got_clumps(get_clumps_partition(test, got_y));
		Partition expect_clumps(point_to_ptr(test, {0, 1, 2, 1, 3, 2}));
		if (DP()) {
			cout << "Vector" << endl;
			show_vector(test);
			cout << "Expected Y equipartition" << endl;
			cout << expect_y;
			cout << "Obtained Y equipartition" << endl;
			cout << got_y;

			cout << "Expected clumps" << endl;
			cout << expect_clumps;
			cout << "Obtained clumps" << endl;
			cout << got_clumps;
		}
		assert(equal(expect_y.begin(), expect_y.end(), got_y.begin()));
		assert(equal(expect_clumps.begin(), expect_clumps.end(), got_clumps.begin()));
	}
}

void
test_get_superclumps_partition()
{
	const int MANY_PARTITIONS = 1000;
	const int FEW_PARTITIONS = 100;
	const int MAX_POINTS_PER_PARTITION = 50;
	int total_points = 0;

	// Create many partitions with a random number of points each up to MAX_POINTS_PER_PARTITION
	Partition many(MANY_PARTITIONS);
	srand(42);		// Ensure deterministic behavior
	// Fill in partitions with a random amount of points
	for (Partition::iterator i = many.begin(); i != many.end(); i++) {
		int npoints = rand() % MAX_POINTS_PER_PARTITION;
		for (int j = 0; j < npoints; j++)
			i->insert(new Point(rand(), rand()));
		total_points += npoints;
	}

	Partition few(get_superclumps_partition(many, total_points, FEW_PARTITIONS));

	if (DP())
		cout << "few.size()=" << few.size() << endl;
	assert(few.size() == FEW_PARTITIONS);
	int points_in_few = 0;
	for (Partition::const_iterator i = few.begin(); i != few.end(); i++) {
		assert(i->size() < total_points / FEW_PARTITIONS + MAX_POINTS_PER_PARTITION - 1);
		points_in_few += i->size();
		for (Partition::value_type::const_iterator j = i->begin(); j != i->end(); j++)
			delete *j;
	}
	assert(points_in_few == total_points);
}

void
test_H()
{
	// Should be an exact result!
	assert(H(vector <double>({1./8, 1./4, 1./8, 1./2})) == 7./4);

	// Above example on partitions
	vector <Point> test({{1, 1}, {1, 1}, {1, 1}, {1, 1}, {2, 2,}, {2, 2,}, {3, 3}, {4, 4}});
	Partition got(equipartition_y_axis(test, 4));
	Partition expect(point_to_ptr(test, {0, 0, 0, 0, 1, 1, 2, 3}));
	assert(equal(expect.begin(), expect.end(), got.begin()));
	assert(H(got) == 7./4);
}

void
test_ExtensiblePartition()
{
	/*
	 * 4             x
	 * 3         x
	 * 2       x
	 * 1   x x
	 * 0 x         x
	 *   0 1 2 3 4 5 6
	 *
	 * Consider the above points.
	 * Their Y axis equipartition will be {{(0,0), (5,0)}, {(1, 1), (2, 1)}, {(3,2), (4,3), (6, 4)}}
	 * Partition ordinals:                    0      0        1       1         2      2     2
	 * The corresponding clumps will be {{(0,0)},  {(1, 1), (2, 1)}, {(3,2), (4, 3)}, {(5,0)}, {(6,4)}}
	 *                                             1                 2                3        4
	 */

	Point p[] = {{0, 0}, {1, 1}, {3, 2}, {2, 1}, {5, 0}, {4, 3}, {6, 4}};

	vector <Point> test(p, p + 7);

	// Six points into three bins
	Partition q(equipartition_y_axis(test, 3));
	Partition clumps(get_clumps_partition(test, q));

	ExtensiblePartition::set_q(&q);
	ExtensiblePartition::set_clumps(&clumps);

	// Test ctors
	ExtensiblePartition a12(1, 2);
	assert(a12.number_of_horizontal_partition_points(1) == 1);
	assert(a12.number_of_horizontal_partition_points(2) == 2);

	ExtensiblePartition a13(1, 3);
	assert(a13.number_of_horizontal_partition_points(1) == 1);
	assert(a13.number_of_horizontal_partition_points(2) == 4);

	ExtensiblePartition a23(2, 3);
	assert(a23.number_of_horizontal_partition_points(1) == 3);
	assert(a23.number_of_horizontal_partition_points(2) == 2);

	ExtensiblePartition a24(2, 4);
	assert(a24.number_of_horizontal_partition_points(1) == 3);
	assert(a24.number_of_horizontal_partition_points(2) == 3);

	ExtensiblePartition a25(2, 5);
	assert(a25.number_of_horizontal_partition_points(1) == 3);
	assert(a25.number_of_horizontal_partition_points(2) == 4);

	// Test add_point
	ExtensiblePartition a234(a23.add_point(4));
	assert(a234.number_of_horizontal_partition_points(1) == 3);
	assert(a234.number_of_horizontal_partition_points(2) == 2);
	assert(a234.number_of_horizontal_partition_points(3) == 1);

	ExtensiblePartition a124(a12.add_point(4));
	assert(a124.number_of_horizontal_partition_points(1) == 1);
	assert(a124.number_of_horizontal_partition_points(2) == 2);
	assert(a124.number_of_horizontal_partition_points(3) == 3);

	// Verify entropy of the partition across the horizontal axis
	assert(a124.hp() == H(vector <double>({1./6, 2./6, 3./6})));

	/*
	 * 4  |   |     |x
	 * 3  |   |  x  |
	 * 2  |   |x    |
         *----+---+-----+-
	 * 1  |x x|     |
         *----+---+-----+-
	 * 0 x|   |    x|
	 *   0|1 2|3 4 5 6
	 *    |   |     |
	 */

	// Verify number_of_cell_points
	assert(a124.number_of_cell_points(1, 1) == 1);
	assert(a124.number_of_cell_points(1, 2) == 0);
	assert(a124.number_of_cell_points(1, 3) == 1);

	assert(a124.number_of_cell_points(2, 1) == 0);
	assert(a124.number_of_cell_points(2, 2) == 2);
	assert(a124.number_of_cell_points(2, 3) == 0);

	assert(a124.number_of_cell_points(3, 1) == 0);
	assert(a124.number_of_cell_points(3, 2) == 0);
	assert(a124.number_of_cell_points(3, 3) == 2);

	// Verify entropy of the points across both partitions
	assert(fabs(a124.hpq() - H(vector <double>({
		0,	0,	2./6,
		0,	2./6,	0./6,
		1./6,	0,	1./6,
	}))) < 1e-10);

	// Test add_point of previously added point
	ExtensiblePartition a1244(a124.add_point(4));
	assert(a1244.number_of_horizontal_partition_points(1) == 1);
	assert(a1244.number_of_horizontal_partition_points(2) == 2);
	assert(a1244.number_of_horizontal_partition_points(3) == 3);

	ExtensiblePartition a122(a12.add_point(2));
	assert(a122.number_of_horizontal_partition_points(1) == 1);
	assert(a122.number_of_horizontal_partition_points(2) == 2);

}

void
test_get_clump_point_ordinals()
{
	/*
	 * 3         x
	 * 2       x
	 * 1   x x
	 * 0 x         x
	 *   0 1 2 3 4 5
	 *
	 * Consider the above points.
	 * Their Y axis equipartition would be {{(0,0), (5,0)}, {(1, 1), (2, 1)}, {(3,2), (4,3)}}
	 * Partition ordinals:                    0      0        1       1         2      2
	 * The corresponding clumps would be {{(0,0)},  {(1, 1), (2, 1)}, {(3,2), (4, 3)}, {(5,0)}}
	 * Partition ordinals:                    0      1        1         2      2         3
	 * Point ordinals                         0      1                  3                5	   6
	 */

	Point p[] = {{0, 0}, {1, 1}, {3, 2}, {2, 1}, {5, 0}, {4, 3}};

	// Six points into three bins
	vector <Point> test(p, p + 6);
	Partition ypartition(equipartition_y_axis(test, 3));
	Partition clumps(get_clumps_partition(test, ypartition));
	vector <int> ordinals(get_clump_point_ordinals(clumps));
	if (DP()) {
		cout << clumps;
		show_vector(ordinals);
	}
	// According to Yakir we must get
	assert(ordinals[0] == 0);
	assert(ordinals[1] == 1);
	assert(ordinals[2] == 3);
	assert(ordinals[3] == 5);
	assert(ordinals[4] == 6);
}

void
test_CounterOutputIterator()
{
	vector<int> v(5, 0);
	int n = 0;
	CounterOutputIterator count_elements(n);
	copy(v.begin(), v.end(), count_elements);
	assert(n == 5);
}
#endif
