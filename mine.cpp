#include <vector>
#include <set>
#include <cassert>
#include <algorithm>	// sort
#include <iostream>	// cout
#include <fstream>	// ifstream
#include <cstdio>	// perror
#include <cstdlib>	// exit
#include <cmath>	// log2
#include <iterator>	// ostream_iterator

#define DP() 0
#define var(x) " " #x "=" << x

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

struct Point {
	double x, y;
	Point() {}
	Point(double px, double py) : x(px), y(py) {}
	friend ostream& operator<<(ostream& o, const Point &p);
	friend ostream& operator<<(ostream& o, const Point *p);
};

// Functor for sorting points on x
struct less_x : public binary_function<const Point *, const Point *, bool> {
	bool operator()(const Point *a, const Point *b) { return a->x < b->x; }
};

// Functor for sorting points on y
struct less_y : public binary_function<const Point *, const Point *, bool> {
	bool operator()(const Point *a, const Point *b) { return a->y < b->y; }
};

ostream&
operator<<(ostream& o, const Point &p)
{
	o << p.x << ',' << p.y << '\n';
	return o;
}

ostream&
operator<<(ostream& o, const Point *p)
{
	o << *p;
	return o;
}

typedef vector<set<const Point *> > Partition;


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

void
show_partition(const Partition &p)
{
	for (Partition::const_iterator i = p.begin(); i != p.end(); i++) {
		for (Partition::value_type::const_iterator j = i->begin(); j != i->end(); j++)
			cout << *j;
		cout << '\n';
	}
	cout << endl;
}

template <typename T>
void
show_matrix(const vector <vector <T> > &m)
{
	for (typename vector <vector <T> >::const_iterator i = m.begin(); i != m.end(); i++)
		show_vector(*i);
}

// Return the mutual information of data in the specified grid
// See http://www.scholarpedia.org/article/Mutual_information
double
mutual_information(vector <Point> &data, vector <double> &xbins, vector <double> ybins)
{
	// Probability mass * n (SOM 2.1)
	vector <vector <int> > pm(ybins.size() + 1, vector<int>(xbins.size() + 1, 0));
	// Marginals * n; see http://en.wikipedia.org/wiki/Marginal_probability
	vector <int> xmar(xbins.size() + 1), ymar(ybins.size() + 1);
	for (vector <Point>::const_iterator i = data.begin(); i != data.end(); i++) {
		int xpos = lower_bound(xbins.begin(), xbins.end(), i->x) - xbins.begin();
		int ypos = lower_bound(ybins.begin(), ybins.end(), i->y) - ybins.begin();
		pm[ypos][xpos]++;
		xmar[xpos]++;
		ymar[ypos]++;
	}
	cout << "Probability mass\n";
	show_matrix(pm);
	cout << "X marginals\n";
	show_vector(xmar);
	cout << "Y marginals\n";
	show_vector(ymar);

	double n = data.size();
	double result = 0;
	for (int y = 0; y < ybins.size() + 1; y++)
		for (int x = 0; x < xbins.size() + 1; x++)
			if (pm[y][x] > 0 && xmar[x] > 0 && ymar[y] > 0) {		// Avoid log(0)
				result += (double)pm[y][x] * (double)pm[y][x] / n * log2((double)pm[y][x] / ((double)xmar[x] * (double)ymar[y] / n));
				// cout << pm[y][x] << '\t' << xmar[x] << '\t' << ymar[y] << '\t' << result << endl;
			}
	cout << "Mutual information: " << result << endl;
	return result;
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
			cout << "same_points=" << same_points << " currently_assigned=" << currently_assigned << " have=" << abs(currently_assigned + same_points - desired_row_size) << " want=" << abs(currently_assigned - desired_row_size) << endl;
		if (currently_assigned == 0 ||
		    // Distance from target to handle tie breaks
		    abs(currently_assigned + same_points - desired_row_size) <= abs(currently_assigned - desired_row_size)) {
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
	Partition::value_type const *current_partition = NULL;
	for (int i = 0; i < data.size(); i++) {
		if (DP())
			cout << "Look at point " << i << ": " << *data[i] << endl;
		// Indirect through data to get correct point ordinals
		if (ypartition_map[data[i] - &*points.begin()] != current_partition) {
			clumps.push_back(Partition::value_type());	// Start a new partition
			current_partition = ypartition_map[data[i] - &*points.begin()];
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

	assert(0);
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
vector <vector <double> >
characteristic_matrix(vector <Point> &data, double b, int clump_factor)
{
	assert(clump_factor > 0);
	assert(b > 3);

	// data2 (D\bot) is (y1, x1), (y2, x2) ...
	vector <Point> data2;
	transform(data.begin(), data.end(), back_inserter(data2), flip());

	// Calculare the information content matrix (lines 2-6)
	vector <vector <double> > mi(2, vector<double>(b / 2, 0));
	vector <vector <double> > mi2(2, vector<double>(b / 2, 0));
	for (int y = 2; y <= b / 2; y++) {
		int x = b / y;
		cout << "x=" << x << " y=" << y << " b=" << b << endl;
		vector <double> mmi(max_mi(data, x, y, clump_factor * x));
		cout << "max_mi" << endl;
		show_vector(mmi);
		mi.push_back(mmi);
		//mi.push_back(max_mi(data, x, y, clump_factor * x));
		mi2.push_back(max_mi(data2, x, y, clump_factor * x));
	}

	// Fill-in the characteristic matrix (lines 7-10)
	vector <vector <double> > cm(b / 2, vector<double>(b / 2, 0));
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
	exit(0);
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

	vector <double> xbins, ybins;
	read_vector(argv[2], xbins);
	read_vector(argv[3], ybins);

	mutual_information(points, xbins, ybins);
	vector <vector <double> > cm(characteristic_matrix(points, pow(points.size(), grid_exponent), clumping));
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
			show_partition(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}

	{	// 3 elements into 3 rows
		vector <Point> test(p, p + 3);
		Partition got(equipartition_y_axis(test, 3));
		Partition expect(point_to_ptr(test, {0, 1, 2}));
		if (DP()) {
			show_vector(test);
			show_partition(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 6 elements into 3 rows
		vector <Point> test(p, p + 6);
		Partition got(equipartition_y_axis(test, 3));
		Partition expect(point_to_ptr(test, {0, 0, 1, 1, 2, 2, }));
		if (DP()) {
			show_vector(test);
			show_partition(expect);
			show_partition(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 3 elements into 2 rows
		vector <Point> test(p, p + 3);
		Partition got(equipartition_y_axis(test, 2));
		Partition expect(point_to_ptr(test, {0, 1, 1}));
		if (DP()) {
			show_vector(test);
			show_partition(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 8 elements into 3 rows
		vector <Point> test(p, p + 8);
		Partition got(equipartition_y_axis(test, 3));
		Partition expect(point_to_ptr(test, {0, 0, 1, 1, 1, 2, 2, 2, }));
		if (DP()) {
			show_vector(test);
			show_partition(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 9 elements into 3 rows with tie
		vector <Point> test(p, p + 9);
		Partition got(equipartition_y_axis(test, 3));
		Partition expect(point_to_ptr(test, {0, 0, 0, 1, 1, 1, 1, 2, 2, }));
		if (DP()) {
			show_vector(test);
			show_partition(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 10 elements into 5 rows with two ties
		vector <Point> test(p, p + 10);
		Partition got(equipartition_y_axis(test, 5));
		Partition expect(point_to_ptr(test, {0, 0, 1, 1, 2, 2, 2, 3, 3, 3, }));
		if (DP()) {
			show_vector(test);
			show_partition(got);
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
			show_partition(expect);
			show_partition(got);
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
			show_partition(expect_y);
			cout << "Obtained Y equipartition" << endl;
			show_partition(got_y);

			cout << "Expected clumps" << endl;
			show_partition(expect_clumps);
			cout << "Obtained clumps" << endl;
			show_partition(got_clumps);
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
#endif
