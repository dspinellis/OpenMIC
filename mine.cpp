#include <vector>
#include <cassert>
#include <algorithm>	// sort
#include <iostream>	// cout
#include <fstream>	// ifstream
#include <cstdio>	// perror
#include <cstdlib>	// exit
#include <cmath>	// log2
#include <iterator>	// ostream_iterator

#define DP() 0

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
};

// Functor for sorting points on x
struct less_x : public binary_function<const Point &, const Point &, bool> {
	bool operator()(const Point &a, const Point &b) { return a.x < b.x; }
};

// Functor for sorting points on y
struct less_y : public binary_function<const Point &, const Point &, bool> {
	bool operator()(const Point &a, const Point &b) { return a.y < b.y; }
};

ostream&
operator<<(ostream& o, const Point &p)
{
	o << p.x << ',' << p.y << '\n';
	return o;
}

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
vector <int>
equipartition_y_axis(const vector <Point> &data, int y)
{
	int n = data.size();
	int i = 0;			// Input position in data
	int desired_row_size = n / y;
	vector <int> q(n, -1);
	int current_row = 0;		// Output position in q
	int currently_assigned = 0;	// Equivalent of #
	if (DP())
		cout << "\nn=" << n << endl;
	do {
		if (DP())
			cout << "current_row=" << current_row << " i=" << i << " currently_assigned=" << currently_assigned << " desired_row_size=" << desired_row_size << endl;
		// Line 6: Cardinality of S is all that is needed; exploit ordering by y
		int same_points = 1;	// Number of points with same y (|S|)
		for (int j = i + 1; j < n && data[j].y == data[i].y; j++)
			same_points++;
		if (DP())
			cout << "same_points=" << same_points << " currently_assigned=" << currently_assigned << " have=" << abs(currently_assigned + same_points - desired_row_size) << " want=" << abs(currently_assigned - desired_row_size) << endl;
		if (currently_assigned == 0 ||
		    // Distance from target to handle tie breaks
		    abs(currently_assigned + same_points - desired_row_size) <= abs(currently_assigned - desired_row_size)) {
			for (int j = 0; j < same_points; j++) {
				q[i + j] = current_row;
				if (DP())
					cout << "Assign point " << i + j << " to row " << current_row << endl;
			}
		    	i += same_points;
			currently_assigned += same_points;
			if (DP())
				cout << "i=" << i << " currently_assigned=" << currently_assigned << " current_row=" << current_row << endl;
			if (y - current_row)
				desired_row_size = (n - i + currently_assigned) / (y - current_row);
			else
				desired_row_size = numeric_limits<int>::max();
		} else {
			current_row++;
			currently_assigned = 0;
		}
	} while(i < n);
	return q;
}

/*
 * Algorithm 2
 * "Returns a list of scores (I_2 ... I_x) such that each I_l is the maximum value of I(P;Q) over all
 * partitions P of size l."
 */
vector <double>
optimize_x_axis(const vector <Point> &data, const vector <int> &q, int x, int clumps)
{
	assert(0);
}

/*
 * Algorithm 4
 * "Returns a set of mutual information scores (0, 0, I_{2,y} ... I_{x,y}) such that I_{i,j} is
 * heuristically close to the highest achievable mutual information score using i rows and j columns"
 * We copy-in data, because we will sort it later on.
 */
vector <double>
max_mi(vector <Point> data, int x, int y, int clumps)
{
	assert(x > 1);
	assert(y > 1);
	assert(clumps > 1);

	sort(data.begin(), data.end(), less_y());	// Do it here to avoid duplication in algorithm 3 line 1
	vector <int> q(equipartition_y_axis(data, y));
	return optimize_x_axis(data, q, x, clumps);
}

// Functor for flipping x, y
struct flip : public unary_function<const Point &, Point> {
	Point operator()(const Point &p) { return Point(p.y, p.x); }
};

/*
 * Algorithm 5
 * Return \forall where x * y < b the matrix with the highest information content
 */
vector <vector <double> >
characteristic_matrix(vector <Point> &data, double b, int clumps)
{
	assert(clumps > 0);
	assert(b > 3);

	// data2 (D\bot) is (y1, x1), (y2, x2) ...
	vector <Point> data2;
	transform(data.begin(), data.end(), back_inserter(data2), flip());
	// Sort by X value
	sort(data.begin(), data.end(), less_x());
	sort(data2.begin(), data2.end(), less_x());

	// Calculare the information content matrix (lines 2-6)
	vector <vector <double> > mi(2, vector<double>(b / 2, 0));
	vector <vector <double> > mi2(2, vector<double>(b / 2, 0));
	for (int y = 2; y <= b / 2; y++) {
		int x = b / y;
		cout << "x=" << x << " y=" << y << " b=" << b << endl;
		vector <double> mmi(max_mi(data, x, y, clumps * x));
		cout << "max_mi" << endl;
		show_vector(mmi);
		mi.push_back(mmi);
		//mi.push_back(max_mi(data, x, y, clumps * x));
		mi2.push_back(max_mi(data2, x, y, clumps * x));
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

int
main(int argc, char *argv[])
{
	vector <Point> points;
	double grid_exponent = 0.6;
	int clumping = 15;

	test_equipartition();
	exit(0);
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
void
test_equipartition()
{
	//           1       2       3       4       5       6       7       8       9       10      11      12      13        14
	Point p[] = {{1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}, {6, 6}, {7, 7}, {7, 7}, {7, 7}, {8, 8}, {9, 9}, {10, 10}, {11, 11}};

	{	// 2 elements into 2 rows
		vector <Point> test(p, p + 2);
		vector <int> got(equipartition_y_axis(test, 2));
		vector <int> expect = {0, 1};
		if (DP()) {
			show_vector(test);
			show_vector(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}

	{	// 3 elements into 3 rows
		vector <Point> test(p, p + 3);
		vector <int> got(equipartition_y_axis(test, 3));
		vector <int> expect = {0, 1, 2};
		if (DP()) {
			show_vector(test);
			show_vector(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 6 elements into 3 rows
		vector <Point> test(p, p + 6);
		vector <int> got(equipartition_y_axis(test, 3));
		vector <int> expect = {0, 0, 1, 1, 2, 2, };
		if (DP()) {
			show_vector(test);
			show_vector(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 3 elements into 2 rows
		vector <Point> test(p, p + 3);
		vector <int> got(equipartition_y_axis(test, 2));
		vector <int> expect = {0, 1, 1};
		if (DP()) {
			show_vector(test);
			show_vector(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 8 elements into 3 rows
		vector <Point> test(p, p + 8);
		vector <int> got(equipartition_y_axis(test, 3));
		vector <int> expect = {0, 0, 1, 1, 1, 2, 2, 2, };
		if (DP()) {
			show_vector(test);
			show_vector(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 9 elements into 3 rows with tie
		vector <Point> test(p, p + 9);
		vector <int> got(equipartition_y_axis(test, 3));
		vector <int> expect = {0, 0, 0, 1, 1, 1, 1, 2, 2, };
		if (DP()) {
			show_vector(test);
			show_vector(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
	{	// 10 elements into 5 rows with two ties
		vector <Point> test(p, p + 10);
		vector <int> got(equipartition_y_axis(test, 5));
		vector <int> expect = {0, 0, 1, 1, 2, 2, 2, 3, 3, 3, };
		if (DP()) {
			show_vector(test);
			show_vector(got);
		}
		assert(equal(expect.begin(), expect.end(), got.begin()));
	}
}
#endif
