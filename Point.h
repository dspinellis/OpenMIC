#ifndef POINT_H
#define POINT_H

#include <iostream>	// ostream

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

#endif /* POINT_H */
