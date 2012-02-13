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
