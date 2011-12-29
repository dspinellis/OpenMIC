#include "point.h"

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
