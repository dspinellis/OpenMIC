#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
#include <set>

#include "Point.h"

using namespace std;

typedef vector<set<const Point *> > Partition;

// Output partition p on o
ostream& operator<<(ostream& o, const Partition &p);

#endif /* PARTITION_H */
