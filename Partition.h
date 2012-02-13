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
