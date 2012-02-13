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

#include "Partition.h"

// Output partition p on o
ostream&
operator<<(ostream& o, const Partition &p)
{
	for (Partition::const_iterator i = p.begin(); i != p.end(); i++) {
		for (Partition::value_type::const_iterator j = i->begin(); j != i->end(); j++)
			o << *j;
		o << '\n';
	}
	o << endl;
}
