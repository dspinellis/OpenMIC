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
