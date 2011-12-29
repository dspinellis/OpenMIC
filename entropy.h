#ifndef ENTROPY_H
#define ENTROPY_H

#include "Partition.h"

// Return the Shannon entropy of a probability vector P
// See http://www.scholarpedia.org/article/Entropy#Shannon_entropy
template <typename T> double H(const T &p);

// Return the Shannon entropy of the specified partition of a set of npoints
double H(const Partition &part);

#endif /* ENTROPY_H */
