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

#ifndef ENTROPY_H
#define ENTROPY_H

#include "Partition.h"

// Return the Shannon entropy of a probability vector P
// See http://www.scholarpedia.org/article/Entropy#Shannon_entropy
template <typename T> double H(const T &p);

// Return the Shannon entropy of the specified partition of a set of npoints
double H(const Partition &part);

#endif /* ENTROPY_H */
