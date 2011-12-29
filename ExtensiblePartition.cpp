#include "ExtensiblePartition.h"

static Partition dummy;

const Partition *ExtensiblePartition::clumps;	// All possible partitions on the horizontal axis
const Partition *ExtensiblePartition::q;	// Partition along the vertical axis
