#if !defined(PROJ2SORTER_IMPL_H)
#define      PROJ2SORTER_IMPL_H

#include "proj2sorter.h"

typedef struct _proj2memlink *Proj2MemLink;

struct _proj2memlink
{
  Proj2MemLink next;
  size_t size;
  char *array;
};

struct _proj2sorter
{
  MPI_Comm     comm;
  Proj2MemLink avail;
  Proj2MemLink inUse;
};

/* A parallel sort implementation based on divide and conquer, by choosing a
 * pivot and dividing into subproblems.  In parallel, each processor has a
 * communicating partner on the other half of the communicator.  Once the
 * pivot is known, each partner sends its entries less than the pivot to the
 * lower of the pair and the others to the greater of the pair. */
int Proj2SorterSort_quicksort(Proj2Sorter sorter, size_t numKeysLocal, uint64_t *keys);

/* A divide and conquer based sorting algorith that that samples the distribution of keys amonst
 * the processors to redistrubute them based on buckets */
//int Proj2SorterSort_samplesort(Proj2Sorter sorter, MPI_Comm comm, size_t numKeysLocal, uint64_t *keys);


#endif
