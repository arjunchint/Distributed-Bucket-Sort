#include <stdlib.h>
#include <string.h>
#include "proj2sorter.h"
#include "proj2sorter_impl.h"
#include <Random123/threefry.h>

 
static int ChoosePivot(Proj2Sorter sorter, MPI_Comm comm, size_t numKeysLocal, uint64_t *keys, uint64_t *pivot_p)
{
  int      err;
  uint64_t pivot = 0;

  if (numKeysLocal) {
    pivot = keys[numKeysLocal / 2];
  }
  err = MPI_Bcast(&pivot, 1, MPI_UINT64_T, 0, comm); MPI_CHK(err);
  *pivot_p = pivot;
  return 0;
}

/* instead of finding a perfect match, use this comparison operation to find
 * when a key fits between two entries in an array */
static int uint64_compare_pair(const void *key, const void *array)
{
  uint64_t Key = *((const uint64_t *) key);
  const uint64_t *pair = (const uint64_t *) array;
  if (Key < pair[0]) return -1;
  if (Key < pair[1]) return 0;
  return 1;
}


static int Proj2SorterSort_samplesort(Proj2Sorter sorter, MPI_Comm comm, size_t numKeysLocal, uint64_t *keys)
{

  uint64_t pivot = 0;
  uint64_t *lower_half, *upper_half;
  size_t   lower_size, upper_size;
  size_t   numKeysLocalNew = 0;
  uint64_t *keysNew = NULL;
  int      color;
  int      size, rank;
  int      equivRank;
  MPI_Comm subcomm;
  int      err;

  err = MPI_Comm_size(comm, &size); PROJ2CHK(err);
  err = MPI_Comm_rank(comm, &rank); PROJ2CHK(err);
  /* sort locally up front */
  qsort(keys, numKeysLocal, sizeof(*keys), uint64_comp_forward);
  //Proj2SorterSortLocal(sorter, numKeysLocal, keys, PROJ2SORT_FORWARD); PROJ2CHK(err);
  if (size == 1) {
    /* base case: nothing to do */
    return 0;
  }

  //Sample unifromly from local proc
  uint64_t v_local[size-1];
  if (numKeysLocal) {
  for(int i=1; i < size; i++){
    v_local[i-1] = keys[i*(numKeysLocal)/size];
  }
  }
  
  // Gather samples at root
  uint64_t *buffer=NULL;
  uint64_t *v_all=NULL;
  err = Proj2SorterGetWorkArray(sorter, size*(size-1), sizeof(uint64_t), &buffer); PROJ2CHK(err);
  err = Proj2SorterGetWorkArray(sorter, size-1, sizeof(uint64_t), &v_all); PROJ2CHK(err);
  MPI_Gather(v_local, size-1, MPI_UINT64_T, buffer, size-1, MPI_UINT64_T, 0, comm);
  
  // Sort and sample locally to find global bucket size
  if(rank == 0){
    err = Proj2SorterSortLocal(sorter, (size_t)(size*(size-1)), buffer, PROJ2SORT_FORWARD); PROJ2CHK(err);
    for(int i=1; i < size; i++){
      v_all[i-1] = buffer[i*(size-1)];
    }}
    err = MPI_Bcast(v_all, size-1, MPI_UINT64_T, 0, comm); MPI_CHK(err);  

////////////////// ALL PROCS KNOW BUCKETS TO SEGMENT BY /////////////////////////////////////////////
  // Resistrubute keys from processors based to new processors based on which bucket key is in
  int bucket_idx[size];
  int bucket_size[size]; 
  int idx=0;
  for(int pid=0; pid<size-1; pid++)
  {
    bucket_idx[pid] = idx;
    bucket_size[pid]=0;
    while((idx < numKeysLocal) && (keys[idx]<=v_all[pid]))
    {
      bucket_size[pid]++;
      idx++;
    }   
  }
  bucket_idx[size-1] = idx;
  bucket_size[size-1] = numKeysLocal - idx;

  //Create Recieve Buffer to store keys sent by other procs
  uint64_t *recvbuffer;
  err = Proj2SorterGetWorkArray(sorter, numKeysLocal*2, sizeof(uint64_t), &recvbuffer); PROJ2CHK(err);
  int recvBucketID[size];    
  int recvBucketLengths[size];           
        

  for(int ip=0; ip<size; ip++)
  { 
    MPI_Gather(&bucket_size[ip], 1, MPI_INT, recvBucketLengths,1,MPI_INT,ip, comm);
    if (rank == ip)
    {
      recvBucketID[0]=0;
      for(int i=1;i<size; i++)
      {
        recvBucketID[i] = recvBucketID[i-1]+recvBucketLengths[i-1];
      }
    }

    MPI_Gatherv(&keys[bucket_idx[ip]], bucket_size[ip],MPI_UINT64_T, recvbuffer,recvBucketLengths,recvBucketID,MPI_UINT64_T, ip, comm);
}

  /* Now the array is sorted, but we have to move it back to its original
   * distribution */
  // Attempted prefix sum to load balance across buckets but ran out of time
  {
    uint64_t myNewCount =0;
    for(int i = 0; i < size; i++){
      myNewCount += recvBucketLengths[i];
    }

    keysNew = recvbuffer;
    
    uint64_t myOldCount = numKeysLocal;
    qsort(recvbuffer, myNewCount, sizeof(*keys), uint64_comp_forward); 
    //err = Proj2SorterSortLocal(sorter, myNewCount, recvbuffer, PROJ2SORT_FORWARD); PROJ2CHK(err);

    uint64_t *oldOffsets;
    uint64_t *newOffsets;
    uint64_t oldOffset, newOffset;
    MPI_Request *recv_reqs, *send_reqs;
    int firstRecv = -1;
    int lastRecv = -2;
    int firstSend = -1;
    int lastSend = -2;
    int nRecvs, nSends;
    uint64_t thisOffset;

    err = Proj2SorterGetWorkArray(sorter, size + 1, sizeof(uint64_t), &oldOffsets); PROJ2CHK(err);
    err = Proj2SorterGetWorkArray(sorter, size + 1, sizeof(uint64_t), &newOffsets); PROJ2CHK(err);
    err = MPI_Allgather(&myOldCount,1,MPI_UINT64_T,oldOffsets,1,MPI_UINT64_T,comm); PROJ2CHK(err);
    err = MPI_Allgather(&myNewCount,1,MPI_UINT64_T,newOffsets,1,MPI_UINT64_T,comm); PROJ2CHK(err);

    oldOffset = 0;
    for (int i = 0; i < size; i++) {
      uint64_t count = oldOffsets[i];
      oldOffsets[i] = oldOffset;
      oldOffset += count;
    }
    oldOffsets[size] = oldOffset;

    newOffset = 0;
    for (int i = 0; i < size; i++) {
      uint64_t count = newOffsets[i];
      newOffsets[i] = newOffset;
      newOffset += count;
    }
    newOffsets[size] = newOffset;

    if (myOldCount) {
      for (int i = 0; i < size; i++) {
        if (newOffsets[i] <= oldOffsets[rank] && oldOffsets[rank] < newOffsets[i + 1]) {
          firstRecv = i;
        }
        if (newOffsets[i] <= oldOffsets[rank + 1] - 1 && oldOffsets[rank + 1] - 1 < newOffsets[i + 1]) {
          lastRecv = i + 1;
          break;
        }
      }
    }

    err = Proj2SorterGetWorkArray(sorter, lastRecv - firstRecv, sizeof(*recv_reqs), &recv_reqs); PROJ2CHK(err);

    thisOffset = oldOffsets[rank];
    nRecvs = 0;
    for (int i = firstRecv; i < lastRecv; i++) {
      size_t recvStart = thisOffset;
      size_t recvEnd   = newOffsets[i + 1];

      if (oldOffsets[rank + 1] < recvEnd) {
        recvEnd = oldOffsets[rank + 1];
      }
      if (recvEnd > recvStart) {
        err = MPI_Irecv(&keys[thisOffset - oldOffsets[rank]],(int) (recvEnd - recvStart), MPI_UINT64_T, i, PROJ2TAG_QUICKSORT, comm, &recv_reqs[nRecvs++]); MPI_CHK(err);
      }
      thisOffset = recvEnd;
    }

    if (myNewCount) {
      for (int i = 0; i < size; i++) {
        if (oldOffsets[i] <= newOffsets[rank] && newOffsets[rank] < oldOffsets[i + 1]) {
          firstSend = i;
        }
        if (oldOffsets[i] <= newOffsets[rank + 1] - 1 && newOffsets[rank + 1] - 1 < oldOffsets[i + 1]) {
          lastSend = i + 1;
          break;
        }
      }
    }

    err = Proj2SorterGetWorkArray(sorter, lastSend - firstSend, sizeof(*send_reqs), &send_reqs); PROJ2CHK(err);

    thisOffset = newOffsets[rank];
    nSends = 0;
    for (int i = firstSend; i < lastSend; i++) {
      size_t sendStart = thisOffset;
      size_t sendEnd   = oldOffsets[i + 1];

      if (newOffsets[rank + 1] < sendEnd) {
        sendEnd = newOffsets[rank + 1];
      }
      if (sendEnd > sendStart) {
        err = MPI_Isend(&keysNew[thisOffset - newOffsets[rank]],(int) (sendEnd - sendStart), MPI_UINT64_T, i, PROJ2TAG_QUICKSORT, comm, &send_reqs[nSends++]); MPI_CHK(err);
      }
      thisOffset = sendEnd;
    }

    err = MPI_Waitall(nRecvs, recv_reqs, MPI_STATUSES_IGNORE); PROJ2CHK(err);
    err = MPI_Waitall(nSends, send_reqs, MPI_STATUSES_IGNORE); PROJ2CHK(err);

    
    err = Proj2SorterRestoreWorkArray(sorter, numKeysLocal*4, sizeof(uint64_t), &recvbuffer); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, lastSend - firstSend, sizeof(*send_reqs), &send_reqs); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, lastRecv - firstRecv, sizeof(*recv_reqs), &recv_reqs); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, size + 1, sizeof(uint64_t), &newOffsets); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, size + 1, sizeof(uint64_t), &oldOffsets); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, size*(size-1), sizeof(uint64_t), &buffer); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, size-1, sizeof(uint64_t), &v_all); PROJ2CHK(err);  
}
  return 0;
}

static int Proj2SorterSort_quicksort_recursive(Proj2Sorter sorter, MPI_Comm comm, size_t numKeysLocal, uint64_t *keys)
{
  uint64_t pivot = 0;
  uint64_t *lower_half, *upper_half;
  size_t   lower_size, upper_size;
  size_t   numKeysLocalNew = 0;
  uint64_t *keysNew = NULL;
  int      color;
  int      size, rank;
  int      equivRank;
  MPI_Comm subcomm;
  int      err;

  err = MPI_Comm_size(comm, &size); PROJ2CHK(err);
  err = MPI_Comm_rank(comm, &rank); PROJ2CHK(err);
  /* sort locally up front */
  err = Proj2SorterSortLocal(sorter, numKeysLocal, keys, PROJ2SORT_FORWARD); PROJ2CHK(err);
  if (size == 1) {
    /* base case: nothing to do */
    return 0;
  }
  err = ChoosePivot(sorter, comm, numKeysLocal, keys, &pivot); PROJ2CHK(err);
  lower_half = NULL;
  upper_half = NULL;
  lower_size = 0;
  upper_size = 0;
  /* split the keys into those less than and those greater than the pivot */
  if (numKeysLocal) {
    if (pivot < keys[0]) { /* if the pivot is less than the first entry, they all are in the upper half */
      upper_half = keys;
      upper_size = numKeysLocal;
    } else if (pivot >= keys[numKeysLocal - 1]) { /* if the pivot is greater than or equal to the last entry, they all are in the lower half */
      lower_half = keys;
      lower_size = numKeysLocal;
    } else { /* otherwise use a binary search to split the keys */
      lower_half = keys;
      upper_half = (uint64_t *) bsearch(&pivot, keys, numKeysLocal - 1, sizeof(uint64_t), uint64_compare_pair);
      if (!upper_half) PROJ2ERR(MPI_COMM_SELF,5,"Could not find pivot\n");
      /* we have found the greatest entry less than the pivot: increase by one
       * to makek upper_half the start of the entries greater than or equal to
       * the pivot */
      upper_half++;
      lower_size = upper_half - lower_half;
      upper_size = numKeysLocal - lower_size;
    }
  }

  /* color the upper half to split the communicator */
  color = (rank >= (size / 2));
  equivRank = color ? (rank - (size / 2)) : (rank + (size / 2));
  if ((size % 2) == 0 || (rank != size - 1)) {
    /* Easy case: matched send and receive */
    MPI_Request sendreq;
    MPI_Status  recvstatus;

    if (equivRank > rank) {
      int numIncoming;

      err = MPI_Isend(upper_half, (int) upper_size, MPI_UINT64_T, equivRank, PROJ2TAG_QUICKSORT, comm, &sendreq); MPI_CHK(err);
      err = MPI_Probe(equivRank, PROJ2TAG_QUICKSORT, comm, &recvstatus); MPI_CHK(err);
      err = MPI_Get_count(&recvstatus, MPI_UINT64_T, &numIncoming); MPI_CHK(err);
      numKeysLocalNew = lower_size + numIncoming;
      err = Proj2SorterGetWorkArray(sorter, numKeysLocalNew, sizeof(uint64_t), &keysNew); PROJ2CHK(err);
      memcpy (keysNew, lower_half, lower_size * sizeof(*keysNew));
      err = MPI_Recv(&keysNew[lower_size], numIncoming, MPI_UINT64_T, equivRank, PROJ2TAG_QUICKSORT, comm, MPI_STATUS_IGNORE); PROJ2CHK(err);
      err = MPI_Wait(&sendreq, MPI_STATUS_IGNORE); MPI_CHK(err);
    } else {
      int numIncoming;

      err = MPI_Isend(lower_half, (int) lower_size, MPI_UINT64_T, equivRank, PROJ2TAG_QUICKSORT, comm, &sendreq); MPI_CHK(err);
      err = MPI_Probe(equivRank, PROJ2TAG_QUICKSORT, comm, &recvstatus); MPI_CHK(err);
      err = MPI_Get_count(&recvstatus, MPI_UINT64_T, &numIncoming); MPI_CHK(err);
      numKeysLocalNew = numIncoming + upper_size;
      err = Proj2SorterGetWorkArray(sorter, numKeysLocalNew, sizeof(uint64_t), &keysNew); PROJ2CHK(err);
      memcpy (&keysNew[numIncoming], upper_half, upper_size * sizeof(*keysNew));
      err = MPI_Recv(keysNew, numIncoming, MPI_UINT64_T, equivRank, PROJ2TAG_QUICKSORT, comm, MPI_STATUS_IGNORE); PROJ2CHK(err);
      err = MPI_Wait(&sendreq, MPI_STATUS_IGNORE); MPI_CHK(err);
    }
  }
  if (size % 2) { /* odd number of processes: size - 1 has been left out */
    MPI_Request sendreq;
    MPI_Status  recvstatus;

    /* (size - 1) sends its lower half to (size / 2) - 1 */
    if (rank == (size / 2) - 1) {
      int numIncoming;
      size_t numKeysLocalNewNew;
      uint64_t *keysNewNew;

      /* concatenate the addition from (size - 1) to the already constructed
       * keysNew */
      err = MPI_Probe(size - 1, PROJ2TAG_QUICKSORT, comm, &recvstatus); MPI_CHK(err);
      err = MPI_Get_count(&recvstatus, MPI_UINT64_T, &numIncoming); MPI_CHK(err);
      numKeysLocalNewNew = numKeysLocalNew + numIncoming;
      err = Proj2SorterGetWorkArray(sorter, numKeysLocalNewNew, sizeof(uint64_t), &keysNewNew); PROJ2CHK(err);
      memcpy (keysNewNew, keysNew, numKeysLocalNew * sizeof(*keysNewNew));
      err = Proj2SorterRestoreWorkArray(sorter, numKeysLocalNew, sizeof(uint64_t), &keysNew); PROJ2CHK(err);
      err = MPI_Recv(&keysNewNew[numKeysLocalNew], numIncoming, MPI_UINT64_T, size - 1, PROJ2TAG_QUICKSORT, comm, MPI_STATUS_IGNORE); PROJ2CHK(err);
      numKeysLocalNew = numKeysLocalNewNew;
      keysNew = keysNewNew;
    } else if (rank == (size - 1)) {
      /* (size - 1) does not receive from anyone */
      err = MPI_Isend(lower_half, (int) lower_size, MPI_UINT64_T, (size / 2) - 1, PROJ2TAG_QUICKSORT, comm, &sendreq); MPI_CHK(err);
      numKeysLocalNew = upper_size;
      err = Proj2SorterGetWorkArray(sorter, numKeysLocalNew, sizeof(uint64_t), &keysNew); PROJ2CHK(err);
      memcpy (keysNew, upper_half, upper_size * sizeof(*keysNew));
      err = MPI_Wait(&sendreq, MPI_STATUS_IGNORE); MPI_CHK(err);
    }
  }
  err = MPI_Comm_split(comm, color, rank, &subcomm); PROJ2CHK(err);
  err = Proj2SorterSort_quicksort_recursive(sorter, subcomm, numKeysLocalNew, keysNew); PROJ2CHK(err);
  err = MPI_Comm_free(&subcomm); PROJ2CHK(err);

  /* Now the array is sorted, but we have to move it back to its original
   * distribution */
  {
    uint64_t myOldCount = numKeysLocal;
    uint64_t myNewCount = numKeysLocalNew;
    uint64_t *oldOffsets;
    uint64_t *newOffsets;
    uint64_t oldOffset, newOffset;
    MPI_Request *recv_reqs, *send_reqs;
    int firstRecv = -1;
    int lastRecv = -2;
    int firstSend = -1;
    int lastSend = -2;
    int nRecvs, nSends;
    uint64_t thisOffset;

    err = Proj2SorterGetWorkArray(sorter, size + 1, sizeof(uint64_t), &oldOffsets); PROJ2CHK(err);
    err = Proj2SorterGetWorkArray(sorter, size + 1, sizeof(uint64_t), &newOffsets); PROJ2CHK(err);
    err = MPI_Allgather(&myOldCount,1,MPI_UINT64_T,oldOffsets,1,MPI_UINT64_T,comm); PROJ2CHK(err);
    err = MPI_Allgather(&myNewCount,1,MPI_UINT64_T,newOffsets,1,MPI_UINT64_T,comm); PROJ2CHK(err);

    oldOffset = 0;
    for (int i = 0; i < size; i++) {
      uint64_t count = oldOffsets[i];
      oldOffsets[i] = oldOffset;
      oldOffset += count;
    }
    oldOffsets[size] = oldOffset;

    newOffset = 0;
    for (int i = 0; i < size; i++) {
      uint64_t count = newOffsets[i];
      newOffsets[i] = newOffset;
      newOffset += count;
    }
    newOffsets[size] = newOffset;

    if (myOldCount) {
      for (int i = 0; i < size; i++) {
        if (newOffsets[i] <= oldOffsets[rank] && oldOffsets[rank] < newOffsets[i + 1]) {
          firstRecv = i;
        }
        if (newOffsets[i] <= oldOffsets[rank + 1] - 1 && oldOffsets[rank + 1] - 1 < newOffsets[i + 1]) {
          lastRecv = i + 1;
          break;
        }
      }
    }

    err = Proj2SorterGetWorkArray(sorter, lastRecv - firstRecv, sizeof(*recv_reqs), &recv_reqs); PROJ2CHK(err);

    thisOffset = oldOffsets[rank];
    nRecvs = 0;
    for (int i = firstRecv; i < lastRecv; i++) {
      size_t recvStart = thisOffset;
      size_t recvEnd   = newOffsets[i + 1];

      if (oldOffsets[rank + 1] < recvEnd) {
        recvEnd = oldOffsets[rank + 1];
      }
      if (recvEnd > recvStart) {
        err = MPI_Irecv(&keys[thisOffset - oldOffsets[rank]],(int) (recvEnd - recvStart), MPI_UINT64_T, i, PROJ2TAG_QUICKSORT, comm, &recv_reqs[nRecvs++]); MPI_CHK(err);
      }
      thisOffset = recvEnd;
    }

    if (myNewCount) {
      for (int i = 0; i < size; i++) {
        if (oldOffsets[i] <= newOffsets[rank] && newOffsets[rank] < oldOffsets[i + 1]) {
          firstSend = i;
        }
        if (oldOffsets[i] <= newOffsets[rank + 1] - 1 && newOffsets[rank + 1] - 1 < oldOffsets[i + 1]) {
          lastSend = i + 1;
          break;
        }
      }
    }

    err = Proj2SorterGetWorkArray(sorter, lastSend - firstSend, sizeof(*send_reqs), &send_reqs); PROJ2CHK(err);

    thisOffset = newOffsets[rank];
    nSends = 0;
    for (int i = firstSend; i < lastSend; i++) {
      size_t sendStart = thisOffset;
      size_t sendEnd   = oldOffsets[i + 1];

      if (newOffsets[rank + 1] < sendEnd) {
        sendEnd = newOffsets[rank + 1];
      }
      if (sendEnd > sendStart) {
        err = MPI_Isend(&keysNew[thisOffset - newOffsets[rank]],(int) (sendEnd - sendStart), MPI_UINT64_T, i, PROJ2TAG_QUICKSORT, comm, &send_reqs[nSends++]); MPI_CHK(err);
      }
      thisOffset = sendEnd;
    }

    err = MPI_Waitall(nRecvs, recv_reqs, MPI_STATUSES_IGNORE); PROJ2CHK(err);
    err = MPI_Waitall(nSends, send_reqs, MPI_STATUSES_IGNORE); PROJ2CHK(err);

    err = Proj2SorterRestoreWorkArray(sorter, lastSend - firstSend, sizeof(*send_reqs), &send_reqs); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, lastRecv - firstRecv, sizeof(*recv_reqs), &recv_reqs); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, size + 1, sizeof(uint64_t), &newOffsets); PROJ2CHK(err);
    err = Proj2SorterRestoreWorkArray(sorter, size + 1, sizeof(uint64_t), &oldOffsets); PROJ2CHK(err);
  }
  err = Proj2SorterRestoreWorkArray(sorter, numKeysLocalNew, sizeof(uint64_t), &keysNew); PROJ2CHK(err);
  return 0;
}

int Proj2SorterSort_quicksort(Proj2Sorter sorter, size_t numKeysLocal, uint64_t *keys)
{
  int      err;
  int      size, rank;
  err = MPI_Comm_size(sorter->comm, &size); PROJ2CHK(err);
  /* initiate recursive call */
  if(numKeysLocal > 20000000){
  err = Proj2SorterSort_samplesort(sorter, sorter->comm, numKeysLocal, keys); PROJ2CHK(err);
  //PROJ2LOG(sorter->comm, "USING SAMPLE SORT");
  } else {
  err = Proj2SorterSort_quicksort_recursive(sorter, sorter->comm, numKeysLocal, keys); PROJ2CHK(err);
 // PROJ2LOG(sorter->comm, "USING QUICKSORT SORT");
  }
  return 0;
}
