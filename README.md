# Project 2: Sorting

**Due ~~Monday, October 23~~ Friday, October 27**

**Your repo transfer request** (name format `cse6230fa17-proj2-gtusername` or `cse6230fa17-proj2-gtusername1-gtusername2` must be submitted by **11:59:55 p.m. on that day**

- If working in pairs, maintain and submit one repo.
- If repo history is entirely from one half of the pair, points may be
  deducted, unless other evidence can be provided of the participation of both
  parties.
- David (Hansol)'s repo is https://bitbucket.org/dhsuh/cse6230fa17-proj2-hsuh7
- As Arjun and David decided to work as a group at a later time, and decided to use
  Arjun's repo as main repo, at that point David's repo was ignored. In order to 
  prove cooperation, David's repo is hereby sited. As merging two repos with
  many similar files proved to be difficult, merging was not done.

## Objectives

- The goal of this project is to *use profiling* to optimize the performance of
  an MPI-based library for distributed memory sorting.
    * The main library interface is (as declared in [proj2sorter.h](proj2sorter.h)):

    /* This is the default implementation of sorting:
     * \param[in] sorter       The sorting context.  Put all of your customizations
     *                         in this object.  Defined in proj2sorter_impl.h, where
     *                         you can change the struct to include more data
     * \param[in] numKeysLocal The number of keys on this process.
     * \param[in/out] keys     The input array.  On output, should be globally
     *                         sorted in ascending order.
     * \return                 Non-zero if an error occured.
     */
    int Proj2SorterSort(Proj2Sorter sorter, size_t numKeysLocal, uint64_t *keys);

    * The small library comes with some logging and error macros (see
      [proj2.h](proj2.h)) as well as an interface for obtaining/restoring workspace
      arrays (see `Proj2SorterGetWorkArray()` and `Proj2SorterRestoreWorkArray()` in
      [proj2sorter.h](proj2sorter.h).  To be memory neutral, restore every
      workspace that you get.

    * Functioning parallel implementations have been provided, one based on
      [quicksort](https://en.wikipedia.org/wiki/Quicksort#Parallelization) and
      one base on [bitonic mergesort](https://en.wikipedia.org/wiki/Bitonic_sorter).

    * A [template library](https://github.com/swenson/sort) originated by Chris
      Swenson has been imported for a quicksort implementation that is faster
      than `qsort` from the standard library.  The template library includes
      other implementations that you are welcome to explore.

    * Indeed, as with previous assignments, the implementation details are up
      to you.  As with Checkpoint 1 of the final, there is a test program
      ([test_proj2.c](test_proj2.c)), which may not be edited, that calls your
      library.  It will test the sorting bandwidth (bytes sorted per second) of
      your code on uniformly generated random data at varying numbers of *keys
      per MPI process* (a *key* in our library is just a `uint64_t`: a large
      integer).  You may incorporate additional files into your library by
      adding a `Makefile.inc` file to your project.

    * The final arbiter of your performance is the submission script for
      Stampede2,
      [slurm-proj2-stampede2-final.sh](slurm-proj2-stampede2-final.sh), which
      may not be edited.  Notice some aspects of the script:

        * `-N 4`: the job runs on 4 nodes, so pure shared memory is impossible.
        * `git rev-parse HEAD`: tells us which version of your code was submitted
          to the job scheduler.
        * `git diff-files`: tells us if you have any unsaved changes in your code.
        * The code tests your performance at 64 tasks per node, at
          granularities of 80, 5120, 327,680, and 20,971,520 keys per MPI
          processes. This results in global sizes of 20,480, 1,310,720,
          83,886,080, and 5,368,709,120.  **Careful**: this last number is too
          big to be represented by `int`, watch your conversions.
        * The [harmonic
          mean](https://en.wikipedia.org/wiki/Harmonic_mean) of the
          sorting bandwidths for these four problems is our metric of performance.
          A harmonic mean emphasizes small outliers more than an arithmetic mean.

    * The final submission script currently runs in 5 minutes on Stampede.
      With four nodes, this equates to 20 SU minutes.  Each *individual* should
      try to use at most 10 SU hours on this project (pairs thus have twice the
      SUs available).  You should have your own submission scripts for testing
      and development:

        * The test program is run like

    ./test_proj2 MIN_KEYS_PER_PROCESS MAX_KEYS_PER_PROCESS MULTIPLIER SEED NUM_REPS

This means that the test program seeds the random number generator with `SEED`,
starts with `MIN_KEYS_PER_PROCESS`, tests `NUM_REPS` times to get an average,
and gets the next problem size by multiplying by `MULTIPLIER`, until at most
`MAX_KEYS_PER_PROCESS`.  In your testing, you will probably want to test **one
problem size** at a time, using **a single node**, if possible.  If you are
having problems with correctness (segmentation faults, hangs/deadlocks, etc.),
it is best to work those out on your workstation/laptop is possible before
using SUs on Stampede2.  You are starting (knock on wood) from a correct
implementation: try to work in small changes, testing for correctness at each change.

- Instead of a report, your `proj2.pdf` will be a *notebook* chronicling your
  work.  Each entry should have the following format:
    * **A script for running ./test_proj2 and its output**:  This can be a
      submission script for Stampede2, or simply a shell script for another
      machine (like your laptop).  It can be for a different set of problem
      sizes than the final submission script; it can be for a single problem
      size.  What it should retain from the final submission script, in
      addition to a call to `./test_proj2`:
        * `get rev-parse HEAD` and `git diff-files`: you should be testing a
          *clean* version of your code.  Points will be deducted if there are
          uncommitted changes.
    * **Evidence from some form of profiling** (such as a screenshot/ terminal
      output from `pprof` or `jumpshot` from `tau`, `hpcprof` , `hpcviewer` or
      `hpctraceviewer` from `hpctoolkit``, `gprof`, statistics from your own
      instrumentation of the code, etc.) that you think indicates a
      hotspot/bottleneck in the code that could be improved.
    * **A short description of the changes you will make**, motivated by the
      evidence from profiling.
    * You don't have to do this for every commit (you can use more than one
      commit to implement a new feature, and your history doesn't have to be
      linear), but you should have a notebook entry for every major
      optimization strategy you undertake.

## Grading

- 0-3 points for hassle-free usage: maximized if the final batch script
  `slurm-proj2-stampede2-final.sh` runs the first time.
    * Points lost if we have to figure out how to reproduce your reported results.
- 0-6 points for correctness:
    * Whether the final script `slurm-proj2-stampede2-final.sh` runs to
      completion (it will abort if a list of keys is not properly sorted).
    * You lose half the points if your code is not correct; subsequent points
      can be lost for poor code organization.
- 0-9 Points for the notebook:
    * 0-3 points for how well the notebook tracks your `git` history: did we find the commits
      used to generate the entries?  Is there an entry for all the major
      aspects of your development?
    * 0-3 points for your profiling evidence: is it present?  Does it seem to
      indicate what you say it indicates?
    * 0-3 points for your planning: do your proposed code changes follow
      logically from the evidence?
- 0-3 **Extra credit** for performance: one point for every 10% speedup over
  the starting code, up to three points.

