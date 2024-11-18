#ifndef _LIBMPI_
#define _LIBMPI_

#define MPI_RANK_ROOT 0
#define MPI_RANK_TASK 1

#include "mpi.h"

typedef struct mpi_cmds{
    static int task;
    static int kill;
} mpi_cmds;

int mpi_cmds::task = 2;
int mpi_cmds::kill = 1;

typedef struct mpi_pmsg{
    static int ready;
    static int error;
    static int warning;
} mpi_pmsg;

int mpi_pmsg::ready   = 1;
int mpi_pmsg::error   = 2;
int mpi_pmsg::warning = 3;

#endif
