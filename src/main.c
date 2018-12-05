#include <assert.h>
#include <stdio.h>
#include <mpi.h>

#define IS_MAIN_THREAD (rank == 0)

int np, rank;

int main(int argc, char** argv) {
    assert(MPI_Init(&argc, &argv) == MPI_SUCCESS);
    assert(MPI_Comm_size(MPI_COMM_WORLD, &np) == MPI_SUCCESS);
    assert(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);

    printf("--- rank: %d\n", rank);

    if (IS_MAIN_THREAD) {
        assert(MPI_Finalize() == MPI_SUCCESS);
    }

    return 0;
}
