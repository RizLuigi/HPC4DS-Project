#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
//#include <omp.h>

int main(int argc, char* argv[]) {
    int NUM_ZOOMS = strtod(argv[1], NULL);
    int NUM_THREADS = strtod(argv[2], NULL);

    int my_rank, comm_sz;
    MPI_Init( NULL , NULL);
    MPI_Comm_size( MPI_COMM_WORLD , &comm_sz);
    MPI_Comm_rank( MPI_COMM_WORLD , &my_rank);

    if (comm_sz <= NUM_ZOOMS) {
        // In questo caso i processi sono indipendenti: ognuno lavora su uno o più zoom
        for(int i=0; i<NUM_ZOOMS/comm_sz; i++) {
            int actual_zoom = my_rank*NUM_ZOOMS/comm_sz+i;
            printf("Process %d working on zoom %d\n", my_rank, actual_zoom);
        }
    }else{
        // In questo caso più processi lavorano su uno stesso zoom, quindi ognuno avrà un master di riferimento
        // Forse è il caso  di creare nuovi communicators in modo da poter utilizzare collective communications?
        // https://mpitutorial.com/tutorials/introduction-to-groups-and-communicators/
        int proc_per_zoom = comm_sz/NUM_ZOOMS;
        int actual_zoom = my_rank/proc_per_zoom;
        int my_master = my_rank - my_rank%proc_per_zoom;
        printf("Process %d working on zoom %d with master %d\n", my_rank, actual_zoom, my_master);
    }

    MPI_Finalize();
    return 0;
}
