#include <cstdio>
#include <cstdlib>

#include "lib_mem.h"

int main(int argc, char *argv[]){

    system("clear");
    fprintf(stdout, "----------------------\n");
    fprintf(stdout, "Test module for libmem\n");
    fprintf(stdout, "----------------------\n");

    const std::size_t xs = 64;
    const std::size_t ys = 64;
    const std::size_t zs = 64;
    const std::size_t ws = 64;

    fprintf(stdout, "(Info):\tAllocating memory\n");
    fprintf(stdout, "\t1D array: ");
    double*    data1D = memory<double>::allocate(xs);
    fprintf(stdout, "[%s]\n", data1D == nullptr ? "Failed" : "Done");

    for(int i = 0; i < xs; i++){
        data1D[i] = 1.0*i;
    }
    memory<double>::deallocate(data1D);

    fprintf(stdout, "\t2D array: ");
    double**   data2D = memory<double>::allocate(xs, ys);
    fprintf(stdout, "[%s]\n", data1D == nullptr ? "Failed" : "Done");
    for(int i = 0; i < xs*ys; i++){
        data2D[0][i] = 1.0*i;
    }
    memory<double>::deallocate(data2D);

    fprintf(stdout, "\t3D array: ");
    double***  data3D = memory<double>::allocate(xs, ys, zs);
    fprintf(stdout, "[%s]\n", data1D == nullptr ? "Failed" : "Done");
    for(int i = 0; i < xs*ys; i++){
        data3D[0][0][i] = 1.0*i;
    }
    memory<double>::deallocate(data3D);

    fprintf(stdout, "\t4D array: ");
    double**** data4D = memory<double>::allocate(xs, ys, zs, ws);
    fprintf(stdout, "[%s]\n", data1D == nullptr ? "Failed" : "Done");
    for(int i = 0; i < xs*ys; i++){
        data4D[0][0][0][i] = 1.0*i;
    }
    memory<double>::deallocate(data4D);

    fprintf(stdout, "(Info):\tMemory test complete\n");
    return(EXIT_SUCCESS);
}
