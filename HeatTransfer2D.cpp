#include <iostream>
#include <mpi.h>
#include <time.h>
#include "src/CommunicationScheme.h"
#include "src/IterativeSolver.h"
#include "src/Temperatures.h"
#include "src/Printer.h"

//The domain is divided in npx*npy processes
#define npx 2   // number of divisions in X
#define npy 2   // number of divisions in Y

#define tnx 200  // number of control volumes in X
#define tny 200 // number of control volumes in Y
#define h 1     // ghost cells added to each side

//Constants from the heat transfer problem
#define K 100.0
#define T1 0.0 // Temperature at the left boundary
#define T2 1.0 // Temperature at the right boundary
#define Lx 1.0 // Length of the domain in X
#define Ly 1.0 // Length of the domain in Y
#define Tini 0.0 // Initial temperature

#define errorParallel 1e-7 // Desired precision in the solution
#define SCHEME 1 // 1: GS 2: Jacobi
#define PRINT_SCREEN 0 // 0: No 1: Yes

int main(int argc,char ** argv){
    
    MPI_Init(&argc, &argv);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Verifying that MPI has been executed with the correct number of processes
    if(size != npx * npy){
        std::cout << "Error: Number of processes is not equal to npx*npy" << std::endl;
        MPI_Finalize();
        return 0;
    }
    //Calculating the size of the cells
    double dx = Lx / (tnx - 1.0);
    double dy = Ly / (tny - 1.0);

    //Dividing the domain in npx*npy processes
    CommunicationScheme commScheme(npx, npy, tnx, tny, h);
    //Setting up the solver
    IterativeSolver solver(commScheme.proc, T1, T2);
    //Setting up the temperatures
    Temperatures temperatures(commScheme.proc, Tini, T1, T2);
    //Setting up the printer to print on screen and tecplot
    Printer printer(commScheme.proc, tnx, tny, npx, npy, PRINT_SCREEN);
    
    time_t initialTime, finalTime;
    double totalTime, averageTime;
    initialTime = time(NULL);

    int i=0,k=0;
    double error,maximumError;

	do {
        temperatures.updateTemperature();
        if(SCHEME == 1){
		    k = k + solver.solve2DGS(temperatures.getTemperaturesT());
        } else {
            temperatures.updateTemperatureJacobi();
            k = k + solver.solve2DJacobi(temperatures.getTemperaturesT(), temperatures.getTemperaturesTjacob());
        }
        commScheme.update(temperatures.getTemperaturesT());
        //calculating the local error
    	error = temperatures.calculateError();
        //calculating the maximum error among processes
		MPI_Allreduce(&error, &maximumError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
        //print the error every 1000 iterations
		if(i%1000 == 0 && commScheme.proc.rank == 0){
	        std::cout << "Iteration:" << i << " Residue:" << maximumError << std::endl;
		}
        i++;
	} while(maximumError > errorParallel);

    finalTime=time(NULL);
    totalTime = difftime(finalTime, initialTime);
    MPI_Allreduce(&totalTime, &averageTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    averageTime = averageTime / (npx * npy);

    //Printing the results on screen and generate a tecplot file
    printer.printOnScreenAndTecplot(temperatures.getTemperaturesT(), commScheme.last(),dx, dy);
        
    if(commScheme.proc.rank==0) {
        std::cout << "Final residue: " << maximumError << " Number of parallel iterations: " << i << " Total number of iterations: " << k << std::endl;
        std::cout << std::endl << "Average time: " << averageTime << " Time of Proc 0: " << totalTime << std::endl;
    }

	MPI_Finalize();
}
