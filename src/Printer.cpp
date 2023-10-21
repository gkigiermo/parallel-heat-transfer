#include <mpi.h>
#include <iostream>
#include "Printer.h"

Printer::Printer(CommunicationScheme::Process& process, int sizeX, int sizeY, int processX, int processY, int hasToPrintOnScreen){
    proc = process;
    tnx = sizeX;
    tny = sizeY;
    npx = processX;
    npy = processY;
    printOnScreen = hasToPrintOnScreen;
    tfinal = new double[tnx * tny];
}

Printer::~Printer(){
        delete[] tfinal;
}   

void Printer::printOnScreenAndTecplot(double *t, int plast, double dx, double dy)
{
    int i,j;
    //When there is only one process
    if(npx*npy==1){
        if (printOnScreen){
		    printOnScreenUnique(t);
        }
        printOnTecplot(t, 1, dx, dy);
            
	}else { // When there are more than one process
        //Communicate the info to a single process 
		joinInfo(t, plast);
        if (printOnScreen){
            printOnScreenMultiple();    
        }
        printOnTecplot(tfinal, 0, dx, dy);
	}
}

// To print the info of process 0 on the screen
void Printer::printOnScreenMultiple(){
    
    int i,j;
    if( proc.rank == 0) {
        for(j = tny-1; j >= 0; j--){
            for(i = 0; i < tnx; i++){
                std::cout << " " << tfinal[i + j*tnx]<<" ";
            }
            std::cout << std::endl;
        }
    }		
}

//Gather the info of all the processes in process 0
void Printer::joinInfo(double *t, int plast)
{
    int i,j,k;;
    long int ctx,cty;
    int tag=0;
    MPI_Status msg;

    if(proc.rank==0) {
        k=0;
        for(j = 0; j < (tny/npy); j++){
            for(i = 0; i < (tnx/npx); i++){
                tfinal[i + j * tnx]=t[i + proc.startx + (j + proc.starty) * proc.lgx];
            }
        }

        MPI_Send(tfinal, tnx * tny, MPI_DOUBLE, proc.rank + 1, tag, MPI_COMM_WORLD); 
        MPI_Recv(tfinal, tnx * tny, MPI_DOUBLE, plast, tag, MPI_COMM_WORLD, &msg);
    } else {
        if(proc.rank == plast) {
            MPI_Recv(tfinal, tnx * tny, MPI_DOUBLE, proc.rank-1, tag,MPI_COMM_WORLD,&msg);


            ctx=proc.rank%npx;
            ctx=ctx*(tnx/npx);
            cty=(long int)proc.rank/npx;
            cty=cty*tnx*(tny/npy);
            for(j=0;j<(tny/npy);j++)
                for( i=0;i<(tnx/npx);i++)
                    tfinal[i+ctx+ j*tnx + cty]=t[i+proc.startx+ (j+proc.starty)*proc.lgx];


            MPI_Send(tfinal,tnx*tny,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);


        }else {
            MPI_Recv(tfinal,tnx*tny,MPI_DOUBLE,proc.rank-1,tag,MPI_COMM_WORLD,&msg);
            ctx = proc.rank % npx;
            ctx = ctx*(tnx/npx);
            cty = (long int)proc.rank/npx;
            cty = cty*tnx*(tny/npy);
            for(j=0;j<(tny/npy);j++){
                for( i=0;i<(tnx/npx);i++){
                    tfinal[i+ctx+ j*tnx + cty]=t[i+proc.startx+ (j+proc.starty)*proc.lgx];
                }
            }
            MPI_Send(tfinal, tnx*tny, MPI_DOUBLE, proc.rank + 1, tag, MPI_COMM_WORLD);
        }
    }
}


void Printer::printOnScreenUnique(double *t){
	for(int j = proc.endy - 1; j >= proc.starty; j--){
		for(int i = proc.startx; i < proc.endx; i++)
			std::cout  << " " << t[XY(i,j,proc)] << " ";
		std::cout  <<  std::endl;
	}
}
void Printer::printOnTecplot(double* temperature, int option, double dx, double dy){
    FILE *fp1;

    int i,j;

    fp1=fopen("temperatures.dat" , "w+");
    fprintf(fp1, "VARIABLES=X,Y,T\n");
    fprintf(fp1, "ZONE I=%d J=%d \n", tnx, tny);

    if(option) {
        for(j = 0; j < tny; j++){
            for(i = 0; i < tnx; i++){
                fprintf(fp1, "%f %f %f \n",i*dx, j*dy, temperature[XY(i + proc.startx, j + proc.starty, proc)]);
            }
        }
    } else {
        for(j = 0; j < tny; j++) {
            for(i = 0; i < tnx; i++) {
                fprintf(fp1, "%f %f %f \n",i*dx,j*dy,temperature[i + j*tnx]);
            }
        }
    }
    fclose(fp1);
}
