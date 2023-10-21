#include "CommunicationScheme.h"
#include <mpi.h>

CommunicationScheme::CommunicationScheme(int partitionsX, int partitionsY, int controlVolumesX, int controlVolumesY, int halo_size) {
    npx = partitionsX;
    npy = partitionsY;
    tnx = controlVolumesX;
    tny = controlVolumesY;
    h   = halo_size;
    initProcess();
    initHalos();
}
CommunicationScheme::~CommunicationScheme() {
    freeHalos();
}

// To findout the neighbors of each process and the lenghts with halos size
void CommunicationScheme::initProcess(){

	MPI_Comm_rank(MPI_COMM_WORLD,&proc.rank);

	//Finding neighbors
	if(proc.rank >= npx) {
		proc.ps  = proc.rank - npx;
		proc.nbs = 1;		
	} else {
		proc.nbs = 0;
		proc.ps  = -1;
	}

	if(proc.rank < npx * (npy - 1)) {
		proc.pn  = proc.rank + npx;
		proc.nbn = 1;
	} else {
		proc.nbn =  0;
		proc.pn  = -1;
	}
	
	if(proc.rank % npx != 0) {
		proc.pw  = proc.rank - 1;
		proc.nbw = 1;
	} else {
		proc.nbw =  0;
		proc.pw  = -1;
	}

	if((proc.rank + 1) % npx != 0 ) {
		proc.pe  = proc.rank + 1;
		proc.nbe = 1;
	} else {
		proc.nbe =  0;
		proc.pe  = -1;
	}
	
	//Calculate lenghts
	proc.lgx = tnx/npx + 2*h; 
	proc.lgy = tny/npy + 2*h; 

	proc.startx = h;

	proc.endx = proc.lgx - h;
	proc.starty = h;
	proc.endy = proc.lgy - h;
}	

//To allocate memory for the halos
void CommunicationScheme::initHalos(){
	if(proc.nbs) {
		halos.rhs = new double[h * (tnx/npx)];
		halos.shs = new double[h * (tnx/npx)];
	}
	if(proc.nbn) {
		halos.rhn = new double[h * (tnx/npx)];
		halos.shn = new double[h * (tnx/npx)];
	}
	if(proc.nbw) {
		halos.rhw = new double[h * (tny/npy)];
		halos.shw = new double[h * (tny/npy)];
	}
	if(proc.nbe) {
		halos.rhe = new double[h * (tny/npy)];
		halos.she = new double[h * (tny/npy)];
	}	
}


//To perform the communication considering packing and unpacking
void CommunicationScheme::update(double *vector) {
    //packing
    extractHalos(vector);
    //communication
    exchangeHalos();
    //unpacking
    insertHalos(vector);
}

// To communicate information between processes  
// A classic 4-neighbors scheme is used with a red-black ordering and blocking communications
void CommunicationScheme::exchangeHalos() {
    int tag=0;
    MPI_Status msg;
    //If the process is even
    if(!proc.rank%2) {
        //Sends to the right and receives from the right
        if(proc.nbe){
            MPI_Send(halos.she,h*(tny/npy),MPI_DOUBLE,proc.pe,tag,MPI_COMM_WORLD);
            MPI_Recv(halos.rhe,h*(tny/npy),MPI_DOUBLE,proc.pe,tag,MPI_COMM_WORLD,&msg);
        }
        //Receives from the left and sends to the left
        if(proc.nbw){
            MPI_Recv(halos.rhw,h*(tny/npy),MPI_DOUBLE,proc.pw,tag,MPI_COMM_WORLD,&msg);
            MPI_Send(halos.shw,h*(tny/npy),MPI_DOUBLE,proc.pw,tag,MPI_COMM_WORLD);
        }
        //Sends to the top and receives from the top
         if(proc.nbn){
            MPI_Send(halos.shn,h*(tnx/npx),MPI_DOUBLE,proc.pn,tag,MPI_COMM_WORLD);
            MPI_Recv(halos.rhn,h*(tnx/npx),MPI_DOUBLE,proc.pn,tag,MPI_COMM_WORLD,&msg);
        }
        //Receives from the bottom and sends to the bottom
        if(proc.nbs){
            MPI_Recv(halos.rhs,h*(tnx/npx),MPI_DOUBLE,proc.ps,tag,MPI_COMM_WORLD,&msg);
            MPI_Send(halos.shs,h*(tnx/npx),MPI_DOUBLE,proc.ps,tag,MPI_COMM_WORLD);
        }
    } else { //If the process is odd
        //Receives from the left and sends to the left
        if(proc.nbw){
            MPI_Recv(halos.rhw,h*(tny/npy),MPI_DOUBLE,proc.pw,tag,MPI_COMM_WORLD,&msg);
            MPI_Send(halos.shw,h*(tny/npy),MPI_DOUBLE,proc.pw,tag,MPI_COMM_WORLD);
        }
        //Sends to the right and receives from the right
        if(proc.nbe){
            MPI_Send(halos.she,h*(tny/npy),MPI_DOUBLE,proc.pe,tag,MPI_COMM_WORLD);
            MPI_Recv(halos.rhe,h*(tny/npy),MPI_DOUBLE,proc.pe,tag,MPI_COMM_WORLD,&msg);
        }
        //Receives from the bottom and sends to the bottom
        if(proc.nbs){
            MPI_Recv(halos.rhs,h*(tnx/npx),MPI_DOUBLE,proc.ps,tag,MPI_COMM_WORLD,&msg);
            MPI_Send(halos.shs,h*(tnx/npx),MPI_DOUBLE,proc.ps,tag,MPI_COMM_WORLD);
        }
        //Sends to the top and receives from the top
        if(proc.nbn){
            MPI_Send(halos.shn,h*(tnx/npx),MPI_DOUBLE,proc.pn,tag,MPI_COMM_WORLD);
            MPI_Recv(halos.rhn,h*(tnx/npx),MPI_DOUBLE,proc.pn,tag,MPI_COMM_WORLD,&msg);
        }
    }
}

//Unpacks the halos into the grid
void CommunicationScheme::insertHalos(double *t) {
    int i,j,k;
    //Top halos
    if(proc.nbn) {
        k=0;
        for(j = proc.endy; j < proc.lgy; j++) {
            for(i = proc.startx; i < proc.endx; i++) {
                t[XY(i,j,proc)] = halos.rhn[k];
                k++;       
            }
        }
    }
    //Bottom halos
    if(proc.nbs) {
        k=0;
        for(j = 0; j < proc.starty; j++) {
            for(i = proc.startx; i < proc.endx; i++) {
                t[XY(i,j,proc)] = halos.rhs[k];
                k++;
            }
        }
    }
    //Left halos
    if(proc.nbw) {
        k=0;
        for(j = proc.starty; j < proc.endy; j++) {
            for(i = 0; i < proc.startx; i++) {
                t[XY(i,j,proc)] = halos.rhw[k];
                k++;
            }
        }
    }
    //Right halos
    if(proc.nbe) {
        k=0;
        for(j = proc.starty; j < proc.endy; j++) {
            for(i = proc.endx;i < proc.lgx; i++) {
                t[XY(i,j,proc)] = halos.rhe[k];
                k++;
            }
        }
    }
}

//Pack the halos from the grid to be communicated
void CommunicationScheme::extractHalos(double *t) {
    int i,j,k;
    //Top halos
    if(proc.nbn) {
        k=0;
        for(j = proc.endy-h ; j < proc.endy; j++)
            for(i = proc.startx; i < proc.endx; i++) {
                halos.shn[k] = t[XY(i,j,proc)];
                k++;
            }
    }
    //Bottom halos
    if(proc.nbs) {
        k=0;
        for(j = proc.starty; j < proc.starty + h; j++){
            for(i = proc.startx; i < proc.endx; i++) {
                halos.shs[k] = t[XY(i,j,proc)];
                k++;
            }
        }
    }
     //Left halos
    if(proc.nbw) {
        k=0;
        for(j = proc.starty; j < proc.endy; j++) {
            for(i = proc.startx; i < proc.startx + h; i++) {
                halos.shw[k] = t[XY(i,j,proc)];
                k++;
            }
        }
    }      
    //Right halos
    if(proc.nbe) {
        k=0;
        for(j = proc.starty; j < proc.endy; j++) {
            for(i = proc.endx - h; i < proc.endx; i++) {
                halos.she[k] = t[XY(i,j,proc)];
                k++;
            }            
        }
    }
}

//Freeing halos memory
void CommunicationScheme::freeHalos(){
    if(proc.nbs) {
        delete[] halos.rhs;
        delete[] halos.shs;
    }
    if(proc.nbn) {
        delete[] halos.rhn;
        delete[] halos.shn;
    }
    if(proc.nbw) {
        delete[] halos.rhw;
        delete[] halos.shw;
    }
    if(proc.nbe) {
        delete[] halos.rhe;
        delete[] halos.she;
    }
}

//To findout the last processor
int CommunicationScheme::last() {
    int last;
    MPI_Comm_size(MPI_COMM_WORLD, &last);
    return (last - 1);
}
