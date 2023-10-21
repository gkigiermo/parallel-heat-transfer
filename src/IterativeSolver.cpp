#include <math.h>
#include "IterativeSolver.h"

IterativeSolver::IterativeSolver(CommunicationScheme::Process& process, int temperatureL, int temperatureR){
    proc = process;
    initMatrix(temperatureL, temperatureR);
}

IterativeSolver::~IterativeSolver(){
    freeMatrix();
}

//Allocating and initializing the matrix with the coefficients of the heat equation
void IterativeSolver::initMatrix(int T1, int T2){
    int i,j;

    matrix.ap= new double[proc.lgx * proc.lgy];
    matrix.an= new double[proc.lgx * proc.lgy];
    matrix.as= new double[proc.lgx * proc.lgy];
    matrix.aw= new double[proc.lgx * proc.lgy];
    matrix.ae= new double[proc.lgx * proc.lgy];
    matrix.b = new double[proc.lgx * proc.lgy];

    //Interior nodes
    for(i = proc.startx; i < proc.endx; i++){
        for(j = proc.starty; j < proc.endy; j++){
            matrix.ae[XY(i,j,proc)] = 1.0;
            matrix.aw[XY(i,j,proc)] = 1.0;
            matrix.an[XY(i,j,proc)] = 1.0;
            matrix.as[XY(i,j,proc)] = 1.0;
            matrix.ap[XY(i,j,proc)] = matrix.as[XY(i,j,proc)] + matrix.an[XY(i,j,proc)] + 
                                      matrix.aw[XY(i,j,proc)] + matrix.ae[XY(i,j,proc)];
            matrix.b[XY(i,j,proc)] = 0.0;
        }
    }
    //Left Boundary
    i = proc.startx;
    if(!proc.nbw){
        for(j = proc.starty; j < proc.endy; j++){
            matrix.ae[XY(i,j,proc)] = 0.0;
            matrix.aw[XY(i,j,proc)] = 0.0;
            matrix.an[XY(i,j,proc)] = 0.0;
            matrix.as[XY(i,j,proc)] = 0.0;
            matrix.ap[XY(i,j,proc)] = 1.0;
            matrix.b[XY(i,j,proc)] = T1;      
        }
    }
    //Right Boundary
    i = proc.endx - 1; 
    if(!proc.nbe){
        for(j=proc.starty; j<proc.endy; j++) {
            matrix.ae[XY(i,j,proc)] = 0.0;
            matrix.aw[XY(i,j,proc)] = 0.0;
            matrix.an[XY(i,j,proc)] = 0.0;
            matrix.as[XY(i,j,proc)] = 0.0;
            matrix.ap[XY(i,j,proc)] = 1.0;
            matrix.b[XY(i,j,proc)] = T2;      
        }
    }
    //Top Boundary
    j = proc.endy - 1;
    if(!proc.nbn){
        for(i = proc.startx; i < proc.endx; i++){
            matrix.ae[XY(i,j,proc)] = 0.0;
            matrix.aw[XY(i,j,proc)] = 0.0;
            matrix.an[XY(i,j,proc)] = 0.0;
            matrix.as[XY(i,j,proc)] = 1.0;
            matrix.ap[XY(i,j,proc)] = 1.0;
            matrix.b[XY(i,j,proc)] = 0.0;      
        }
    }
    //Bottom Boundary
    j = proc.starty;
    if(!proc.nbs){
        for(i = proc.startx; i < proc.endx; i++) {
            matrix.ae[XY(i,j,proc)] = 0.0;
            matrix.aw[XY(i,j,proc)] = 0.0;
            matrix.an[XY(i,j,proc)] = 1.0;
            matrix.as[XY(i,j,proc)] = 0.0;
            matrix.ap[XY(i,j,proc)] = 1.0;
            matrix.b[XY(i,j,proc)] = 0.0;      
        }
    }
}

//Freeing the memory allocated for the matrix
void IterativeSolver::freeMatrix()
{
    delete[] matrix.an;
    delete[] matrix.as;
    delete[] matrix.ap;
    delete[] matrix.aw;
    delete[] matrix.ae;
    delete[] matrix.b;
}

// To solve the system of equations using Gauss Seidel
int IterativeSolver::solve2DGS(double *t){
    int i,j;
    int k=0;
    double max_diff,nbs;   
    max_diff=0;

    for(j = proc.starty; j < proc.endy; j++){
        for(i = proc.startx; i < proc.endx; i++) {
            nbs = (1.0/(double)matrix.ap[XY(i,j,proc)])*(matrix.an[XY(i,j,proc)]*t[XY(i,j+1,proc)] + 
                               matrix.as[XY(i,j,proc)]*t[XY(i,j-1,proc)] + matrix.aw[XY(i,j,proc)]*t[XY(i-1,j,proc)] +
                               matrix.ae[XY(i,j,proc)]*t[XY(i+1,j,proc)] + matrix.b[XY(i,j,proc)]);
            if(fabs(t[XY(i,j,proc)] - nbs) > max_diff) {
                max_diff = fabs(t[XY(i,j,proc)] - nbs);
            }
            t[XY(i,j,proc)] = nbs;
        }
    }
    k++;        
    return k;
}

//To solve the system of equations using Jacobi
int IterativeSolver::solve2DJacobi(double *t,double *tjacob){

    int i,j, k = 0;
    double max_diff,nbs;
    max_diff=0;

    for(j = proc.starty; j < proc.endy; j++) {
        for(i = proc.startx; i < proc.endx; i++) {
            nbs = (1.0/(double) matrix.ap[XY(i,j,proc)]) * (matrix.an[XY(i,j,proc)] * tjacob[XY(i,j+1,proc)] +
                                matrix.as[XY(i,j,proc)]  * tjacob[XY(i,j-1,proc)] + matrix.aw[XY(i,j,proc)] * tjacob[XY(i-1,j,proc)] +
                                matrix.ae[XY(i,j,proc)]  * tjacob[XY(i+1,j,proc)] + matrix.b[XY(i,j,proc)]);
            if(fabs(t[XY(i,j,proc)] - nbs) > max_diff) {
                max_diff = fabs(t[XY(i,j,proc)] - nbs);
            }
            t[XY(i,j,proc)] = nbs;
        }
    }

    for(j = proc.starty; j < proc.endy; j++) {
        for( i = proc.startx; i < proc.endx; i++) {
            tjacob[XY(i,j,proc)] = t[XY(i,j,proc)];    
        }
    }
    k++;

    return k;

}