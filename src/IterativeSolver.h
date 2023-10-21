#ifndef ITERATIVESOLVER_H
#define ITERATIVESOLVER_H
#include "CommunicationScheme.h"

class IterativeSolver {
    struct DiagonalMatrix{
        double *ap,*an,*as,*ae,*aw,*b; // 2-order discretization in cartessian coordinates
    };

    public:
    DiagonalMatrix matrix;
    CommunicationScheme::Process proc;

    IterativeSolver(CommunicationScheme::Process& process, int temperatureL, int temperatureR);
    ~IterativeSolver();
    int solve2DGS(double *);
    int solve2DJacobi(double *,double *);
    

    private:
    void initMatrix(int T1, int T2);
    void freeMatrix();
};
#endif