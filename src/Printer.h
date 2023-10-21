#ifndef PRINTER_H
#define PRINTER_H
#include "CommunicationScheme.h"

class Printer{
public:
    Printer(CommunicationScheme::Process& process, int sizeX, int sizeY, int processX, int processY, int hasToPrintOnScreen);
    ~Printer();
    void printOnScreenAndTecplot(double *t, int plast, double dx, double dy);

private:
    void printOnScreenMultiple();
    void printOnScreenUnique(double *t);
    void printOnTecplot(double* temperature, int plast, double dx, double dy);
    void joinInfo(double *t, int plast);

    CommunicationScheme::Process proc;
    double* tfinal;
    int tnx, tny;
    int npx, npy;
    int printOnScreen;
};
#endif