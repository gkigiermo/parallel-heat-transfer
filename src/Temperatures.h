#ifndef TEMPERATURES_H
#define TEMPERATURES_H
#include "CommunicationScheme.h"

class Temperatures{

public:
    Temperatures(CommunicationScheme::Process& process, double initialTemperature, double temperatureLeft, double temperatureRight);
    ~Temperatures();
    void updateTemperature();
    void updateTemperatureJacobi();
    double calculateError();
    double* getTemperaturesT(){return t;};
    double* getTemperaturesTant(){return tant;};
    double* getTemperaturesTjacob(){return tjacob;};
private:
    void initTemperatures(double Tini, double T1, double T2);
    void freeTemperatures();
    double *t, *tant, *tjacob;
    CommunicationScheme::Process proc;

};
#endif