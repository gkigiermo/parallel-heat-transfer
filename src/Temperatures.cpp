#include <math.h>
#include "Temperatures.h"

Temperatures::Temperatures(CommunicationScheme::Process& process, double initialTemperature, double temperatureLeft, double temperatureRight){
    proc = process;
    t = new double[proc.lgx * proc.lgy];
    tant = new double[proc.lgx * proc.lgy];
    tjacob = new double[proc.lgx * proc.lgy];
    initTemperatures(initialTemperature, temperatureLeft, temperatureRight);
}

Temperatures::~Temperatures(){
    freeTemperatures();
}

// To initialize the temperatures
void Temperatures::initTemperatures(double Tini, double T1, double T2){

    int i,j;
	//Initialize the temperatures with Tini value
    for(j=0;j< proc.lgy;j++){
		for(i=0;i< proc.lgx;i++){
		    t[XY(i,j,proc)]=Tini;
            tant[XY(i,j,proc)]=Tini;        
        }
    }
    //Initialize the temperatures with T1 and T2 values
    i= proc.startx;
    if(! proc.nbw){
        for(j= proc.starty;j< proc.endy;j++){
            t[XY(i,j,proc)]=T1;
            tant[XY(i,j,proc)]=T1;
        }
    }

    i= proc.endx-1;
    if(! proc.nbe){
        for(j= proc.starty;j< proc.endy;j++){
            t[XY(i,j,proc)]=T2;
            tant[XY(i,j,proc)]=T2;
        }
    }
}


//Freeing the temperatures vectors
void Temperatures::freeTemperatures() {
    delete[] t;
    delete[] tant;
    delete[] tjacob;
}

void Temperatures::updateTemperature(){

    for(int j = 0; j < proc.lgy; j++) {
	    for(int i = 0; i < proc.lgx; i++) {
		    tant[XY(i,j,proc)] = t[XY(i,j,proc)];
        }
    }
}

void Temperatures::updateTemperatureJacobi(){

    for(int j = 0; j < proc.lgy; j++) {
	    for(int i = 0; i < proc.lgx; i++) {
		    tjacob[XY(i,j,proc)] = t[XY(i,j,proc)];
        }
    }
}

//Calculates the error on each process
double Temperatures::calculateError(){
    double err=0;
    for(int j = 0; j < proc.lgy; j++){
	    for(int i = 0; i < proc.lgx; i++){
		    if(fabs(t[XY(i,j,proc)] - tant[XY(i,j,proc)]) > err) {
			    err = fabs(t[XY(i,j,proc)] - tant[XY(i,j,proc)]);
            }
        }
    }
    return err;
}
