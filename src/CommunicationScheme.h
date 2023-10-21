#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#define XY(i, j, p)  (i + (j) * p.lgx)
#endif

#ifndef COMMUNICATIONSCHEME_H
#define COMMUNICATIONSCHEME_H

class CommunicationScheme{
public:
    struct Process{
	
	    int rank; //process id
        int startx, endx; //start and end of x corresponding to the processor
        int starty, endy; //start and end of y corresponding to the processor
        int pe,pw,pn,ps;//who are the neighbors
        int lgx,lgy; // length of the grid in x and y
        bool nbe,nbw,nbn,nbs; // flags to know if the process has neighbors in the corresponding direction
    };
    struct Halos{

        double *rhw,*rhe,*rhn,*rhs;// auxiliary vectors to receive the halos
        double *shw,*she,*shn,*shs;// auxiliary vectors to send the halos
    };


    Process proc;
    Halos halos;
    int npx, npy, tnx, tny, h;

    CommunicationScheme(int partitionsX, int partitionsY, int totalX, int totalY, int haloSize);
    ~CommunicationScheme();
    
    void update(double* vector);
    int last();

private:
    void initProcess();
    void initHalos();
    void extractHalos(double* vector);
    void insertHalos(double* vector);
    void exchangeHalos();
    void freeHalos();
};
#endif