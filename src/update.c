#include "update.h"
#include "basics.h"
#include "wall.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void update(double* x, double* y, double* q, double* fx, double* fy, double* tq,double* fxold, double* fyold, double* tqold, int*** cellList, int** numCell, int* firstRun, int NCx, int NCy, double* RPlist, int sample){
    int i;
    placeInCell(cellList, numCell, x, y, NCx, NCy);
    forceUpdate(x, y, q, fx, fy, tq, cellList, numCell, NCx, NCy, RPlist, sample);
    if (method == 0){
        for (i=0; i<NP; i++) {
            x[i] += fx[i]*DT;
            y[i] += fy[i]*DT;
            q[i] += tq[i]*DT+sqrtf(DIFF*DT)*12*(drand48()-0.5);
        }
    }
    if (method == 1) {
        for (i=0; i<NP; i++) {
            if (*firstRun==1) {
                x[i] += (fx[i])*DT;
                y[i] += (fy[i])*DT;
                q[i] += (tq[i])*DT;
                *firstRun = 0;
            }
            else {
                x[i] += (1.5*fx[i]-0.5*fxold[i])*DT;
                y[i] += (1.5*fy[i]-0.5*fyold[i])*DT;
                q[i] += (1.5*tq[i]-0.5*tqold[i])*DT;
            }
            fxold[i] = fx[i];
            fyold[i] = fy[i];
            tqold[i] = tq[i];
        }
    }
}

void forceUpdate(double* x, double* y, double* q, double* fx, double* fy, double* tq, int*** cellList, int** numCell, int NCx, int NCy, double* RPlist, int sample) {
    int i, j, k, l, NCxi, NCyi;
    double fljix, fljiy, fwallx, fwally, dv;

    for (i=0; i<NP; i++) {
        fwallx=0; 
        fwally=0;
        fljix=0;
        fljiy=0;

        //finding the cell of the particle
        findNCxiyi(&NCxi, &NCyi, x, y, i);
        //The 8 cells around the cell of the particle + the cell itself
        int neighCell[9][2] = {{NCxi-1, NCyi-1}, {NCxi, NCyi-1}, {NCxi+1, NCyi-1}, {NCxi-1, NCyi}, {NCxi, NCyi}, {NCxi+1, NCyi}, {NCxi-1, NCyi+1}, {NCxi, NCyi+1}, {NCxi+1, NCyi+1}};
        //Interaction with the other particles in the cell
        for (l=0; l<9; l++) {
            for (j=0; j<numCell[neighCell[l][0]][neighCell[l][1]]; j++) {
                k=cellList[neighCell[l][0]][neighCell[l][1]][j];
                if (i!=k) {
                    dv = FHertz(x[i], y[i], x[k], y[k], RPlist[i], RPlist[k]);
                    fljix += dv*(x[i]-x[k]);
                    fljiy += dv*(y[i]-y[k]);

                }
            }
        }
        //Interaction with the wall
        if (wall==0) {
            circWall2(x, y, &fwallx, &fwally, i, RPlist);
        }
        else if (wall==1) {
            squareWall(x, y, &fwallx, &fwally, i, RPlist);
        }
        else if (wall==2) {
            rectWall(x, y, &fwallx, &fwally, i, RPlist);
        }
        else if (wall==3) {
            triWall(x, y, &fwallx, &fwally, i, RPlist);
        }
        else if (wall==4) {
            superEllipse2(x, y, &fwallx, &fwally, i, RPlist);
        }
        double sampleFrac =(1.0 - ((double )sample)/((double )totalSamples));
        fx[i]=sampleFrac*FP*cosf(q[i])+fljix+fwallx;
        fy[i]=sampleFrac*FP*sinf(q[i])+fljiy+fwally;
        tq[i]=XI*(-sinf(q[i])*(fljix+fwallx)+cosf(q[i])*(fljiy+fwally));
        
    }
}

void updateFile(double* x, double* y, double* q,  double* xtilde, double* ytilde, double* qtilde, double* fx, double* fy, double* tq, double* fxold, double* fyold, double* tqold, int* firstRun, double* RPlist, int** M, double** r){
    int i;
    forceUpdateFile(x, y, q, fx, fy, tq, RPlist, M, r);
    if (method == 0){
        for (i=0; i<NP; i++) {
            x[i] += fx[i]*DT;
            y[i] += fy[i]*DT;
            // q[i] += tq[i]*DT;
            q[i] += tq[i]*DT+sqrtf(DIFF*DT)*(drand48()-0.5);
        }
    }
    if (method == 1) {
        for (i=0; i<NP; i++) {
            if (*firstRun==1) {
                x[i] += (fx[i])*DT;
                y[i] += (fy[i])*DT;
                q[i] += (tq[i])*DT;
                *firstRun = 0;
            }
            else {
                x[i] += (1.5*fx[i]-0.5*fxold[i])*DT;
                y[i] += (1.5*fy[i]-0.5*fyold[i])*DT;
                q[i] += (1.5*tq[i]-0.5*tqold[i])*DT;
            }
            fxold[i] = fx[i];
            fyold[i] = fy[i];
            tqold[i] = tq[i];
        }
    }
    if (method ==2) {
        for (i=0; i<NP; i++) {
            xtilde[i] = x[i]+fx[i]*DT*0.5;
            ytilde[i] = y[i]+fy[i]*DT*0.5;
            qtilde[i] = q[i]+tq[i]*DT*0.5;
        }
        forceUpdateFile(xtilde, ytilde, qtilde, fx, fy, tq, RPlist, M, r);
        for (i=0; i<NP; i++) {
            x[i] = x[i]+fx[i]*DT;
            y[i] = y[i]+fy[i]*DT;
            q[i] = q[i]+tq[i]*DT;
        }
    }
}

void forceUpdateFile(double* x, double* y, double* q, double* fx, double* fy, double* tq, double* RPlist, int** M, double** r){
    int i, j, k;
    double fljix, fljiy, fwallx, fwally, dv;
    for (i=0; i<NP; i++) {
        fwallx=0; 
        fwally=0;
        fljix=0;
        fljiy=0;
        if (sqrt(x[i]*x[i]+y[i]*y[i]) + RPlist[i] >= R - RPlist[i]/2) {
            fwallx=0; 
            fwally=0;
            fljix=0;
            fljiy=0;
        }
        else{
        for (j=0; j<9; j++) {
            k = M[i][j]; //Neighbouring particle of particle i
            if (k==-1) {
                break;
            }
            dv = FSpring(x[i], y[i], x[k], y[k], r[i][j]);
            fljix += dv*(x[i]-x[k]);
            fljiy += dv*(y[i]-y[k]);
            }

        //Interaction with the wall
        if (wall==0) {
            circWall(x, y, &fwallx, &fwally, i, RPlist);
        }
        else if (wall==1) {
            squareWall2(x, y, &fwallx, &fwally, i, RPlist);
        }
        else if (wall==2) {
            rectWall(x, y, &fwallx, &fwally, i, RPlist);
        }
        else if (wall==3) {
            triWall(x, y, &fwallx, &fwally, i, RPlist);
        }
        else if (wall==4) {
            superEllipse2(x, y, &fwallx, &fwally, i, RPlist);
        }
        fx[i]=FP*cosf(q[i])+fljix+fwallx;
        fy[i]=FP*sinf(q[i])+fljiy+fwally;
        tq[i]=XI*(-sinf(q[i])*(fljix+fwallx)+cosf(q[i])*(fljiy+fwally)); }
    }
}


//Lennard-Jones force
double DV(double x1, double y1, double x2, double y2, double* RPlist, int index1, int index2){
    double r2 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
    double u2, v2;
    double cutoff2 = 2.9*2.9*(RPlist[index1]+RPlist[index2])*(RPlist[index1]+RPlist[index2]);
    if (r2<=cutoff2) {
        v2=r2/(cutoff2);
        // u2 = (2*RPev*2*RPev)/r2;
        u2 = ((RPlist[index1]+RPlist[index2])*(RPlist[index1]+RPlist[index2]))/r2;
        return 12*EPS*u2*u2*(u2-1)*(1/r2)*3.71*(1-v2*v2*v2);
    }
    else {
        return 0;
    }
    
}

//Hertz contact elastic force
double FHertz(double x1, double y1, double x2, double y2, double R1, double R2){
    double r = sqrtf((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    double h = R1+R2-r;
    if (h>0){
        return C*sqrtf(h*h*h*R1*R2/(R1+R2))/r;}
    else {
        return 0;}
}

//Spring force
double FSpring(double x1, double y1, double x2, double y2, double r){
    double rij = sqrtf((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return (K*(r-rij)+K4*(r-rij)*(r-rij)*(r-rij))/rij;
}



//Find the cell of the particle
void findNCxiyi(int* NCxi, int* NCyi, double* x, double* y, int index){
    if ((wall==0) | (wall== 1)) {
        *NCxi = (int) floor((x[index]+R+DELTA)/DELTA);                 
        *NCyi = (int) floor((y[index]+R+DELTA)/DELTA);
    }
    else if (wall==2) {
        *NCxi = (int) floor((x[index]+rectWidthHalf+DELTA)/DELTA);  
        *NCyi = (int) floor((y[index]+rectHeightHalf+DELTA)/DELTA);
        if (x[index] > rectWidthHalf) {
                    printf("Particle %d is too far out \n", index);
                    printf("x: %f\n", x[index]);
                    printf("y: %f\n", y[index]);
                }
    }
    else if (wall==3) {
        *NCxi = (int) floor(x[index]/DELTA);
        *NCyi = (int) floor(y[index]/DELTA);
    }
    else if (wall==4) {
        *NCxi = (int) floor((x[index]+a+DELTA)/DELTA);                 
        *NCyi = (int) floor((y[index]+b+DELTA)/DELTA);
    }                       
}
