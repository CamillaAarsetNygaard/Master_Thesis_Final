#include "basics.h"
#include "update.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int checkCollision(double* x, double* y, int index, double* RPlist){
    int j;
    for (j=0; j<index; j++){
        if (wall==0){
            if ((sqrtf((x[index]-x[j])*(x[index]-x[j])+(y[index]-y[j])*(y[index]-y[j]))<(RPlist[index]+RPlist[j])) || ((sqrtf(x[index]*x[index]+y[index]*y[index])+RPlist[index])>=R)){
            return 1;
        }}
        else if (wall==1){
            if ((sqrtf((x[index]-x[j])*(x[index]-x[j])+(y[index]-y[j])*(y[index]-y[j]))<(RPlist[index]+RPlist[j])) || (x[index]+RPlist[index]>=R) || (x[index]-RPlist[index]<=-R) || (y[index]+RPlist[index]>=R) || (y[index]-RPlist[index]<=-R)){
            return 1;
        }} 
        else if (wall==2){
            if ((sqrtf((x[index]-x[j])*(x[index]-x[j])+(y[index]-y[j])*(y[index]-y[j]))<(RPlist[index]+RPlist[j])) || (y[index]-RPlist[index]<=-R*sn30) || (y[index]-((-(sn30+1))/(cs30))*x[index]>=R-RPlist[index]/(sn30)) || (y[index]-((sn30+1)/(cs30))*x[index]>=R-RPlist[index]/(sn30))){
            return 1;
        }}
        else if (wall==4) {
            double theta = atan2(y[index], x[index]);
            double Rellipse = pow(pow((fabs(cosf(theta)/a)), N)+pow((fabs(sinf(theta)/b)), N), -1/N);
            // double Rellipse = pow(abs(x[index]/a), N) + pow(abs(y[index]/b), N);
            if ((sqrtf((x[index]-x[j])*(x[index]-x[j])+(y[index]-y[j])*(y[index]-y[j]))<(RPlist[index]+RPlist[j])) || ((sqrtf(x[index]*x[index]+y[index]*y[index])+RPlist[index])>=Rellipse)){
            return 1;
        }}
    }
    return 0;
}


void resetCell(int*** cellList, int** numCell, int NCx, int NCy) { //make cellList "empty" so we can fill it up again
    int i, j, k;
    for (i=0; i<NCx; i++) {
        for (j=0; j<NCy; j++) {
            for (k=0; k<20; k++) {
            cellList[i][j][k] = -1;}
            numCell[i][j] = 0;
        } 
    }
}

void placeInCell(int*** cellList, int** numCell, double* x, double* y, int NCx, int NCy){
    int i;
    int NCxi, NCyi;
    resetCell(cellList, numCell, NCx, NCy);
    for (i=0; i<NP; i++) {
        findNCxiyi(&NCxi, &NCyi, x, y, i);
        numCell[NCxi][NCyi] += 1;
        cellList[NCxi][NCyi][numCell[NCxi][NCyi]-1] = i;
    }
}

void makeRadiuses(double* RPlist, double* totAreaPart){ //function to make the radiuses of the particles
    int i;
    if (diffRadius==0){ //if all particles have the same radius
        for (i=0; i<NP; i++){
            RPlist[i] = RP;
            *totAreaPart += PI*RP*RP;
        }
    }
    if (diffRadius==1){ //if the particles have different radiuses
        for (i=0; i<NP; i++){
            RPlist[i] = RP-radDev + (double)rand() / RAND_MAX * (2*radDev);
            *totAreaPart += PI*RPlist[i]*RPlist[i];
        }
    }
}
