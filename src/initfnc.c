#include "initfnc.h"
#include "basics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void initxyq(double* x, double* y, double* q, double *RPlist){
    int i;
    srand((unsigned int) time(NULL));
    if (wall==0){
        x[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
        y[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
        q[0]=((float) rand() / (float)(RAND_MAX))*2*PI;
        while (sqrt(x[0]*x[0]+y[0]*y[0])>=R-RPlist[0]){
            x[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
            y[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
        }
        for (i=1; i<NP; i++) {
            x[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
            y[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
            q[i]=((float) rand() / (float)(RAND_MAX))*2*PI;
            while (checkCollision(x, y, i, RPlist)==1){
                x[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
                y[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
    }}} 
    else if (wall==1){
        for (i=1; i<NP; i++) {
            x[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(R-RPlist[i]);
            y[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(R-RPlist[i]);
            q[i]=((float) rand() / (float)(RAND_MAX))*2*PI;
            while (checkCollision(x, y, i, RPlist)==1){
                x[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(R-RPlist[i]);
                y[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(R-RPlist[i]);
    }}}
    else if (wall==2){ //First particle always in the middle to make my life easier
        x[0]=0;
        y[0]=0;
        q[0]=((float) rand() / (float)(RAND_MAX))*2*PI;
        for (i=1; i<NP; i++) {
            x[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(R-RPlist[i]);
            y[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(R-RPlist[i]);
            q[i]=((float) rand() / (float)(RAND_MAX))*2*PI;
            while (checkCollision(x, y, i, RPlist)==1){
                x[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(R-RPlist[i]);
                y[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(R-RPlist[i]);
    }}}
    else if (wall==4) {
        x[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
        y[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
        q[0]=((float) rand() / (float)(RAND_MAX))*2*PI;
        double theta = atan2(y[0], x[0]);
        double Rellipse = pow(pow((fabs(cosf(theta)/a)), N)+pow((fabs(sinf(theta)/b)), N), -1/N);
        while (sqrt(x[0]*x[0]+y[0]*y[0])>=Rellipse-RPlist[0]){
            x[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
            y[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
        }
        for (i=1; i<NP; i++) {
            x[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
            y[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
            q[i]=((float) rand() / (float)(RAND_MAX))*2*PI;
            while (checkCollision(x, y, i, RPlist)==1){
                x[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
                y[i]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
    }}}
}

void initfile(double* x, double* y, double* q, double* RPlist, char *fileName, char * tesselationName, int **M, double** r){
    FILE *fptr;
    FILE *tessptr2;
    int lenTess;
    int p1, p2;
    int i, j, k;
    double rij;
    fptr = fopen(fileName,"r");
    tessptr2 = fopen(tesselationName,"r");

    //Initializing the M and r arrays
    for (i=0; i<NP; i++) {
        for (j=0; j<9; j++){
            M[i][j] = -1;
            r[i][j] = 0;}
    }
    
    //Make the arrays x, y, q and RPlist
    for (i=0; i<NP; i++) {
        fscanf(fptr, "%lf %lf %lf %lf", &x[i], &y[i], &q[i], &RPlist[i]);
    }

    //Find the length of the tessellation file
    fscanf(tessptr2, "%d %d", &lenTess, &p1);

    //Make the M array
    for (i=0; i<lenTess; i++) {
        fscanf(tessptr2, "%d %d", &p1, &p2);
        for (j=0; j<9; j++) {
            rij = sqrtf((x[p1]-x[p2])*(x[p1]-x[p2])+(y[p1]-y[p2])*(y[p1]-y[p2]));
            if (M[p1][j] == -1 ) {
                M[p1][j] = p2;
                r[p1][j] = rij;
                for (k=0; k<7; k++) {
                    if (M[p2][k] == -1) {
                        M[p2][k] = p1;
                        r[p2][k] = rij;
                        break;
                        }
                    }
                break;}
        }
    }


    fclose(fptr);
    fclose(tessptr2);
    //randomize the q values
    for (i=1; i<NP; i++) {
        q[i]=((float) rand() / (float)(RAND_MAX))*2*PI;}
}



void inithex(double* x, double* y, double* q, double *RPlist){
    int i, j, k;
    k = 0;
    for (i=0; i<SQNP; i++) {
        for (j=0; j<SQNP; j++) {
            x[k] = j*RPlist[k]*2+(i % 2)*RPlist[k]-SQNP*RPlist[k]+RPlist[k]/2;
            y[k] = i*RPlist[k]*SQ3-((SQ3*0.5)*(SQNP-1)*RPlist[k]+RPlist[k])+RPlist[k];
            q[k] = ((float) rand() / (float)(RAND_MAX))*2*PI;
            k += 1;
        }
    }}

