#include <stdio.h>
#include "basics.h"
#include "old.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>


double xd(double theta, const double Fp, const double F12x){
    double deriv = Fp*cos(theta)+F12x;
    return deriv;
}

double yd(double theta, const double Fp, const double F12x){
    double deriv = Fp*sin(theta)+F12x;
    return deriv;
}

double thetad(double theta, const double F12[], const double beta) {
    double deriv = beta*2;
    return deriv;
}

void eulerStep(double* yn, double h, double fd[], double tn, int sizey){
    // for (int i=0; i<sizey; i++) {
    //     yn[i] = yn[i] + h*fd[i](tn, yn);
    // }
}

void evol(double t0, int n, double Fp, double xyt10[], double xyt20[]) { // The xyt0s are the start position and orientation of the particles
    double y[6] = {xyt10[0], xyt20[0], xyt10[1], xyt20[1], xyt10[2], xyt20[2]};
    printArray2(y, 6);
}


void printArray2(double arr[], int size) {
    for (int i=0; i<size; i++) {
        printf("%f \n", arr[i]);
    }
}

void findNeighbour(double* x, double* y, int* M, int index, double RPev){
    int j;
    for (j=0; j<NP; j++){
        if (j!=index){
            if (sqrt((x[index]-x[j])*(x[index]-x[j])+(y[index]-y[j])*(y[index]-y[j]))<2.5*RPev){
                M[index*7]+=1;
                M[index*7+M[index*7]]=j;
            }
        }
    }
}

//function to find the cells of the particles
void findCell(double* x, double* y, int index, int* M, double RPev, int Ncol, int Nrow){
    int i;
    double xi, yi;
    for (i=0; i<NP; i++){
        if ((wall==0) | (wall==1)){
            xi=(x[i]+R)/DELTA;
            yi=(y[i])/DELTA;
        }
        else if (wall==2) {
            xi=(x[i])/DELTA;
            yi=(y[i])/DELTA;
        }
        
        printf("Kolonne %f\n", xi);
        printf("Rad %f\n", yi);}
    
}

double Fljx2(double x1, double y1, double x2, double y2, double RPev){
    double r = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return 12*EPS*(pow(2*RPev/r, 6)-pow(2*RPev/r, 4))*(x1-x2)/pow(r, 2)*3.71*exp(-pow(r/(2*RPev), 2));
}

double Fljy2(double x1, double y1, double x2, double y2, double RPev){
    double r = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return 12*EPS*(pow(2*RPev/r, 6)-pow(2*RPev/r, 4))*(y1-y2)/pow(r, 2)*3.71*exp(-pow(r/(2*RPev), 2));}


double Fljx(double x1, double y1, double x2, double y2, double RPev){
    double r = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    if (r<=CUTOFF) {
        return 12*EPS*(pow(2*RPev/r, 6)-pow(2*RPev/r, 4))*((x1-x2)/pow(r, 2))*3.71*exp(-pow(r/(2*RPev), 2))*(1-pow((r/CUTOFF), 6));
    }
    else {
        return 0;
    }
    
}

double Fljy(double x1, double y1, double x2, double y2, double RPev){
    double r = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    if (r<=CUTOFF) {
        return 12*EPS*(pow(2*RPev/r, 6)-pow(2*RPev/r, 4))*((y1-y2)/pow(r, 2))*3.71*exp(-pow(r/(2*RPev), 2))*(1-pow((r/CUTOFF), 6));
    }
    else {
        return 0;
    }
    // return 12*EPS*(pow(2*RPev/r, 6)-pow(2*RPev/r, 4))*((y1-y2)/pow(r, 2))*3.71*exp(-pow(r/(2*RPev), 2));
}

void printArray(double arr[], int size) {
    printf("[ ");
    for (int i=0; i<size; i++) {
        printf("%f ", arr[i]);
    }
    printf("]\n");
}
void printArrayInt(int arr[], int size) {
    printf("[ ");
    for (int i=0; i<size; i++) {
        printf("%d ", arr[i]);
    }
    printf("]\n");
}

void printMatrix(int mat[], int n, int M){
    int i, j;
    for (i=0; i<N; i++) {
        for (j=0; j<M; j++) {
            printf("%d ", mat[i*M+j]);
        }
        printf("\n");
    }
}

void addArray(double* arr1, double* arr2, int size){
    for (int i=0; i<size; i++) {
        arr1[i] += arr2[i];
    }
}

void inithex2(double* x, double* y, double* q){
    srand((unsigned int) time(NULL));
    int i = 0; // to count number of particles placed
    int n = 0; // Check how many particles have been added to the system
    int pos[6]; //0 = not occupied, 1 = occupied. array with if neighbours in different positions exist or not
    // [right up, right down, down, left down, left up, up]
    x[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
    y[0]=((float) (-1+2*rand()) / (float)(RAND_MAX))*(SG/2);
    q[0]=((float) rand() / (float)(RAND_MAX))*2*PI;
    while (i<NP){
        checkNeighbour(x, y, i, pos, n);
        if (pos[0]==0) {
            x[n+1]=x[i]+RHEX*cs30;
            y[n+1]=y[i]+RHEX*sn30;
            q[n+1]=q[0];//((float) rand() / (float)(RAND_MAX))*2*PI;
            n+=1;
            if (n==NP-1) {
                break;
            }
        }
        if (pos[1]==0) {
            x[n+1]=x[i]+RHEX*cs30;
            y[n+1]=y[i]-RHEX*sn30;
            q[n+1]=q[0];//((float) rand() / (float)(RAND_MAX))*2*PI;
            n+=1;
            if (n==NP-1) {
                break;
            }
        }
        if (pos[2]==0) {
            x[n+1]=x[i];
            y[n+1]=y[i]-RHEX;
            q[n+1]=q[0];//((float) rand() / (float)(RAND_MAX))*2*PI;
            n+=1;
            if (n==NP-1) {
                break;
            }
        }
        if (pos[3]==0) {
            x[n+1]=x[i]-RHEX*cs30;
            y[n+1]=y[i]-RHEX*sn30;
            q[n+1]=q[0];//((float) rand() / (float)(RAND_MAX))*2*PI;
            n+=1;
            if (n==NP-1) {
                break;
            }
        }
        if (pos[4]==0) {
            x[n+1]=x[i]-RHEX*cs30;
            y[n+1]=y[i]+RHEX*sn30;
            q[n+1]=q[0];//((float) rand() / (float)(RAND_MAX))*2*PI;
            n+=1;
            if (n==NP-1) {
                break;
            }
        }
        if (pos[5]==0) {
            x[n+1]=x[i];
            y[n+1]=y[i]+RHEX;
            q[n+1]=q[0];//((float) rand() / (float)(RAND_MAX))*2*PI;
            n+=1;
            if (n==NP-1) {
                break;
            }
        }
        i+=1;
    }
}

//Chech neighbour in the hexagonal lattice
void checkNeighbour(double* x, double* y, int index, int* pos, int n){ //n is the number of particles placed in the system
    int j;
    pos[0]=0;
    pos[1]=0;
    pos[2]=0;
    pos[3]=0;
    pos[4]=0;
    pos[5]=0; //initialize to zero
    for (j=0; j<n; j++){
        if (x[j]-x[index]-RHEX*cs30>-0.00001 && y[j]-y[index]-RHEX*sn30>-0.00001){
            pos[0]=1;
        }
        if (x[j]-x[index]-RHEX*cs30>-0.00001 && y[j]-y[index]+RHEX*sn30<0.00001){
            pos[1]=1;
        }
        if ((x[j]==x[index]) && (y[j]-y[index]+RHEX<0.00001)){
            pos[2]=1;
        }
        if (x[j]-x[index]+RHEX*cs30<0.00001 && y[j]-y[index]+RHEX*sn30<0.00001){
            pos[3]=1;
        }
        if (x[j]-x[index]+RHEX*cs30<0.00001 && y[j]-y[index]-RHEX*sn30<-0.00001){
            pos[4]=1;
        }
        if (x[j]==x[index] && y[j]-y[index]-RHEX>-0.00001){
            pos[5]=1;
        }
    }
}
