#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "update.h"
#include "basics.h"
#include "initfnc.h"




int main() {
    // srand(17); //Seed for random number generator when wanting the same starting polarization
    int i, j, k;
    double *x, *y, *q, *fx, *fy, *tq, *fxold, *fyold, *tqold, *RPlist, *RPlistOld, *xtilde, *ytilde, *qtilde;
    int firstrun = 1; //1 = first run, 0 = not first run
    int **M; //Neighbor list for the particles
    int ***cellList; //The first and second dimension is the cell position, the third dimension is the particles in the cell. -1 means empty 
    int **numCell; //The number of particles in each cell
    double **r; //The ideal distance between the particles when using a file as the initial configuration
    int NCx, NCy; //Number of cells in cell list in the x and y direction
    double totAreaPart = 0.0; //The total area of the particles
    double areaWall = 0.0; //The area inside the wall
    FILE *fptr;
    FILE *fptrForce;
    FILE *inforptr;
    FILE *configptr;
    FILE *infoconfigptr;

    //Check if particles have enough space to be placed
    if (wall==0) {
        areaWall = PI*R*R;
        if (NP*PI*RP*RP > areaWall) {
            printf("The particles can not fit in the circle\n");
            return 0;
        }
    }
    else if (wall==1) {
        areaWall = 4*R*R;
        if (NP*PI*RP*RP > areaWall) {
            printf("The particles can not fit in the square\n");
            return 0;
        }
    }

    else if (wall==4) {
        areaWall = 4*a*b*(tgammaf(1+1/N))*(tgammaf(1+1/N))/(tgammaf(1+2/N));
        if (NP*PI*RP*RP > areaWall) {
            printf("The particles can not fit in the superellipse\n");
            return 0;
        }
    }

   
    fptr = fopen("Utvikling.txt","w"); 
    fptrForce = fopen("Kraft.txt","w");
    inforptr = fopen("Info.txt","w");

    double RPev = RP;

     // find number of cells in the grid, haven't made a case for the triangle wall 
     //+2 so that we have empty cells at the ends and do not need to check if a cell is at the boundary
    if ((wall==0) | (wall==1)){ //circle or square
        NCx = ((int) ceil( (2*R+2*DELTA)/DELTA))+2;
        NCy = NCx;
    }

    else if (wall==2) {
        NCx = ((int) ceil( (2*SQNP*RP+RP)/DELTA))+2;
        NCy = ((int)  ceil ((SQ3*(SQNP-1)*RP+2*RP)/DELTA))+2;
    }

    else if (wall==4){
        NCx = ((int) ceil( (2*a+2*DELTA)/DELTA))+2;
        NCy = ((int) ceil( (2*b+2*DELTA)/DELTA))+2;
    }
    
    //Allocate memory for the arrays
    x = malloc(NP*sizeof(double));
    y = malloc(NP*sizeof(double));
    q = malloc(NP*sizeof(double));
    xtilde = malloc(NP*sizeof(double));
    ytilde = malloc(NP*sizeof(double));
    qtilde = malloc(NP*sizeof(double));
    fx = malloc(NP*sizeof(double));
    fy = malloc(NP*sizeof(double));
    tq = malloc(NP*sizeof(double));
    fxold = malloc(NP*sizeof(double));
    fyold = malloc(NP*sizeof(double));
    tqold = malloc(NP*sizeof(double));
    RPlist = malloc(NP*sizeof(double));
    RPlistOld = malloc(NP*sizeof(double));
    
    // Allocate memory for the first dimension
    cellList = (int ***)malloc(NCx * sizeof(int **));
    numCell = (int **)malloc(NCx * sizeof(int *));
    M = (int **)malloc(NP*sizeof(int*));
    r = (double **)malloc(NP*sizeof(double*));
    for (i=0; i<NP; i++) {
        M[i] = (int *)malloc(9*sizeof(int));
        r[i] = (double *)malloc(9*sizeof(double));
    }
    
    // Allocate memory for the remaining dimensions
    for (i = 0; i < NCx; i++) {
        cellList[i] = (int **)malloc(NCy * sizeof(int *));
        numCell[i] = (int *)malloc(NCy * sizeof(int));
        for (j = 0; j < NCy; j++) {
            cellList[i][j] = (int *)malloc(20 * sizeof(int));
        }
    }

    makeRadiuses(RPlist, &totAreaPart);


    if (initMethod==0) {//random initial configuration
        initxyq(x, y, q, RPlist);
    }
    else if (initMethod==1) { // use config file
        //Config 3, 4 and 5 are close packed 50 particles within a circular wall
        //Config 6 and 7 is closed packed 49 particles within a square wall
        initfile(x, y, q, RPlist, initFile, tessConfigName, M, r);
        totAreaPart = 0;
        for (i=0; i<NP; i++) {
            totAreaPart += PI*RPlist[i]*RPlist[i];
        }
        RPev=RPlist[0];
    }
    else if (initMethod==2) {//hexagonal initial configuration with 
        inithex(x, y, q, RPlist); 
    }
    else if (initMethod==3) { //make a config file
        configptr = fopen(configName,"w");
        infoconfigptr = fopen(infoConfigName,"w");
        // initxyq(x, y, q, RPlist);
        initfile(x, y, q, RPlist, initFile, tessConfigName, M, r);
        totAreaPart = 0;
        for (i=0; i<NP; i++) {
            totAreaPart += PI*RPlist[i]*RPlist[i];
        }
        RPev=RPlist[0];
    }

    fprintf(inforptr,"%d %d %f %d %f %d %f %f %f %f %f\n", NP, R, DT, samples+1, XI, wall, EPS, RPev, a, b, N);

    for (j=0; j<NP; j++) {
            fprintf(fptr,"%.4f %.4f %.4f %.4f\n", x[j], y[j], q[j], RPlist[j]);
            fprintf(fptrForce,"%.4f %.4f\n", fx[j], fy[j]);
        }
    fprintf(fptr,"\n");
    fprintf(fptrForce,"\n");

    if (initMethod == 1) {
        for (i=0; i<samples; i++) {
            if (i%50==0) {
                printf("sample %d\n", i);}
            for (j=0; j<stepsInSample; j++) {
                updateFile(x, y, q, xtilde, ytilde, qtilde, fx, fy, tq, fxold, fyold, tqold, &firstrun, RPlist, M, r);
            }
            for (j=0; j<NP; j++) {
                fprintf(fptr,"%.4f %.4f %.4f %.4f\n", x[j], y[j], q[j], RPlist[j]);
                fprintf(fptrForce,"%.4f %.4f\n", fx[j], fy[j]);
            }
            if (i<samples-1) { //to avoid newline at end of file to make it easier to read into array in python
                fprintf(fptr,"\n");
                fprintf(fptrForce,"\n");
            }
        }
    }

    else{
    for (i=0; i<samples; i++) {
        if (i%50==0) {
            printf("sample %d\n", i);}
        if ((i==samples/2) && (reOrder==1) && (initMethod==3)) {
                    for (k=0; k<NP; k++) {
                        RPlistOld[k] = RPlist[k];
                    }
                    for (k=0; k<NP; k++) {
                        int l = rand() % NP;
                        if (RPlistOld[l]!=0){
                            RPlist[k] = RPlistOld[l];
                            RPlistOld[l] = 0;
                        }
                    }
                }
        for (j=0; j<stepsInSample; j++) {
            if ((totAreaPart > slowDownDen*areaWall && slowDown==1)){
                update(x, y, q, fx, fy, tq, fxold, fyold, tqold, cellList, numCell, &firstrun, NCx, NCy, RPlist, i*stepsInSample+j);
            }
            else {
                update(x, y, q, fx, fy, tq, fxold, fyold, tqold, cellList, numCell, &firstrun, NCx, NCy, RPlist, 0);
            }
            if (initMethod==3) {
                for (k=0; k<NP; k++) {
                    if (totAreaPart<density*areaWall){
                        RPlist[k] += growthRate*DT;
                        if (2*RPlist[k] > DELTA) {
                            printf("Particle %d is too big\n", k);
                            printf("Radius: %f\n", RPlist[k]);
                            printf("Area: %f\n", PI*RPlist[k]*RPlist[k]);
                            printf("Total area: %f\n", totAreaPart);
                            printf("Time: %f\n", i*DT);
                            return 0;
                        }
                        totAreaPart += growthRate*DT*PI*(2*RPlist[k]+3*DT*growthRate);
                    }
                }
            }
            
        }
        for (j=0; j<NP; j++) {
            fprintf(fptr,"%.4f %.4f %.4f %.4f\n", x[j], y[j], q[j], RPlist[j]);
            fprintf(fptrForce,"%.4f %.4f\n", fx[j], fy[j]);
            if (initMethod==3) {//make a config file
                if (i==samples-1){
                fprintf(configptr,"%.4f %.4f %.4f %.4f\n", x[j], y[j], q[j], RPlist[j]);}
            }
        }
        if (i<samples-1) { //to avoid newline at end of file to make it easier to read into array in python
            fprintf(fptr,"\n");
            fprintf(fptrForce,"\n");
        }
    }
    }
    
    if (initMethod==3) { //make a config file
        fprintf(infoconfigptr,"%d %d %f %f %f %f %f\n", NP, R, DT, RPev, a, b, N);
        fclose(configptr);
        fclose(infoconfigptr);
    }

    fclose(fptr);
    fclose(fptrForce);
    fclose(inforptr);

    return 0;
}
