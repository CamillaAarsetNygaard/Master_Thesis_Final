#ifndef _BASICS_H_
#define _BASICS_H_

#define PI 3.14159265
#define NP 1000 //number of particles, OBS: must be a square number for hex config
#define FP 1.0 //force on each particle
#define DT 0.0001 //time step
#define SG 10000 //the size of the grid of where the particles can be placed inside
#define EPS 50.0 //Potential well depth lennard jones
#define C 100 //constant in the hertz contact elastic force
#define K 300.0 //constant in the spring interaction force
#define K4 0//K/(0.1*0.1) //second constant in the spring interaction force
#define L 0.5//0.95 //constant in exponential in wall force
#define XI 15.0//rotational dampning coefficient. OBS! This must be a double
#define R 50//157 //radius of the wall/2xlength of the square wall/size from origo to corner of triangle 160 = 10000 particles
#define a 50.0 //parameter in superellipse, semi-major axis
#define b 50.0 //parameter in superellipse, semi-minor axis
#define N 4.0 //parameter in superellipse
#define A 40.0 //strength of the wall force/potential
#define RP 1.1 //radius of the particles, OBS: it has to be double
#define samples 1000 //number of samples written to file
#define stepsInSample 3000//number of steps of code run between each sample
#define totalSamples samples*stepsInSample //Total number of samples
#define DELTA 4*RP //size of cell side in cell list OBS: must be larger than 2*RP and the cutoff of the potential
#define DIFF 0.0 //thing in noise
#define wall 0 //0 = circle, 1 = square, 2 = rectangle,  3 = triangle, 4 = superellipse
#define method  0//0 = Euler, 1 = Adams-Bashforth, 2 = RK2
#define initMethod 1 //0 = random, 1 = file, 2 = hex, 3 = make config file
#define diffRadius 1 //0 = same radius, 1 = different radius, OBS: DO NOT USE DIFFERENT RADIUS FOR HEX CONFIG
#define slowDown  0//0 = no slowdown, 1 = slowdown
#define slowDownDen 0 //density at which the slowdown starts
#define radDev 0.1 //deviation in radius as in RP+-radDev*RP
#define reOrder 0 //0 = no reordering, 1 = reordering
#define growthRate 0.0027 //rate of growth of the particles
#define density 0.91//density of the particles
#define initFile "startconfigs/Config1.txt" //name of the file to read from
#define tessConfigName "startconfigs/tesselationConfig1.txt" //name of the file with the tesselation
#define configName "startconfigs/Config16.txt" //name of the file to write to
#define infoConfigName "startconfigs/InfoConfig16.txt" //name of the file to write to
#define cs30 cos(PI/6)
#define sn30 sin(PI/6)
#define SQNP sqrt(NP)
#define SQ3 sqrt(3)
#define rectWidth 2.0*SQNP*RP+RP
#define rectWidthHalf SQNP*RP + RP*0.5
#define rectHeight sqrt(3.0)*(sqrt(NP)-1.0)*RP+2.0*RP
#define rectHeightHalf (sqrt(3.0)*(sqrt(NP)-1.0)*RP+2.0*RP)*0.5
#define CUTOFF 2.9*RP //cutoff for the potential
#define CUTOFF2 2.9*RP*2.9*RP //cutoff squared
#define RHEX 2*RP //distance between particles in old hex config

int checkCollision(double* x, double* y, int index, double* RPlist); //function to check if the particles collide

void resetCell(int*** cellList, int** numCell, int NCx, int NCy); //function to reset the cellList so we can fill it up again

void placeInCell(int*** cellList, int** numCell, double* x, double* y, int NCx, int NCy); //function to place the particles in the cellList

void makeRadiuses(double* RPlist, double* totAreaPart); //function to make the radiuses of the particles

#endif