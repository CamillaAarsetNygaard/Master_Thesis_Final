#ifndef _OLD_H_
#define _OLD_H_

double xd(double theta, const double Fp, const double F12x);

double yd(double theta, const double Fp, const double F12y);

double thetad(double theta, const double F12[], const double beta);

double euler(double* yn, double h, double fd[], double tn, int sizey);

void evol(double t0, int n, double Fp, double xyt10[], double xyt20[]);

void printArray2(double arr[], int size);

double Fljx2(double x1, double y1, double x2, double y2, double RPev);

double Fljy2(double x1, double y1, double x2, double y2, double RPev);

double Fljx(double x1, double y1, double x2, double y2, double RPev);

double Fljy(double x1, double y1, double x2, double y2, double RPev);

void printArrayInt(int arr[], int size);
void printMatrix(int mat[], int n, int M);
void printArray(double arr[], int size);

void addArray(double* arr1, double* arr2, int size);

void inithex2(double* x, double* y, double* q);

void checkNeighbour(double* x, double* y, int index, int* pos, int n); 

#endif