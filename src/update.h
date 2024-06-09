#ifndef _UPDATE_H_
#define _UPDATE_H_

void update(double* x, double* y, double* q, double* fx, double* fy, double* tq, double* fxold, double* fyold, double* tqold, int*** cellList, int** numCell, int* firstRun, int NCx, int NCy, double* RPlist, int sample);

void forceUpdate(double* x, double* y, double* q, double* fx, double* fy, double* tq, int*** cellList, int** numCell, int NCx, int NCy, double* RPlist, int sample);

void updateFile(double* x, double* y, double* q, double* xtilde, double* ytilde, double* qtilde, double* fx, double* fy, double* tq, double* fxold, double* fyold, double* tqold, int* firstRun, double* RPlist, int** M, double** r);

void forceUpdateFile(double* x, double* y, double* q, double* fx, double* fy, double* tq, double* RPlist, int** M, double** r);

void findNCxiyi(int* NCxi, int* NCyi, double* x, double* y, int index);

double DV(double x1, double y1, double x2, double y2, double* RPlist, int index1, int index2);

double FHertz(double x1, double y1, double x2, double y2, double R1, double R2);

double FSpring(double x1, double y1, double x2, double y2, double r);

#endif