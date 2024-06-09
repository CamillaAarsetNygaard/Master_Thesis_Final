#ifndef _WALL_H_
#define _WALL_H_

void circWall(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void circWall2(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void circWall3(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void squareWall(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void squareWall2(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void rectWall(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void rectWall2(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void triWall(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void superEllipse(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

void superEllipse2(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist);

#endif