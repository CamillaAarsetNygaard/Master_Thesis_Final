#include "wall.h"
#include "basics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void circWall(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist){
    double r = sqrtf(x[i]*x[i]+y[i]*y[i]);
    double f = -2*A*(r-(R-RPlist[i]))/r;
    if ((r+RPlist[i])<R) {//The interaction between the wall and the particles
        *fwallx=f*x[i]*exp((r-(R-RPlist[i]))/L);
        *fwally=f*y[i]*exp((r-(R-RPlist[i]))/L);
        }
    else{
    *fwallx=f*x[i];
    *fwally=f*y[i];
    }
}

void circWall2(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist){
    double r = sqrtf(x[i]*x[i]+y[i]*y[i]);
    if ((r+RPlist[i])>=R) {//The interaction between the wall and the particles
            *fwallx=-2*A*(r-(R-RPlist[i]))*x[i]/sqrtf(x[i]*x[i]+y[i]*y[i]);
            *fwally=-2*A*(r-(R-RPlist[i]))*y[i]/sqrtf(x[i]*x[i]+y[i]*y[i]);
        }
}

void circWall3(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist){
    double r = sqrtf(x[i]*x[i]+y[i]*y[i]);
    *fwallx=-2*A*(r-(R-RPlist[i]))*x[i]/sqrtf(x[i]*x[i]+y[i]*y[i]);
    *fwally=-2*A*(r-(R-RPlist[i]))*y[i]/sqrtf(x[i]*x[i]+y[i]*y[i]);
}

void squareWall(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist){
    if (x[i]+RPlist[i]>=R) {
        *fwallx=-2*A*(x[i]-(R-RPlist[i]));
    }
    else if ((x[i]+RPlist[i]<R) && (x[i]>=0)) {
        *fwallx=-2*A*(x[i]-(R-RPlist[i]))*(1.0/(1.0+expf(-(x[i]-50.0*R/51.0)/0.9)));//exp((x[i]*x[i]-(R-RPlist[i])*(R-RPlist[i]))/L);
    }
    if (x[i]-RPlist[i]<=-R) {
        *fwallx=-2*A*(x[i]-(-R+RPlist[i]));
    }
    else if ((x[i]-RPlist[i]>-R) && (x[i]<=0)) {
        *fwallx=-2*A*(x[i]-(-R+RPlist[i]))*(1.0/(1.0+expf((x[i]+50.0*R/51.0)/0.9)));//*exp((x[i]*x[i]-(-R+RPlist[i])*(-R+RPlist[i]))/L);
    }
    if (y[i]+RPlist[i]>=R) {
        *fwally=-2*A*(y[i]-(R-RPlist[i]));
    }
    else if ((y[i]+RPlist[i]<R) && (y[i]>=0)) {
        *fwally=-2*A*(y[i]-(R-RPlist[i]))*(1.0/(1.0+expf(-(y[i]-50.0*R/51.0)/0.9)));//exp((y[i]*y[i]-(R-RPlist[i])*(R-RPlist[i]))/L);
    }
    if (y[i]-RPlist[i]<=-R) {
        *fwally=-2*A*(y[i]-(-R+RPlist[i]));
    }
    else if ((y[i]-RPlist[i]>-R) && (y[i]<=0)){
        *fwally=-2*A*(y[i]-(-R+RPlist[i]))*(1.0/(1.0+expf((y[i]+50.0*R/51.0)/0.9)));//exp((y[i]*y[i]-(-R+RPlist[i])*(-R+RPlist[i]))/L);
    }
}

void squareWall2(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist){
    if (x[i]+RPlist[i]>=R) {
        *fwallx=-2*A*(x[i]-(R-RPlist[i]));
    }
    if (x[i]-RPlist[i]<=-R) {
        *fwallx=-2*A*(x[i]-(-R+RPlist[i]));
    }
    if (y[i]+RPlist[i]>=R) {
        *fwally=-2*A*(y[i]-(R-RPlist[i]));
    }
    if (y[i]-RPlist[i]<=-R) {
        *fwally=-2*A*(y[i]-(-R+RPlist[i]));
    }
}



void rectWall ( double * x , double * y , double * fwallx , double * fwally , int i , double* RPlist) {
    if ( x [ i ]+ RP > rectWidthHalf) {
        * fwallx = -2* A *( x [ i ] -( rectWidthHalf - RP) ) ;
        }
    if ( x [ i ] - RP + rectWidthHalf < 0) {
        * fwallx = -2* A *( x [ i ] +( rectWidthHalf - RP ) ) ;
        }
    if ( y [ i ]+ RP >rectHeightHalf) {
        * fwally = -2* A *( y [ i ] -( rectHeightHalf - RP ) ) ;
        }
    if ( y [ i ] - RP + rectHeightHalf < 0 ) {
        * fwally = -2* A *( y [ i ] +( rectHeightHalf - RP ) ) ;
        }
}

void rectWall2( double * x , double * y , double * fwallx , double * fwally , int i , double* RPlist) {
    * fwallx += -2* A *( x [ i ] -( rectWidthHalf - RP) ) ;
    * fwallx += -2* A *( x [ i ] +( rectWidthHalf - RP ) ) ;
    * fwally += -2* A *( y [ i ] -( rectHeightHalf - RP ) ) ;
    * fwally += -2* A *( y [ i ] +( rectHeightHalf - RP ) ) ;
}



void triWall(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist){
    if (y[i]-RPlist[i]<=-R*sn30){
        *fwally+=-2*A*(y[i]-(-R*sn30+RPlist[i]));
    }
    if (y[i]-((-(sn30+1))/(cs30))*x[i]>=R-RPlist[i]/sn30){
        *fwallx+=-2*A*(y[i]+((sn30+1)/(cs30))*x[i]-R+RPlist[i]/sn30)*(((sn30+1))/(cs30));
        *fwally+=-2*A*(y[i]+((sn30+1)/(cs30))*x[i]-R+RPlist[i]/sn30);
    }
    if (y[i]-((sn30+1)/(cs30))*x[i]>=R-RPlist[i]/sn30){
        *fwallx+=+2*A*(y[i]-((sn30+1)/(cs30))*x[i]-R+RPlist[i]/sn30)*(((sn30+1))/(cs30));
        *fwally+=-2*A*(y[i]-((sn30+1)/(cs30))*x[i]-R+RPlist[i]/sn30);
    }

}

void superEllipse(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist){
    double r = sqrtf(x[i]*x[i]+y[i]*y[i]);
    double theta = atan2(y[i],x[i]);
    double Rellipse = pow(pow((fabs(cosf(theta)/a)), N)+pow((fabs(sinf(theta)/b)), N), -1/N);
    double f = -2*A*(r-(Rellipse-RPlist[i]))/r;
    if ((r+RPlist[i])>=Rellipse) {//The interaction between the wall and the particles
            *fwallx=f*x[i];
            *fwally=f*y[i];
        }
    else{
        *fwallx=f*x[i]*exp((r-(Rellipse-RPlist[i]))/L);
        *fwally=f*y[i]*exp((r-(Rellipse-RPlist[i]))/L);
    }
}

void superEllipse2(double* x, double* y, double* fwallx, double* fwally, int i, double* RPlist){
    double r = sqrtf(x[i]*x[i]+y[i]*y[i]);
    double theta = atan2(y[i],x[i]);
    double Rellipse = pow(pow((fabs(cosf(theta)/a)), N)+pow((fabs(sinf(theta)/b)), N), -1/N);
    double f = -2*A*(r-(Rellipse-RPlist[i]))/r;
    if ((r+RPlist[i])>=Rellipse) {//The interaction between the wall and the particles
            *fwallx=f*x[i];
            *fwally=f*y[i];
        }
}