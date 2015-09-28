#include <stdio.h>
#include <math.h>

#define R_EDGE 1
#define M 100
#define N 100
#define r(i,j) (j*1/100)
#define SIGMA .4

int main(void){
    
    double y = 1;
    double dt = .1;

    for(double t = 0; t <= .2; t = t+dt){
       y = y + dt*( y + 2*t*exp(2*t) ); 
       printf("Y(%lf) = %lf\n", t, y);
    }
    
    return 0;
}
