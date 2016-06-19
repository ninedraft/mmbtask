#include "stdio.h" 
#include <stdlib.h> 
#include "iostream" 
#include "clocale" 
#include "math.h" 
#include "iomanip" 
#include "fstream" 
#define pi 3.1415926535 
 
using namespace std; 
 
int Nx, Nt; 
 
double delta = 0.0001; 
double **U; 
double dx, dt; 
 
void base() { 
    U=(double**)malloc(Nt*sizeof(double*)); 
    for (int i=0; i<=Nt; i++) { 
        U[i]=(double*)malloc(Nx*sizeof(double)); 
    } 
} 
 
void bounds() { 
    for(int n=0; n<=Nx; n++) { 
        U[0][n]=sin(pi*(n*dx)); 
    }   
    for(int m=0; m<=Nt; m++) { 
        U[m][0]=exp(-m*dt)-1; 
    } 
} 
 
double f(double u) { 
    return (atan(exp(u))); 
} 
 
double f(int i, int j) { 
    return f(U[i][j]); 
} 
 
double df(double u) { 
    return exp(u)/(1+exp(2*u)); 
} 
 
double equation(double u, int i, int j) { 
    return ((u-U[i][j+1])/dt+(f(u)-f(U[i+1][j]))/dx); 
} 
 
double dif_eq (double u) { 
    return 1/dt + df(u)/dx; 
} 
 
double solution(int i, int j) { 
    double u_prev=1, u=0; 
    while (abs(u-u_prev)>delta) { 
        u_prev = u; 
        u = u_prev - equation(u_prev, i ,j)/dif_eq (u_prev); 
    } 
return u; 
} 
 
void maker() { 
    for(int i = 0; i < Nt-1; i++) { 
        for(int j = 0; j < Nx-1; j++) { 
            U[i+1][j+1]=solution(i, j); 
        } 
    } 
} 
 
int OMM1 (void) { 
    Nx = 100; Nt = 100; 
    dx = (1./Nx), dt = (1.8/Nt); 
    base();
    bounds(); 
    maker(); 
    return 0; 
} 
   
