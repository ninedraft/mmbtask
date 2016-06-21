
#include "stdio.h"
#include "cmath"
#include "stdlib.h"
#include "iostream"
#include "clocale"
#include "iomanip"
#include "fstream"
#define pi 3.1415926535
using namespace std;

const double X = 1, Y = 1, T = 20;
double *A, *B, *C, *F, *U, *alpha, *beta;
double **Umiddle;
double ***u;
double dx, dy, tau;
int Nx, Ny, Nt;

void base() {
  u = (double ***)malloc((Nx + 1) * sizeof(double **));
  for (int i = 0; i <= Nx; i++) {
    u[i] = (double **)malloc((Ny + 1) * sizeof(double *));
    for (int j = 0; j <= Ny; j++)
      u[i][j] = (double *)malloc((Nt + 1) * sizeof(double));
  }
  Umiddle = (double **)malloc((Nx + 1) * sizeof(double *));
  for (int i = 0; i <= Nx; i++)
    Umiddle[i] = (double *)malloc((Ny + 1) * sizeof(double));
  cout « "data field has created" « endl;
}

void delete_base(double **L) {
  for (int i = 0; i <= Nx; i++)
    free(L[i]);
  free(L);
}

void delete_base_all() {
  for (int i = 0; i <= Nx; i++) {
    for (int j = 0; j <= Ny; j++)
      free(u[i][j]);
    free(u[i]);
  }
  delete_base(Umiddle);
}

void clear_base_all() {
  for (int i = 0; i <= Nx; i++)
    for (int j = 0; j <= Ny; j++) {
      Umiddle[i][j] = 1.0;
      for (int s = 0; s <= Nt; s++)
        u[i][j][s] = 1.0;
    }
}

void base_coeff(int N) {
  A = (double *)malloc((N + 1) * sizeof(double));
  B = (double *)malloc((N + 1) * sizeof(double));
  C = (double *)malloc((N + 1) * sizeof(double));
  F = (double *)malloc((N + 1) * sizeof(double));
  U = (double *)malloc((N + 1) * sizeof(double));

  alpha = (double *)malloc((N + 1) * sizeof(double));
  beta = (double *)malloc((N + 1) * sizeof(double));
}

void delete_base_coeff_all(int N) {
  free(A);
  free(B);
  free(C);
  free(F);
  free(U);
  free(alpha);
  free(beta);
}

void clear_base_coeff(double *M, int N) {
  for (int n = 0; n <= N; n++)
    M[n] = 0;
}

void clear_base_coeff_all(int N) {
  clear_base_coeff(A, N);
  clear_base_coeff(B, N);
  clear_base_coeff(C, N);
  clear_base_coeff(F, N);
  clear_base_coeff(U, N);
  clear_base_coeff(alpha, N);
  clear_base_coeff(beta, N);
}

void initial() {
  for (int i = 0; i <= Nx; i++)
    for (int j = 0; j <= Ny; j++)
      u[i][j][0] = 0;
}

void bounds(int s) {
  for (int j = 0; j <= Ny; j++) {
    u[0][j][s] = 0;
    u[Nx][j][s] = 0;
  }
  for (int i = 0; i <= Nx; i++) {
    u[i][0][s] = 0;
    u[i][Ny][s] = 0;
  }
}

void sweep(double kapa1, double kapa2, double mu1, double mu2, int N) {
  alpha[1] = kapa1;
  beta[1] = mu1;

  for (int n = 1; n <= N - 1; n++) {
    alpha[n + 1] = B[n] / (C[n] - alpha[n] * A[n]);
    beta[n + 1] = (A[n] * beta[n] + F[n]) / (C[n] - alpha[n] * A[n]);
  }

  U[N] = (mu2 + beta[N] * kapa2) / (1 - alpha[N] * kapa2);

  for (int n = N - 1; n >= 0; n--)
    U[n] = alpha[n + 1] * U[n + 1] + beta[n + 1];
}

void half_layer_maker(int s) {
  double kapa1 = 0;
  double kapa2 = 0;
  double mu1 = 0;
  double mu2 = 0;
  base_coeff(Nx);
  for (int j = 1; j <= Ny - 1; j++) {
    for (int i = 1; i <= Nx - 1; i++) {
      A[i] = 0.5 * tau / (dx * dx);
      B[i] = 0.5 * tau / (dx * dx);
      C[i] = 1.0 + tau / (dx * dx);
      F[i] = (0.5 * tau / (dy * dy)) * u[i][j + 1][s - 1] +
             (1.0 - tau / (dy * dy)) * u[i][j][s - 1] +
             (0.5 * tau / (dy * dy)) * u[i][j - 1][s - 1] +
             (0.5 * tau) * sin(pi * dx * i) * (tau * (s - 0.5)) * (dy * j) *
                 (dy * j - 1);
    }
    sweep(kapa1, kapa2, mu1, mu2, Nx);

    for (int i = 0; i <= Nx; i++)
      Umiddle[i][j] = U[i];
  }
  delete_base_coeff_all(Nx);
}

void layer_maker(int s) {
  double kapa1 = 1.0;
  double kapa2 = 1.0;
  double mu1 = 0;
  double mu2 = 0;
  base_coeff(Ny);
  for (int i = 1; i <= Nx - 1; i++) {
    for (int j = 1; j <= Ny - 1; j++) {
      A[j] = 0.5 * tau / (dy * dy);
      B[j] = 0.5 * tau / (dy * dy);
      C[j] = 1.0 + tau / (dy * dy);
      F[j] = (0.5 * tau / (dx * dx)) * Umiddle[i + 1][j] +
             (1.0 - tau / (dx * dx)) * Umiddle[i][j] +
             (0.5 * tau / (dx * dx)) * Umiddle[i - 1][j] +
             (0.5 * tau) * sin(pi * dx * i) * (tau * (s - 0.5)) * (dy * j) *
                 (dy * j - 1);
    }
    sweep(kapa1, kapa2, mu1, mu2, Ny);
    for (int j = 0; j <= Ny; j++)
      u[i][j][s] = U[j];
  }
  delete_base_coeff_all(Ny);
}

void maker() {
  for (int s = 1; s <= Nt; s++) {
    bounds(s);
    half_layer_maker(s);
    layer_maker(s);
  }
}

int OMM2() {
  int r, x, y;
  Nx = 100;
  Ny = 100;
  Nt = 100;
  tau = T / Nt;
  dx = X / Nx;
  dy = Y / Ny;
  base();
  clear_base_all();
  initial();
  maker();
  delete_base_all();
  return (0);
}