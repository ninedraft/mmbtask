from math import sqrt, pow, atan, sin, cos, ceil, pi
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import progressbar

hx = 0.01
MinX = 0.0
MaxX = pi 
Nx = ceil((MaxX - MinX)/hx)

hy = 0.01
MinY = 0.0
MaxY = 2
Ny = ceil((MaxY - MinY)/hy)

tau = 0.01
MinT = 0.0
MaxT = 2
Nt = ceil((MaxT - MinT)/tau)

w = np.full((Nx + 1, Ny + 1, Nt + 1), 1.0)
middle = np.full((Nx + 1, Ny + 1), 1.0)

#init
for i in range(Nx + 1):
    for j in range(Ny + 1):
        w[i][j][0] = cos(hx*i)*sin(pi*hy*j)

def set_boundary(s):
    for j in range(0, Ny + 1):
        w[0][j][s] = 0
        w[Nx][j][s] = 0

    for i in range(0, Nx + 1):
	    w[i][0][s] = 0
	    w[i][Ny][s] = 0

def sweep(kapa, mu, Coeff, N):
	alpha = np.zeros((N + 1))
	alpha[1] = kapa[0] 
	beta = np.zeros((N + 1))
	beta[1] = mu[0]
	A = Coeff[0]
	B = Coeff[1]
	C = Coeff[2]
	U = Coeff[3]
	F = Coeff[4]
	for n in range(N):
	    alpha[n + 1] = (B[n] / (C[n]) - alpha[n] * A[n])
	    beta[n + 1] = (A[n] * beta[n] + F[n]) / (C[n] - alpha[n] * A[n])
	U[N] = (mu[1] + beta[N] * kapa[1]) / (1 - alpha[N] * kapa[1])
	for n in range(N - 1, 0, -1):
	    U[n] = alpha[n + 1] * U[n + 1] + beta[n + 1]

def calculate_layer(s):
	kapa = (0, 0)
	mu = (0, 0)
	A = np.zeros((Ny + 1))
	B = np.zeros((Ny + 1))
	C = np.zeros((Ny + 1))
	U = np.zeros((Ny + 1))
	F = np.zeros((Ny + 1))
	for j in range(1, Nx):
	    for i in range(1, Ny):
	        try:
	            A[i] = 0.5 * tau / (hx * hx)    
	            B[i] = 0.5 * tau / (hx * hx)
	            C[i] = 1.0 + tau / (hx * hx)
	            F[i] = (0.5 * tau / (hy * hy)) * w[i][j + 1][s - 1] + (1.0 - tau / (hy * hy)) * w[i][j][s - 1] + (0.5 * tau / (hy * hy)) * w[i][j - 1][s - 1] + (0.5 * tau) * sin(pi * hx * i) * (tau * (s - 0.5)) * (hy * j) * (hy * j - 1)
	        except Exception as err:
	            print('\ni = {0}, j = {1}, Nx = {2}, Ny = {3}, err = {4}'.format(i, j, Nx, Ny, err))
	            exit(1)
	    sweep(kapa, mu, (A, B, C, U, F), Ny)
	    for j in range(0, Ny + 1) :
	        w[i][j][s] = U[j];

def calculate_middle_layer(s):
	kapa = (0, 0)
	mu = (0, 0)
	A = np.zeros((Nx + 1))
	B = np.zeros((Nx + 1))
	C = np.zeros((Nx + 1))
	U = np.zeros((Nx + 1))
	F = np.zeros((Nx + 1))
	for j in range(1, Ny - 1):
	    for i in range(1, Nx - 1):
	        try:
	            A[i] = 0.5 * tau / (hx * hx)    
	            B[i] = 0.5 * tau / (hx * hx)
	            C[i] = 1.0 + tau / (hx * hx)
	            F[i] = (0.5 * tau / (hy * hy)) * w[i][j + 1][s - 1] + (1.0 - tau / (hy * hy)) * w[i][j][s - 1] + (0.5 * tau / (hy * hy)) * w[i][j - 1][s - 1] + (0.5 * tau) * sin(pi * hx * i) * (tau * (s - 0.5)) * (hy * j) * (hy * j - 1)
	        except Exception as err:
	            print('\ni = {0}, j = {1}, Nx = {2}, Ny = {3}, err = {4}'.format(i, j, Nx, Ny, err))
	            exit(1)
	    sweep(kapa, mu, (A, B, C, U, F), Nx)
	    for i in range(Nx + 1):
	        middle[i][j] = U[i];

def calculate():
    for s in range(1, Nt):
	    set_boundary(s)
	    calculate_layer(s)
	    calculate_middle_layer(s)

calculate()
