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

def progonka_neiman(A, B, C, F, v, n):
    alpha = np.zeros(n + 1)
    beta = np.zeros(n + 1)
    alpha[0] = 1.0
    beta[0] = 0.0
    for i in range(1, n):
        try:
            alpha[i] = C/(B - A*alpha[i-1])
            beta[i] = (A * beta[i-1] + F[i-1])/(B - A * alpha[i-1])
            v[n] = beta[n] / (1 - alpha[n])
        except Exception as err:
            print('n = {0} i = {1} alpha.shape = {2}, beta.shape = {3} v.shape = {4} err: {5}'.format(n, i, alpha.shape, beta.shape , v.shape, err))
            exit(1)
    for i in range(n - 1, 0, -1):
        v[i] = alpha[i] * v[i+1] + beta[i]
    v[0] = v[1]

def calculate():
    A1 = tau/(2*hx*hx)
    B1 = 1 + tau/(hx*hx)
    C1 = tau/(2*hx*hx)
    A2 = tau/(2*hy*hy)
    B2 = 1 + tau/(hy*hy)
    C2 = tau/(2*hy*hy)
    u = w
    temp = np.zeros((Nx + 1, Ny + 1))
    temp1 = np.zeros((Nx + 1))
    temp2 = np.zeros((Ny + 1))
    F1 = np.zeros((Nx - 1))
    F2 = np.zeros((Ny - 1))
    print('Nx = {0} Ny = {1}'.format(Nx, Ny))
    print("setting initial values...")
    bar = progressbar.ProgressBar(max_value = (Nx + 1)*(Ny + 1))
    n = 0
    for i in range(0, Nx + 1):
        for j in range(0, Ny + 1):
            u[i][j][0] =  cos(i*hx)*cos(j*hy*pi)
            try:
                n =  (i + 1)*(Nx + 1) + j
                bar.update(value = n)
            except Exception as err:
                print('\nn = {0}, max = {1} i ={2} j = {3}'.format(n, bar.max_value, i, j))
                exit(1)
    print("\n")


    for k in range(0, Nt):
        for i in range(0, Nx + 1):
            for j in range(0, Ny + 1):
                try:
	                temp[i][j] = w[i][j][k]
                except Exception as err:
                    print('i = {0} j = {1} k = {3} temp.shape = {4}, w.shape = {5} err = err'.format(i, j, k, temp.shape, w.shape, err))
                    exit(1)

        # first half step
        for j in range(1, Ny):
            for i in range(0, Nx + 1):
                try:
                    temp1[i] = temp[i][j]
                except Exception as err:
                    print('i = {0} j = {1} k = {3} temp1.shape = {4}, w.shape = {5} err = err'.format(i, j, k, temp4.shape, w.shape, err))
                    exit(1)
            for i in range(0, Nx - 1):
                try:
                    F1[i] = A2 * (w[i + 1][j - 1][k] + w[i + 1][j + 1][k]) + (1 - 2 * C2) * w[i + 1][j][k]
                except Exception as err:
                    print('i = {0} j = {1} k = {3} F1.shape = {4}, w.shape = {5} err = err'.format(i, j, k, F1.shape, w.shape, err))
                    exit(1)

            progonka_neiman(A1, B1, C1, F1, temp1, Nx)
            for i in range(0, Nx + 1):
                temp[i][j] = temp1[i]
        # second half step
        for i in range(1, Nx):
            for j in range(0, Ny + 1):
                temp2[j] = temp[i][j]
            for j in range(0, Ny - 1):
                F2[j] = A1 * (temp[i - 1][j + 1] + temp[i + 1][j + 1]) + (1 - 2 * C1) * temp[i][j + 1]
            for j in range(0, Ny + 1):
                w[i][j][k + 1] = temp2[j]
calculate()