from math import sqrt, pow, atan, sin, ceil
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import progressbar

h = 0.005
tau = 0.005
MaxX = 1.0
MinX = 0.0
NX = ceil((MaxX - MinX)/h)
MinT = 0.0
MaxT = 2.0
NT = ceil((MaxT - MinT)/tau)
delta = 0.01

print(NX, NT)

tick = time.clock()

w = np.zeros((NX, NT))

def F(u):
    return atan(1 + pow(u, 4))

def dF(u):
    return (4 *np.power(u, 3))/(1 + np.power(1 + np.power(u, 4), 2))

def ZeroTime(x):
    return pow(x, 2)

def ZeroX(t):
    return 0

def Mesh(i, j):
	return (i*h, j*tau, w[i][j])

def equation(u, i, j):
    uu = 0
    try:
        uu = ((u - w[i][j]) / (tau) + (F(u) - F(w[i - 1][j + 1])) /(h))
    except Exception as err:
        print('\ni = {0}, j = {1}, err = {2}'.format(i, j, err))
        exit(1)
    return uu

def diffeq(u):
    return 1/tau + dF(u)/h

def solution(i, j, delta):
    u_prev = 1
    u = 0 
    while abs(u - u_prev) > delta: 
	    u_prev = u
	    u = u_prev - equation(u_prev, i, j)/diffeq(u_prev)
    return u 

Xn = NX - 2
Tn = NT - 2

def calculate(delta):
    max_val = (Xn)*(Tn)
    bar = progressbar.ProgressBar(max_value = max_val)
    n = 1
    for i in range(1, Xn): 
         for j in range(0, Tn):
            w[i][j+1] = solution(i, j, delta)
            try:
                bar.update(value = n)
            except Exception as err:
                print('\nn = {0}, max_val = {1}, err = {2}'.format(n, max_val, err))
                exit(1)
            n = n + 1

for i in range(NX):
    w[i][0] = np.power(i*h, 2)

print("calculating...")
calculate(delta)
print("\n")

print("building X")
X = np.zeros((Xn+2))
for i in range(len(X)):
    X[i] = i*h

print("building Y")
Y = np.zeros((Tn+2))
for i in range(len(Y)):
    Y[i] = i*tau

print("building grid")
Y, X = np.meshgrid(Y, X)

print("building surface")
print('w shape: {0}'.format(w.shape))    
print('x shape: {0}'.format(X.shape))
print('y shape: {0}'.format(Y.shape))


fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, w, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)
ax.set_xlim(MinX - 0.01, MaxX + 0.01)
ax.set_ylim(MinT - 0.01, MaxT + 0.01)
ax.zaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
print(time.clock() - tick)
plt.show()