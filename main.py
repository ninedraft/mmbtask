from math import sqrt, pow, atan
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

h = 0.001
tau = 0.001
MaxX = 1
MinX = 0
NX = int((MaxX - MinX)/h)
MinT = 0
MaxT = 1
NT = int((MaxT - MinT)/tau)

tick = time.clock()

w = [[0.0 for t in range(NT)] for x in range(NX)] 


def F(u):
    return atan(1 + pow(u, 4))

def dF(u):
    return (4 *pow(u, 3))/(1 + pow(1 + pow(u, 4), 2))

def ZeroTime(x):
    return pow(x, 2)

def ZeroX(t):
    return 0

def Mesh(i, j):
	return (i*h, j*tau, w[i][j])

def equation(u, i, j):
    return ((u - w[i][j + 1]) / tau + (F(u) - F(w[i + 1][j])) / h)

def diffeq(u):
    return 1/tau + dF(u)/h

def solution(i, j, delta):
    u_prev = 1
    u = 0 
    while abs(u - u_prev) > delta: 
	    u_prev = u
	    u = u_prev - equation(u_prev, i, j)/diffeq(u_prev)
    return u 

def calculate(delta):
    for i in range(NX - 1): 
         for j in range(NT - 1):
            w[i+1][j+1] = solution(i, j, delta) 

for i in range(NX):
    w[i][0] = pow(i*h, 2)

for i in range(NT):
    w[0][i] = 0.0

calculate(0.01)

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
print( time.clock() - tick)
plt.show()