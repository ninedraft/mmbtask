from math import sqrt, pow, atan, sin, ceil
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import progressbar

h = 0.01
tau = 0.01
MaxX = 1.0
MinX = 0.0
NX = ceil((MaxX - MinX)/h)
MinT = 0.0
MaxT = 0.3
NT = ceil((MaxT - MinT)/tau)

print(NX, NT)

tick = time.clock()

w = np.zeros((NX, NT))

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
    u = 0
    try:
        uu = ((u - w[i][j + 1]) / tau + (F(u) - F(w[i - 1][j + 1])) / h)
    except Exception as err:
        print('\ni = {0}, j = {1}, err = {2}'.format(i, j, err))
        exit(1)
    return u

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
    X = NX - 2
    T = NT - 2
    bar = progressbar.ProgressBar(max_value=X*T - X - T)
    n = 0
    for i in range(1, X): 
         for j in range(1, T):
            w[i][j+1] = solution(i, j, delta)
            bar.update(value = n)
            n = n + 1

for i in range(NX):
    w[i][0] = np.power(i*h, 2)

print("calculating...")
calculate(0.01)
print("calculated")

fig = plt.figure()
ax = fig.gca(projection='3d')
print("building X")
X = np.arange(MinX, MaxX, h)[:-2]
print("building Y")
Y = np.arange(MinT, MaxT, tau)[:-2]
print("building grid")
X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
print("building surface")    
surf = ax.plot_surface(X, Y, w, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
print( time.clock() - tick)
plt.show()