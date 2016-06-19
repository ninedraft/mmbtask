from math import sqrt, pow, atan
h = 0.001
tau = 0.001
MaxX = 1
MinX = 0
NX = int((MaxX - MinX)/h)
MinT = 0
MaxT = 1
NT = int((MaxT - MinT)/tau)


class Point:
    def __init__(self, i, j, val):
        self.x = i*h 
        self.t = j*tau
        self.val = val 
     
w = [[0.0 for t in range(NT)] for x in range(NX)] 


def F(u):
    return atan(1 + pow(u, 4))

def dF(u):
    return (4 *pow(u, 3))/(1 + pow(1 + pow(u, 4), 2))

def ZeroTime(x):
    return pow(x, 2)

def ZeroX(t):
    return 0

def Projectile(i, j):
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
	    u = u_prev - equation(u_prev, i, j)/difeq (u_prev)
    return u 

for i in range(NX):
    w[i][0] = pow(i*h, 2)

for i in range(NT):
    w[0][i] = 0.0
