from math import sqrt, pow, atan
h = 0.001
tau = 0.001
MaxX = 1
MinX = 0
NX = int((MaxX - MinX)/h)
MinT = 0
MaxT = 1
NT = int((MaxT - MinT)/tau)
w = [[0.0 for t in range(NT)] for x in range(NX)] 

def F(u):
    return atan(1 + pow(u, 4))

def dF(u):
    return (4 *pow(u, 3))/(1 + pow(1 + pow(u, 4), 2))
