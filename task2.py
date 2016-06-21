from math import sqrt, pow, atan, sin, ceil, pi
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

w = np.zeros((Nx, Ny))