import numpy as np



class Mesh:
    def __init__(selfm h, tau):
        self.h = h 
        self.tau = tau
    def __getitem__(self, i, j):
        return (self.h*i, self.tau*j)