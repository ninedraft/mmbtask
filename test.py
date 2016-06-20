import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np


class Mesh:
    def __init__(self, h, tau):
        self.h = h 
        self.tau = tau
    def __getitem__(self, i, j):
        return (self.h*i, self.tau*j)

X = np.arange(0.0, 5.0, 0.1)
Y = np.arange(11.0, 20.0, 0.1)
print("building grid")
X, Y = np.meshgrid(X, Y)
Z = np.sin(X*Y - X)

print("building surface")
print('w shape: {0}'.format(Z.shape))    
print('x shape: {0}'.format(X.shape))
print('y shape: {0}'.format(Y.shape))

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()