
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas
import numpy as np
from matplotlib import cm


points = pandas.read_csv('test.csv')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_ylim(-3,3)
ax.set_xlim(-3,3)
ax.set_zlim(-5,1)

x = points['x'].values
y = points['y'].values
z = points['z'].values


ax.plot(x, y, z, color="red", linewidth=3)

# surface
X = np.arange(-3, 3, 0.05)
Y = np.arange(-3, 3, 0.05)
X, Y = np.meshgrid(X, Y)
R = -1.0 / (X**2 + Y**2)

ax.plot_surface(X, Y, R, linewidth=0, antialiased=False, alpha=0.3, color = "green")


plt.show()
