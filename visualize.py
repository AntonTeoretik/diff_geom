from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas

points = pandas.read_csv('test.csv')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_ylim(-10,10)
ax.set_xlim(-10,10)

x = points['x'].values
y = points['y'].values
z = points['z'].values


ax.scatter(x, y, z, c='r', marker='o')

plt.show()