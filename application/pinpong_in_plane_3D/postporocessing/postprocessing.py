import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import shlex, subprocess


x = np.loadtxt("../data/x_finish.txt").T
plt.scatter(x[0], x[1], s=1)
plt.show()

x = np.loadtxt("../data/ellipses.txt").T
data = x
N = x[0].shape[0]
print(N)
for i in range(1, N//100, 2):
    data[2][i * 100:i * 100 +100] *= -1
    data[2][i * 100:i * 100 +100] += 2

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax.scatter(x[0], x[1], x[2], s=1)
ax.scatter(data[0], data[1], data[2], s=1)

ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[])

plt.show()

# data = np.array([x[0][4710 * 100:], x[1][4710 * 100:], x[2][4710 * 100:]])


fig = plt.figure()
ax = fig.add_subplot(projection='3d')


def update(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])


line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])

# Setting the axes properties
ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel('X')

ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, 2.0])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)
#ani.save('matplot003.gif', writer='imagemagick')
plt.show()