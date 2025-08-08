import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


def hyperbolla(x, p, e):
    return p / (e * np.sin(x) - 1)


def circle(x, R, rho):
    return R - np.sqrt(rho * rho - x * x)


def parabolla(x, a, b, c, d):
    return b + a * x * x / 2 + c * x * x * x * x / 24 + d * (x ** 6) / 720

x = np.loadtxt("../data/x_finish.txt").T
plt.scatter(x[0][4700:], x[1][4700:], s=1)
plt.show()

x = np.loadtxt("../data/ellipses.txt").T

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.scatter(x[0], x[1], x[2], s=1)

ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[])

plt.show()

data = np.array([x[0][4710 * 100:], x[1][4710 * 100:], x[2][4710 * 100:]])
fig = plt.figure()
ax = fig.add_subplot(projection='3d')


def update(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])

N = x[0].shape[0]
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