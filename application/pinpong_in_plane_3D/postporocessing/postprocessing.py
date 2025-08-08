import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.animation as animation


def hyperbolla(x, p, e):
    return p / (e * np.sin(x) - 1)


def circle(x, R, rho):
    return R - np.sqrt(rho * rho - x * x)


def parabolla(x, a, b, c, d):
    return b + a * x * x / 2 + c * x * x * x * x / 24 + d * (x ** 6) / 720

x = np.loadtxt("../data/x_finish.txt")
plt.hist(x, bins=2000, density=True, )
plt.ylim(0, 1)
plt.show()



# max_Rad = 1 / 0.4926357436
# plt.plot(max_Rad * np.cos(angle), max_Rad * np.sin(angle))
x = np.loadtxt("../data/ellipses.txt").T
plt.plot(x[0], x[1], 'o', markersize=1.)

plt.plot([-2, 2], [1, 1], c="tab:red")
plt.plot([0], [0], 'o', c="tab:red")
plt.show()






x = np.loadtxt("../data/ellipses.txt").T
N = x[0].shape[0]
print(N)

fig, ax = plt.subplots()

scat = ax.scatter(x[0][0], x[1][0], c="b", s=1)
ax.set(xlim=[-2, 2], ylim=[-0.01, 5])


def update(frame):
    # for each frame, update the data stored on each artist.
    x_ = x[0][:frame]
    y_ = x[1][:frame]
    # update the scatter plot:
    data = np.stack([x_, y_]).T
    scat.set_offsets(data)

    return (scat)


ani = animation.FuncAnimation(fig=fig, func=update, frames=N, interval=30)
plt.show()