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

x = np.loadtxt("x_finish.txt")
plt.hist(x, bins=2000, density=True, )
plt.ylim(0, 1)
plt.show()


# p = 10
# e = 10
# r_el = -p / (1 - e * np.sin(angle))
# plt.plot(r_el * np.cos(angle), r_el * np.sin(angle))

x = np.loadtxt("max.txt").T
R = np.loadtxt("max_polar.txt").T
angle = np.linspace(np.max(R[1]), np.min(R[1]), 100)
xlin = np.linspace(np.max(x[0]), np.min(x[0]), 100)

plt.plot(x[0][:1000], x[1][:1000], 'o', markersize=1.)
popt, pcov = curve_fit(parabolla, x[0][:1000], x[1][:1000], p0=[0.05093547645933048, 1.1263674126786873, -0.0028627479070263504, 0.000001])
print(*popt)
# plt.plot(hyperbolla(angle, *popt) * np.cos(angle), hyperbolla(angle, *popt) * np.sin(angle))
plt.plot(xlin, parabolla(xlin, *popt))



max_Rad = 1 / 0.4926357436
plt.plot(max_Rad * np.cos(angle), max_Rad * np.sin(angle))
# x = np.loadtxt("ellipses.txt").T
# plt.plot(x[0], x[1], 'o', markersize=1.)

plt.plot([-2, 2], [1, 1], c="tab:red")
plt.plot([0], [0], 'o', c="tab:red")
plt.show()






x = np.loadtxt("ellipses.txt").T
N = x[0].shape[0]
print(N)

fig, ax = plt.subplots()

scat = ax.scatter(x[0][0], x[1][0], c="b", s=1)
ax.set(xlim=[-2, 2], ylim=[-0.01, 2])


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