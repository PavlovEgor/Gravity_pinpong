import matplotlib.pyplot as plt
import numpy as np

gamma = 1
m = 1
tol = 1e-6

def ellipse_param(vx, vy, phi):
    r = R(phi)
    Energy = m * (vx**2 + vy**2) / 2 - gamma * m / r
    Moment = m * r * ((vx**2 + vy**2) ** 0.5) * np.sin(np.arctan(vy/vx) - phi)
    
    p = Moment ** 2 / gamma
    e = np.sqrt(1 + (2 * Energy * (Moment ** 2)) / (gamma ** 2))

    theta = phi - np.arccos(((p / r) - 1) / e)
    if abs(tangent_el_eq(p, e, theta, R(phi) * np.cos(phi), R(phi) * np.sin(phi)) - np.arctan(vy/vx)) > tol:
        theta = phi + np.arccos(((p / r) - 1) / e)

    # print(tangent_el_eq(p, e, theta, R(phi) * np.cos(phi), R(phi) * np.sin(phi)), np.arctan(vy/vx))
    return (p, e, theta)


def general_eq_param(p, e, theta):
    a = p / (1 - e ** 2)
    b = np.sqrt(a * p)
    c = e * a

    # Ax^2 + 2Bxy + Cy^2 + Dx + Ey + F = 0
    A = (np.cos(theta) ** 2 / (a ** 2)) + (np.sin(theta) ** 2 / (b ** 2))
    C = (np.sin(theta) ** 2 / (a ** 2)) + (np.cos(theta) ** 2 / (b ** 2))
    B = -np.cos(theta) * np.sin(theta) * (c ** 2 / ((a ** 2) * (b ** 2)))
    D = 2 * np.cos(theta) * c / (a ** 2)
    E = 2 * np.sin(theta) * c / (a ** 2)
    F = (c**2 / a**2) - 1

    return (A, B, C, D, E, F)

def R(phi):
    return 1 / np.sin(phi)


def angle_from_xy(x, y):
    if x > 0 and y > 0:
        return np.arctan(y / x)
    elif (x < 0 and y > 0) or (x < 0 and y < 0):
        return np.pi + np.arctan(y / x)  
    elif x > 0 and y < 0:  
        return 2 * np.pi + np.arctan(y / x)
    

def tangent_el_eq(p, e, theta, x, y):
    A, B, C, D, E, F = general_eq_param(p, e, theta)
    return np.arctan(- (D + 2 * A * x + 2 * B * y) / (E + 2 * C * y + 2 * B * x))

def find_collision(vx, vy, phi):

    p, e, theta = ellipse_param(vx, vy, phi)
    A, B, C, D, E, F = general_eq_param(p, e, theta)
    

    # y = 1: Ax^2 + (2B + D)x + C + E + F = 0 => a_1x^2 + b_1x + c_1 = 0
    a_1 = A
    b_1 = 2 * B + D
    c_1 = C + E + F

    d = b_1**2 - 4 * a_1 * c_1
    x1 = (-b_1 + d ** 0.5) / (2 * a_1)
    x2 = (-b_1 - d ** 0.5) / (2 * a_1)

    a1 = angle_from_xy(x1, 1)

    if abs(a1 - phi) < tol:
        new_x = x2
        new_phi = angle_from_xy(x2, 1)
    else:
        new_x = x1
        new_phi = a1
    v = np.sqrt((vx**2 + vy**2) - 2 * gamma * ((1 / R(phi)) - (1 / R(new_phi))))
    alpha = tangent_el_eq(p, e, theta, new_x, 1)
    new_vx = -v * np.cos(alpha) * np.sign(alpha)
    new_vy = abs(v * np.sin(alpha))


    return (p, e, theta, new_vx, new_vy, new_phi)

vx_0 = 0.2
vy_0 = 0.2
phi_0 = 2.475859696
r_0 = R(phi_0)

plt.plot([-2, 2], [1, 1], c="tab:red")
plt.plot([0], [0], 'o', c="tab:red")

vx = vx_0
vy = vy_0
phi_old = phi_0
Phi = [phi_0]
for _ in range(10000):

    p, e, teta, vx, vy, phi_new = find_collision(vx, vy, phi_old)
    # plt.plot([R(phi_new) * np.cos(phi_new), R(phi_new) * np.cos(phi_new) + vx ], [R(phi_new) * np.sin(phi_new), R(phi_new) * np.sin(phi_new) + vy ])
    angle = np.linspace(phi_old, phi_new, 100)
    r_el = p / (1 + e * np.cos(angle - teta))
    plt.plot(r_el * np.cos(angle), r_el * np.sin(angle))

    phi_old = phi_new
    Phi.append(phi_new)


# plt.xlim(-2, 2)
# plt.ylim(-0.1, 2)
plt.show()

Phi = np.array(Phi)
X = 1 / np.tan(Phi)

plt.hist(X, bins=1000)
plt.show()
