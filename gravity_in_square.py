import matplotlib.pyplot as plt
import numpy as np


def sqEq(a, b, c):
    d = b**2 - 4 * a * c
    if d < 0:
        return False
    else:
        x1 = (-b + d ** 0.5)/ (2 * a)
        x2 = (-b - d ** 0.5)/ (2 * a)
        if x1 > 0.5 or x1< -0.5:
            if x2 > 0.5 or x2< -0.5:
                return False
            else:
                return [x2]
        
        if x2 > 0.5 or x2< -0.5:
            if x1 > 0.5 or x1< -0.5:
                return False
            else:
                return [x1]
            
        return [x1, x2]
    

def R(phi):
    if 3 * np.pi/4 > phi_0 > np.pi/4 or 7 * np.pi/4 > phi_0 > 5 * np.pi/4:
        r_ = abs(0.5 / np.sin(phi))
    else:
        r_ = abs(0.5 / np.cos(phi))
    return r_


def angle_from_xy(x, y):
    if x > 0 and y > 0:
        return np.arctan(y / x)
    elif (x < 0 and y > 0) or (x < 0 and y < 0):
        return np.pi + np.arctan(y / x)  
    elif x > 0 and y < 0:  
        return 2 * np.pi + np.arctan(y / x)


def ellipse_param(vx, vy, phi):
    r = R(phi)
    v2 = (vx**2 + vy**2)
    E = v2 / 2 - gamma / r
    L = -r * (v2 ** 0.5) * np.sin(phi - np.arctan(vy/vx))

    p = L ** 2 / gamma
    e = np.sqrt(1 + (2 * E * (L ** 2)) / (gamma ** 2))
    teta = phi + np.arccos(((p / r) - 1) / e)

    return (p, e, teta)


def clock_sort(X, phi, L):
    if X == []:
        exit("Dont find any collisions")
    A = []
    for x in X:
        A.append([angle_from_xy(*x), x])
    
    A.sort()
    for i, a in enumerate(A):
        if abs(a[0] - phi) < 1e-6:
            if L > 0:
                return A[i+1][1]
            else:
                return A[i-1][1]
    


E_0 = 1
gamma = 1
L_0 = 1
teta_0 = np.pi/4 + 4/10

p = lambda L, gamma: L ** 2 / gamma
e = lambda E, L, gamma: np.sqrt(1 + (2 * E * (L ** 2)) / (gamma ** 2))

def find_collision(vx, vy, phi):
    r = R(phi)
    v2 = (vx**2 + vy**2)
    L = -r * (v2 ** 0.5) * np.sin(phi - np.arctan(vy/vx))

    p, e, teta = ellipse_param(vx, vy, phi)
    a = p / (1 - e ** 2)
    b = np.sqrt(a * p)
    c = e * a
    teta = teta + np.pi
    # Ax^2 + 2Bxy + Cy^2 + Dx + Ey + F = 0
    A = (np.cos(teta) ** 2 / (a ** 2)) + (np.sin(teta) ** 2 / (b ** 2))
    C = (np.sin(teta) ** 2 / (a ** 2)) + (np.cos(teta) ** 2 / (b ** 2))
    B = -np.cos(teta) * np.sin(teta) * (c ** 2 / ((a ** 2) * (b ** 2)))
    D = -2 * np.cos(teta) * c / (a ** 2)
    E = -2 * np.sin(teta) * c / (a ** 2)
    F = (c**2 / a**2) - 1

    # x = 0.5: Cy^2 + (B + E)y + A/4 + D/2 + F = 0 => a_1y^2 + b_1y + c_1 = 0
    a_1 = C
    b_1 = B + E
    c_1 = A/4 + D/2 + F

    # x = -0.5: Cy^2 + (-B + E)y + A/4 - D/2 + F = 0 => a_2y^2 + b_2y + c_2 = 0
    a_2 = C
    b_2 = -B + E
    c_2 = A/4 - D/2 + F

    # y = 0.5: Ax^2 + (B + D)x + C/4 + E/2 + F = 0 => a_3x^2 + b_3x + c_3 = 0
    a_3 = A
    b_3 = B + D
    c_3 = C/4 + E/2 + F

    # y = -0.5: Ax^2 + (-B + D)x + C/4 - E/2 + F = 0 => a_4x^2 + b_4x + c_4 = 0
    a_4 = A
    b_4 = -B + D
    c_4 = C/4 - E/2 + F

    X = []
    tmp = sqEq(a_1, b_1, c_1)
    if tmp:
        X.append([0.5, tmp[0]])
        if len(tmp) == 2:
            X.append([0.5, tmp[1]])

    tmp = sqEq(a_2, b_2, c_2)
    if tmp:
        X.append([-0.5, tmp[0]])
        if len(tmp) == 2:
            X.append([-0.5, tmp[1]])
    
    tmp = sqEq(a_3, b_3, c_3)
    if tmp:
        X.append([tmp[0], 0.5])
        if len(tmp) == 2:
            X.append([tmp[1], 0.5])

    tmp = sqEq(a_4, b_4, c_4)
    if tmp:
        X.append([tmp[0], -0.5])
        if len(tmp) == 2:
            X.append([tmp[1], -0.5])
    x = clock_sort(X, phi, L)
    return x

def find_vel_before(E, x, y, vx, vy, phi):
    p, e, teta = ellipse_param(vx, vy, phi)

    r = R(phi)
    E = (vx**2 + vy**2) / 2 - gamma / r
    vel_before = np.sqrt(2 * (E + gamma / R(angle_from_xy(x, y))))

    



vx_0 = 1
vy_0 = 1
phi_0 = 1 * np.pi/4 + 4/10

r_0 = R(phi_0)

r = np.linspace(0, r_0, 100)

plt.plot([0.5, -0.5, -0.5, 0.5, 0.5], [0.5, 0.5, -0.5, -0.5, 0.5], c="tab:red")
plt.plot(r * np.cos(phi_0), r * np.sin(phi_0))
plt.plot([0], [0], 'o', c="tab:red")
plt.plot([r_0 * np.cos(phi_0), r_0 * np.cos(phi_0) + vx_0 * 0.1], [r_0 * np.sin(phi_0), r_0 * np.sin(phi_0) + vy_0 * 0.1])

angle = np.linspace(0, 2 * np.pi, 1000)
p, e, teta = ellipse_param(vx_0, vy_0, phi_0)
r_el = p / (1 + e * np.cos(angle - teta))
x = find_collision(vx_0, vy_0, phi_0)










plt.plot([x[0]], [x[1]], 'o', c="tab:red")
plt.plot(r_el * np.cos(angle), r_el * np.sin(angle))
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.show()
