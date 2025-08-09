import numpy as np
import subprocess

prom_path = "../build/pinpong_in_plane_3D"

Nvx = Nvy = Nvz = Nx = Ny = 3

VX0 = np.linspace(0.1, 0.4, Nvx)
VY0 = np.linspace(0.1, 0.4, Nvy)
VZ0 = np.linspace(0.1, 0.4, Nvz)


X0 = np.linspace(0.5, 2, Nx)
Y0 = np.linspace(0.5, 2, Ny)
Z0 = 1

NUM_OF_ELLIPSE = 20
NUM_OF_POINT = 2

R_minmax = np.empty(shape=(Nvx * Nvy * Nvz * Nx * Ny, 8))
i = 0
for vx in VX0:
    for vy in VY0:
        for vz in VZ0:
            for x in X0:
                for y in Y0:
                    subprocess.run([prom_path, 
                                    str(vx), 
                                    str(vy), 
                                    str(vz), 
                                    str(x), 
                                    str(y), 
                                    str(Z0), 
                                    str(NUM_OF_ELLIPSE), 
                                    str(NUM_OF_POINT)])
                    
                    R = np.loadtxt("../data/x_finish.txt").T[2]

                    R_minmax[i] = np.array([vx, vy, vz, x, y, Z0, np.min(R), np.max(R)])
                    i += 1


np.savetxt("init_param_to_R_minmax.txt", R_minmax)


