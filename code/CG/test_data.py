import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from matplotlib.gridspec import GridSpec


np.random.seed(2022)

def intL63_noisy(y, s, r, b, sgmx, sgmy, sgmz, dt, N, Nens):
    sol_save = np.zeros((3,Nens,N+1))
    der_save = np.zeros((3,Nens,N+1))
    sol_save[:,:,0] = y
    der_save[:,:,0] = y
    x0 = y[0,:]
    y0 = y[1,:]
    z0 = y[2,:]

    for i in range(N):
        der_save[0,:,i] = (-s*x0 + s*y0)*dt + sgmx*np.random.normal(0,1,Nens)*dt**0.5
        der_save[1,:,i] = (-x0*z0 + r*x0 - y0)*dt + sgmy*np.random.normal(0,1,Nens)*dt**0.5
        der_save[2,:,i] = (x0*y0 - b*z0)*dt + sgmz*np.random.normal(0,1,Nens)*dt**0.5
        x1 = x0 + der_save[0,:,i]
        y1 = y0 + der_save[1,:,i]
        z1 = z0 + der_save[2,:,i]
        sol_save[0,:,i] = x1
        sol_save[1,:,i] = y1
        sol_save[2,:,i] = z1

        x0 = x1
        y0 = y1
        z0 = z1

    return sol_save, der_save


par = (10, 28, 8/3)
sigma = (5, 5, 5)
dt = 0.005
N = 20000
Nens = 1
y0 = np.random.normal(1,1,(3,Nens))  # different IC center: (1,1) from training set's (0,1)
sol, der = intL63_noisy(y0, par[0], par[1], par[2], sigma[0], sigma[1], sigma[2], dt, N, Nens)

np.save('sol_test.npy', sol)
np.save('der_test.npy', der)