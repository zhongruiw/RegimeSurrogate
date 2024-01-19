from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from matplotlib.gridspec import GridSpec


def intL63_noisy(y, s, r, b, sgmx, sgmy, sgmz, dt, N, Nens):
    sol_save = np.zeros((3,Nens,N+1))
    sol_save[:,:,0] = y
    x0 = y[0,:]
    y0 = y[1,:]
    z0 = y[2,:]

    for i in range(N):
        x1 = x0 + (-s*x0 + s*y0)*dt + sgmx*np.random.normal(0,1,Nens)*dt**0.5
        y1 = y0 + (-x0*z0 + r*x0 - y0)*dt + sgmy*np.random.normal(0,1,Nens)*dt**0.5
        z1 = z0 + (x0*y0 - b*z0)*dt + sgmz*np.random.normal(0,1,Nens)*dt**0.5
        sol_save[0,:,i] = x1
        sol_save[1,:,i] = y1
        sol_save[2,:,i] = z1
        x0 = x1
        y0 = y1
        z0 = z1

    return sol_save


def intL63_noisy_fix(y, s, r, b, sgmx, sgmy, sgmz, dt, N, Nens, perts):
    sol_save = np.zeros((3,Nens,N+1))
    sol_save[:,:,0] = y
    x0 = y[0,:]
    y0 = y[1,:]
    z0 = y[2,:]

    for i in range(N):
        x1 = x0 + (-s*x0 + s*y0)*dt + sgmx*perts[0, i]*dt**0.5
        y1 = y0 + (-x0*z0 + r*x0 - y0)*dt + sgmy*perts[1, i]*dt**0.5
        z1 = z0 + (x0*y0 - b*z0)*dt + sgmz*perts[2, i]*dt**0.5
        sol_save[0,:,i] = x1
        sol_save[1,:,i] = y1
        sol_save[2,:,i] = z1
        x0 = x1
        y0 = y1
        z0 = z1

    return sol_save


if __name__ == '__main__':
    np.random.seed(2022)

    par = (10, 28, 8/3)
    sigma = (5, 5, 5)
    dt = 0.005
    N = 20000
    Nens = 2
    y0 = np.random.normal(0,1,(3,Nens))
    sol = intL63_noisy(y0, par[0], par[1], par[2], sigma[0], sigma[1], sigma[2], dt, N, Nens)

    # Plot
    sel0 = 10000; sel1 = 20000 # plot time range
    interv = 10 # plot interval
    xaxis = np.arange(sel0*dt, sel1*dt, interv*dt)

    fig = plt.figure()
    plt.subplots_adjust(wspace=0.2, hspace=0.5)     # Adjust the overall spacing of the figure
    gs0 = GridSpec(1, 2, figure=fig)

    gs00 = gs0[0].subgridspec(3, 6)
    ax1 = fig.add_subplot(gs00[0, :])
    ax2 = fig.add_subplot(gs00[1, :])
    ax3 = fig.add_subplot(gs00[2, :])

    ax1.plot(xaxis, sol[0,0,sel0:sel1:interv])
    ax1.set_xlim(50, 100)
    ax1.set_title('(a) Sample trajectory of x')

    ax2.plot(xaxis, sol[1,0,sel0:sel1:interv])
    ax2.set_xlim(50, 100)
    ax2.set_title('(b) Sample trajectory of y')

    ax3.plot(xaxis, sol[2,0,sel0:sel1:interv])
    ax3.set_xlim(50, 100)
    ax3.set_title('(c) Sample trajectory of z')
    ax3.set_xlabel('t')

    # create a new GridSpec for the second column
    gs01 = gs0[1].subgridspec(2, 7)

    ax4 = fig.add_subplot(gs01[0, :3])
    ax5 = fig.add_subplot(gs01[0, 4:])
    ax6 = fig.add_subplot(gs01[1, :3])
    ax7 = fig.add_subplot(gs01[1, 4:], projection='3d')

    ax4.plot(sol[0, 0, :], sol[1, 0, :], lw=0.5)
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.set_title('(d) 2D trajectory of x and y')

    ax5.plot(sol[1, 0, :], sol[2, 0, :], lw=0.5)
    ax5.set_xlabel('y')
    ax5.set_ylabel('z')
    ax5.set_title('(e) 2D trajectory of y and z')

    ax6.plot(sol[2, 0, :], sol[0, 0, :], lw=0.5)
    ax6.set_xlabel('z')
    ax6.set_ylabel('x')
    ax6.set_title('(f) 2D trajectory of z and x')

    ax7.plot(sol[0,0,:], sol[1,0,:], sol[2,0,:], lw=0.5)
    ax7.set_xlabel("x")
    ax7.set_ylabel("y")
    ax7.set_zlabel("z")
    ax7.set_title('(g) 3D trajectory of x,y and z')
    ax7.grid(False)

    plt.show()