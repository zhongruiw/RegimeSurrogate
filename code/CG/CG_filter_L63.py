"""
Using the conditional Gaussian framework to filter the noisy Lorenz 63 system
"""

from L63_noisy import intL63_noisy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde


# fix the random seed
np.random.seed(2022)

# set experiment parameters
dtda = 0.005  # time between analyses
dt   = 0.005  # model integration time step, dtda should be divisible to dt
Nassim = 20000  # number of assimilation times
N = int(Nassim * dtda / dt)  # integration time steps
# Nens = 1  # number of ensembles, ensembles are not supported for now

# set L63 model parameters
par = (10, 28, 8 / 3)  # sigma, rho, beta
sigma = (5, 5, 5)

# create the nature run
X0 = np.random.normal(0, 1, (3, 1))
truth = intL63_noisy(X0, par[0], par[1], par[2], sigma[0], sigma[1], sigma[2], dt, N, 1)

# assume x is observed, y and z are unobserved variables
xt = truth[0, 0, :]

# initialize a background
perts = np.array([0, np.random.normal(0,sigma[1]), np.random.normal(0,sigma[1])])[:, None]
Xb = X0 + perts

# initialize the uII mean and covariance series
u2_means = np.zeros((2, Nassim+1))
u2_covs = np.zeros((2, 2, Nassim+1))
u2_means[:, 0] = np.array([Xb[1, 0], Xb[2, 0]])
u2_covs[:, :, 0] = np.diag(np.array([sigma[1], sigma[2]]))
Sigma1 = np.array(sigma[0])[None, None]
Sigma2 = np.diag(np.array([sigma[1], sigma[2]]))

# filter the unobserved variables when an observation is available
u2_mean0 = u2_means[:, 0][:, None]
u2_cov0 = u2_covs[:, :, 0]
invBoB = np.linalg.inv(Sigma1 @ Sigma1.T)

for i in range(1, Nassim+1):
    x0 = xt[i-1]  # x at previous step
    x = xt[i]   # x at current step
    A0 = -par[0] * x0
    A1 = np.array([par[0], 0])[None, :]
    a0 = np.array([par[1], 0])[:, None] * x0
    a1 = np.array([[-1, -x0], [x0, -par[2]]])

    # Update posterior mean and covariance
    u2_mean = u2_mean0 + (a0 + a1 @ u2_mean0) * dt + (u2_cov0 @ A1.T) @ invBoB @ (
                x - x0 - A0 * dt - A1 @ u2_mean0 * dt)
    u2_cov = u2_cov0 + (a1 @ u2_cov0 + u2_cov0 @ a1.T + Sigma2 @ Sigma2.T - (u2_cov0 @ A1.T) @ invBoB @ (
                    u2_cov0 @ A1.T).T) * dt

    u2_mean0 = u2_mean
    u2_cov0 =u2_cov

    # save the time series
    u2_means[:, i] = u2_mean[:, 0]
    u2_covs[:, :, i] = u2_cov

# Plot
sel0 = 10000; sel1 = 20000 # plot time range
interv = 10 # plot interval
xaxis = np.arange(sel0*dt, sel1*dt, interv*dt)

fig = plt.figure()
widths = [8, 2]
heights = [2, 2, 2]
spec = fig.add_gridspec(ncols=2, nrows=3, width_ratios=widths, height_ratios=heights)

plt.subplots_adjust(wspace=0.2, hspace=0.5)     # Adjust the overall spacing of the figure
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[1, 0])
ax3 = fig.add_subplot(spec[2, 0])
ax4 = fig.add_subplot(spec[0, 1])
ax5 = fig.add_subplot(spec[1, 1])
ax6 = fig.add_subplot(spec[2, 1])

# plot time series
ax1.plot(xaxis, truth[0,0,sel0:sel1:interv])
ax1.set_xlim(sel0*dt, sel1*dt)
ax1.set_ylabel('x')
ax1.set_title('time series')

ax2.plot(xaxis, truth[1,0,sel0:sel1:interv], label='truth')
ax2.plot(xaxis, u2_means[0,sel0:sel1:interv], label='filter')
ax2.set_xlim(sel0*dt, sel1*dt)
ax2.set_ylabel('y')
ax2.legend()

ax3.plot(xaxis, truth[2,0,sel0:sel1:interv], label='truth')
ax3.plot(xaxis, u2_means[1,sel0:sel1:interv], label='filter')
ax3.set_xlim(sel0*dt, sel1*dt)
ax3.set_ylabel('z')
ax3.set_xlabel('t')

# plot pdf
samples = truth[0, 0, :]
kde = gaussian_kde(samples)
xticks = np.linspace(samples.min(), samples.max(), 100)
p = kde.evaluate(xticks)
ax4.plot(xticks, p)
ax4.set_title('PDF')

samples = truth[1, 0, :]
kde = gaussian_kde(samples)
xticks = np.linspace(samples.min(), samples.max(), 100)
p = kde.evaluate(xticks)
ax5.plot(xticks, p, label='truth')
samples = u2_means[0, :]
kde = gaussian_kde(samples)
p = kde.evaluate(xticks)
ax5.plot(xticks, p, label='filter')

samples = truth[2, 0, :]
kde = gaussian_kde(samples)
xticks = np.linspace(samples.min(), samples.max(), 100)
p = kde.evaluate(xticks)
ax6.plot(xticks, p, label='truth')
samples = u2_means[1, :]
kde = gaussian_kde(samples)
p = kde.evaluate(xticks)
ax6.plot(xticks, p, label='filter')

plt.show()