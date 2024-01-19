import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz
from sklearn.cluster import KMeans
from L63_noisy import intL63_noisy, intL63_noisy_fix
from mpl_toolkits.mplot3d import Axes3D
from sklearn.metrics import silhouette_score


# fix the random seed
np.random.seed(2022)

# set experiment parameters
dt = 0.005  # model integration time step, dtda should be divisible to dt
N = 20000  # integration time steps

# set L63 model parameters
par = (10, 28, 8 / 3)  # sigma, rho, beta
sigma = (5, 5, 5)

# # create the nature run
# # X0 = np.random.normal(0, 1, (3, 1))
# # truth = intL63_noisy(X0, par[0], par[1], par[2], sigma[0], sigma[1], sigma[2], dt, N, 1)
# truth = np.load('nature_run1.npy')

# # -------------------------------------------------
# # K-means clustering
# # Create feature matrix (x, y, z) for clustering
# data_matrix = np.squeeze(truth).T
# x_data = data_matrix[:, 0]
# y_data = data_matrix[:, 1]
# z_data = data_matrix[:, 2]
# colors = ['b', 'orange', 'g', 'r']
#
# # Set up the loop and plot
# fig1, axes1 = plt.subplots(1, 3, figsize=(9, 3), subplot_kw=dict(projection='3d'))
#
# for ncenters, ax in enumerate(axes1.reshape(-1), 2):
#     kmeans = KMeans(n_clusters=ncenters)
#     kmeans.fit(data_matrix)
#     labels = kmeans.labels_
#
#     np.save('Kmeans_label_{0:d}center.npy'.format(ncenters), labels)
#
#     # Print the cluster centers
#     cluster_centers = kmeans.cluster_centers_
#     for i, center in enumerate(cluster_centers):
#         print(f"Cluster {i + 1} Center: {center}")
#
#     # The silhouette_score gives the average value for all the samples.
#     # This gives a perspective into the density and separation of the formed
#     # clusters
#     sil_avg = silhouette_score(data_matrix, labels)
#
#     # Plot assigned clusters, for each data point in training set
#     for j in range(ncenters):
#         ax.scatter(x_data[labels == j], y_data[labels == j], z_data[labels == j],
#                    label=f'Regime {j+1}', s=.5, color=colors[j])
#         # ax.scatter(cluster_centers[j, 0], cluster_centers[j, 1], cluster_centers[j, 2], marker='*', s=200, color='black',zorder=1)
#         ax.set_xlabel('X')
#         ax.set_ylabel('Y')
#         ax.set_zlabel('Z')
#
#     ax.set_title('Centers = {0}; Avg Sil = {1:.2f}'.format(ncenters, sil_avg))
#
# plt.subplots_adjust(wspace=0.25)  # Increased the width of intervals among subplots
# plt.show()

# # -------------------------------------------------
# # Fuzzy C-Means
# # Create feature matrix (x, y, z) for clustering
# data_matrix = np.squeeze(truth).T
#
# # Set up the loop and plot
# fig1, axes1 = plt.subplots(1, 3, figsize=(9, 3), subplot_kw=dict(projection='3d'))
# fpcs = []
#
# x_data = data_matrix[:, 0]
# y_data = data_matrix[:, 1]
# z_data = data_matrix[:, 2]
# colors = ['b', 'orange', 'g', 'r']
#
# for ncenters, ax in enumerate(axes1.reshape(-1), 2):
#     # # Perform Fuzzy C-means clustering
#     cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(
#         data_matrix.T, ncenters, 2, error=0.005, maxiter=1000, init=None)
#
#     # Store fpc values for later
#     fpcs.append(fpc)
#
#     # Plot assigned clusters, for each data point in training set
#     cluster_membership = np.argmax(u, axis=0)
#     np.save('FuzKmeans_label_{0:d}center.npy'.format(ncenters), cluster_membership)
#
#     for j in range(ncenters):
#         ax.scatter(x_data[cluster_membership == j], y_data[cluster_membership == j], z_data[cluster_membership == j],
#                    label=f'Regime {j+1}', s=.5, color=colors[j])
#         ax.set_xlabel('X')
#         ax.set_ylabel('Y')
#         ax.set_zlabel('Z')
#
#     # Mark the center of each fuzzy cluster
#     for pt in cntr:
#         ax.scatter(pt[0], pt[1], pt[2], marker='s', s=100, color='red')
#
#     ax.set_title('Centers = {0}; FPC = {1:.2f}'.format(ncenters, fpc))
#
# plt.subplots_adjust(wspace=0.25)  # Increased the width of intervals among subplots
# plt.show()

# # -------------------------------------------------
# label test data
# Create feature matrix (x, y, z) for clustering
truth = np.load('sol_test.npy')
data_matrix = np.squeeze(truth).T
x_data = data_matrix[:, 0]
y_data = data_matrix[:, 1]
z_data = data_matrix[:, 2]

# # # Perform Fuzzy C-means clustering
# ncenters = 2
# cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(
#     data_matrix.T, ncenters, 2, error=0.005, maxiter=1000, init=None)
#
# # save labeled data
# cluster_membership = np.argmax(u, axis=0)
# np.save('FuzKmeans_label_{0:d}center_test.npy'.format(ncenters), cluster_membership)
#
# fig1, axes1 = plt.subplots(1, 1, figsize=(3, 3), subplot_kw=dict(projection='3d'))
# colors = ['b', 'orange', 'g', 'r']
# ax = axes1
# for j in range(ncenters):
#     ax.scatter(x_data[cluster_membership == j], y_data[cluster_membership == j], z_data[cluster_membership == j],
#                label=f'Regime {j+1}', s=.5, color=colors[j])
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
# plt.show()

# # # -------------------------------------------------
# # classify by doubling time
# epsilon0 = 1e-10
# N_new = 1000
# doubling_times = []
# doubling_indexs = []
# for n in range(0, N + 1):
#     x0 = truth[0, :, n]
#     y0 = truth[1, :, n]
#     z0 = truth[2, :, n]
#
#     # Perturb the initial condition
#     perturbed_x0 = np.array([x0 + epsilon0, y0 + epsilon0, z0 + epsilon0])
#
#     # Generate random perturbations for integration
#     perts = np.random.normal(0, 1, [3, N_new])
#
#     # Integrate the Lorenz '63 system for both initial conditions
#     s0 = intL63_noisy_fix(truth[:, :, n], par[0], par[1], par[2], sigma[0], sigma[1], sigma[2], dt, N_new, 1, perts)
#     s1 = intL63_noisy_fix(perturbed_x0, par[0], par[1], par[2], sigma[0], sigma[1], sigma[2], dt, N_new, 1, perts)
#
#     # Calculate the Euclidean distance between the trajectories at each time step
#     distances = np.linalg.norm(s1 - s0, axis=1)
#
#     # Find the time index when the distance doubles
#     doubling_time_index = np.argmax(distances >= 2 * epsilon0)
#     doubling_indexs.append(doubling_time_index)
#     doubling_times.append(dt * doubling_time_index)
#
# np.save('doubling_times.npy', doubling_times)
# np.save('doubling_time_index.npy', doubling_time_index)
#
# # doubling_times = np.load('doubling_times.npy')
# # doubling_time_index = np.load('doubling_time_index.npy')
# # truth = np.load('nature_run1.npy')
#
# p_doubling = 0.25
# doubling_times = np.array(doubling_times)
# unstable = truth[:, 0, doubling_times <= p_doubling]
# stable = truth[:, 0, doubling_times > p_doubling]
#
# colors = ['b', 'r', 'g', 'orange']
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='3d')
#
# # Scatter plot for unstable points
# ax.scatter(unstable[0, :], unstable[1, :], unstable[2, :],
#            label=r'$\tau_2 > 0.25$', s=.5, color=colors[0])
#
# # Scatter plot for stable points
# ax.scatter(stable[0, :], stable[1, :], stable[2, :],
#            label=r'$\tau_2 \leq 0.25$', s=.5, color=colors[1])
#
# # Set labels and title
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title('Classification by doubling time')
#
# # Show the legend
# ax.legend()
# plt.show()

# # # -------------------------------------------------
# # label test data
# # classify by doubling time
# epsilon0 = 1e-10
# N_new = 1000
# doubling_times = []
# doubling_indexs = []
# for n in range(0, N + 1):
#     x0 = truth[0, :, n]
#     y0 = truth[1, :, n]
#     z0 = truth[2, :, n]
#
#     # Perturb the initial condition
#     perturbed_x0 = np.array([x0 + epsilon0, y0 + epsilon0, z0 + epsilon0])
#
#     # Generate random perturbations for integration
#     perts = np.random.normal(0, 1, [3, N_new])
#
#     # Integrate the Lorenz '63 system for both initial conditions
#     s0 = intL63_noisy_fix(truth[:, :, n], par[0], par[1], par[2], sigma[0], sigma[1], sigma[2], dt, N_new, 1, perts)
#     s1 = intL63_noisy_fix(perturbed_x0, par[0], par[1], par[2], sigma[0], sigma[1], sigma[2], dt, N_new, 1, perts)
#
#     # Calculate the Euclidean distance between the trajectories at each time step
#     distances = np.linalg.norm(s1 - s0, axis=1)
#
#     # Find the time index when the distance doubles
#     doubling_time_index = np.argmax(distances >= 2 * epsilon0)
#     doubling_indexs.append(doubling_time_index)
#     doubling_times.append(dt * doubling_time_index)
#
# np.save('doubling_times_test.npy', doubling_times)
# np.save('doubling_time_index_test.npy', doubling_time_index)
#
# p_doubling = 0.25
# doubling_times = np.array(doubling_times)
# unstable = truth[:, 0, doubling_times <= p_doubling]
# stable = truth[:, 0, doubling_times > p_doubling]
#
# colors = ['b', 'r', 'g', 'orange']
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='3d')
#
# # Scatter plot for unstable points
# ax.scatter(unstable[0, :], unstable[1, :], unstable[2, :],
#            label=r'$\tau_2 > 0.25$', s=.5, color=colors[0])
#
# # Scatter plot for stable points
# ax.scatter(stable[0, :], stable[1, :], stable[2, :],
#            label=r'$\tau_2 \leq 0.25$', s=.5, color=colors[1])
#
# # Set labels and title
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title('Classification by doubling time - test')
#
# # Show the legend
# ax.legend()
# plt.show()

doubling_times = np.load('doubling_times_test.npy')
doubling_time_index = np.load('doubling_time_index_test.npy')
label_doubling_time = np.zeros(np.shape(doubling_times)[0])
p_doubling = 0.25
label_doubling_time[doubling_times <= p_doubling] = 1

np.save('Doublingtime_label_test.npy', label_doubling_time)