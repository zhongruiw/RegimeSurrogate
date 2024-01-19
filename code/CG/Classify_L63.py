from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from sklearn.model_selection import cross_validate
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score
from FKNN import FuzzyKNN, calculate_distances
import joblib


# Define a custom function to evaluate the predictions with multiple metrics
def evaluate_predictions(truth_label, predicted_label):
    accuracy = accuracy_score(truth_label, predicted_label)
    precision = precision_score(truth_label, predicted_label, average='macro')
    recall = recall_score(truth_label, predicted_label, average='macro')
    return accuracy, precision, recall,

np.random.seed(2022)

# # Load L63 time series data X and labels y
# # X should be a 2D array with shape (n_samples, n_features)
# sol = np.squeeze(np.load('nature_run1.npy'))
# der = np.squeeze(np.load('derivative.npy'))
# X = np.concatenate([sol, der], axis=0).T
# # y should be a 1D array with shape (n_samples,)
# y = np.load('FuzKmeans_label_2center.npy')

# # ------------------------------------------------------------
# # train KNN
# # Create a KNN classifier object
# k_neighbors = 5  # You can choose the number of neighbors (k)
# knn_classifier = KNeighborsClassifier(n_neighbors=k_neighbors)
#
# # Perform k-fold cross-validation and get the accuracy scores
# k_fold = 5  # Number of folds for cross-validation
# cv_results = cross_validate(knn_classifier, X, y, cv=k_fold,
#                         scoring=('accuracy', 'precision', 'recall', 'roc_auc'),
#                         return_estimator=True)
#                         # ,return_indices=True)
#
# # Calculate the mean accuracy and standard deviation of accuracy scores
# mean_accuracy = cv_results['test_accuracy'].mean()
# mean_precision = cv_results['test_precision'].mean()
# mean_recall = cv_results['test_recall'].mean()
# mean_rocauc = cv_results['test_roc_auc'].mean()
# print(f"Mean accuracy: {mean_accuracy:.2f}")
# print(f"Mean precision: {mean_precision:.2f}")
# print(f"Mean recall: {mean_recall:.2f}")
# print(f"Mean ROC_AUC: {mean_rocauc:.2f}")
#
# # # Save the trained model to a file
# # model_filename = './trained_model/knn_model.joblib'
# # joblib.dump(knn_classifier, model_filename)
#
# # ------------------------------------------------------------
# # train Fuzzy KNN
# k_neighbors = 5
# fknn_classifier = FuzzyKNN(k=k_neighbors)
#
# # # Perform k-fold cross-validation and get the accuracy scores
# # k_fold = 5  # Number of folds for cross-validation
# # cv_results = cross_validate(fknn_classifier, X, y, cv=k_fold,
# #                         scoring=('accuracy', 'precision', 'recall', 'roc_auc'),
# #                         return_estimator=True)
# #                         # ,return_indices=True)
# #
# # # Calculate the mean accuracy and standard deviation of accuracy scores
# # mean_accuracy = cv_results['test_accuracy'].mean()
# # mean_precision = cv_results['test_precision'].mean()
# # mean_recall = cv_results['test_recall'].mean()
# # mean_rocauc = cv_results['test_roc_auc'].mean()
# # print(f"Mean accuracy: {mean_accuracy:.2f}")
# # print(f"Mean precision: {mean_precision:.2f}")
# # print(f"Mean recall: {mean_recall:.2f}")
# # print(f"Mean ROC_AUC: {mean_rocauc:.2f}")
#
# # # Save the trained model to a file
# # model_filename = './trained_model/fknn_model.joblib'
# # joblib.dump(fknn_classifier, model_filename)
#
# # ------------------------------------------------------------
# # scatter truth and predict on test data
# # use Fuzzy C-means clusetered results as truth
# truth_label = np.load('FuzKmeans_label_2center_test.npy')
# sol = np.squeeze(np.load('sol_test.npy'))
# der = np.squeeze(np.load('der_test.npy'))
# test_data = np.concatenate([sol, der], axis=0).T
#
# # KNN train and predict
# knn_classifier.fit(X, y)
# knn_label = knn_classifier.predict(test_data)
# # reverse label 0 - 1
# knn_label = - knn_label + 1
# accuracy_knn, precision_knn, recall_knn = evaluate_predictions(truth_label, knn_label)
# roc_auc_score_knn = roc_auc_score(truth_label, 1 - knn_classifier.predict_proba(test_data)[:, 1], average='macro')
#
# print("Accuracy KNN:", accuracy_knn)
# print("Precision KNN:", precision_knn)
# print("Recall KNN:", recall_knn)
# print("ROC_AUC KNN:", roc_auc_score_knn)
#
# # FKNN train and predict
# fknn_classifier.fit(X, y)
# fknn_label = fknn_classifier.predict(test_data)
# # reverse label 0 - 1
# fknn_label = - fknn_label + 1
# accuracy_knn, precision_knn, recall_knn = evaluate_predictions(truth_label, fknn_label)
# roc_auc_score_knn = roc_auc_score(truth_label, 1 - knn_classifier.predict_proba(test_data)[:, 1], average='macro')
#
# print("Accuracy FKNN:", accuracy_knn)
# print("Precision FKNN:", precision_knn)
# print("Recall FKNN:", recall_knn)
# print("ROC_AUC FKNN:", roc_auc_score_knn)
#
# fig1, axes1 = plt.subplots(1, 3, figsize=(9, 3), subplot_kw=dict(projection='3d'))
# colors = ['b', 'r', 'g', 'orange']
# ncenters = 2
# x_data = sol[0, :]
# y_data = sol[1, :]
# z_data = sol[2, :]
#
# for j in range(ncenters):
#     axes1[0].scatter(x_data[truth_label == j], y_data[truth_label == j], z_data[truth_label == j],
#                label=f'Regime {j+1}', s=.5, color=colors[j])
# axes1[0].set_xlabel('X')
# axes1[0].set_ylabel('Y')
# axes1[0].set_zlabel('Z')
# axes1[0].set_title('Truth')
#
# for j in range(ncenters):
#     axes1[1].scatter(x_data[knn_label == j], y_data[knn_label == j], z_data[knn_label == j],
#                label=f'Regime {j+1}', s=.5, color=colors[j])
# axes1[1].set_xlabel('X')
# axes1[1].set_ylabel('Y')
# axes1[1].set_zlabel('Z')
# axes1[1].set_title('KNN')
#
# for j in range(ncenters):
#     axes1[2].scatter(x_data[fknn_label == j], y_data[fknn_label == j], z_data[fknn_label == j],
#                label=f'Regime {j+1}', s=.5, color=colors[j])
# axes1[2].set_xlabel('X')
# axes1[2].set_ylabel('Y')
# axes1[2].set_zlabel('Z')
# axes1[2].set_title('Fuzzy KNN')
#
# plt.subplots_adjust(wspace=0.25)  # Increased the width of intervals among subplots
# plt.savefig('Figure_6.png', dpi=150)
# plt.show()

# ------------------------------------------------------------
# ------------------------------------------------------------
# classify based on error growth
# Load L63 time series data X and labels y
# X should be a 2D array with shape (n_samples, n_features)
sol = np.squeeze(np.load('nature_run1.npy'))
der = np.squeeze(np.load('derivative.npy'))
X = np.concatenate([sol, der], axis=0).T
# y should be a 1D array with shape (n_samples,)
y = np.load('Doublingtime_label_test.npy')

# ------------------------------------------------------------
# train KNN
# Create a KNN classifier object
k_neighbors = 5  # You can choose the number of neighbors (k)
knn_classifier = KNeighborsClassifier(n_neighbors=k_neighbors)

# Perform k-fold cross-validation and get the accuracy scores
k_fold = 5  # Number of folds for cross-validation
cv_results = cross_validate(knn_classifier, X, y, cv=k_fold,
                        scoring=('accuracy', 'precision', 'recall', 'roc_auc'),
                        return_estimator=True)
                        # ,return_indices=True)

# Calculate the mean accuracy and standard deviation of accuracy scores
mean_accuracy = cv_results['test_accuracy'].mean()
mean_precision = cv_results['test_precision'].mean()
mean_recall = cv_results['test_recall'].mean()
mean_rocauc = cv_results['test_roc_auc'].mean()
print(f"Mean accuracy: {mean_accuracy:.2f}")
print(f"Mean precision: {mean_precision:.2f}")
print(f"Mean recall: {mean_recall:.2f}")
print(f"Mean ROC_AUC: {mean_rocauc:.2f}")

# # Save the trained model to a file
# model_filename = './trained_model/knn_model.joblib'
# joblib.dump(knn_classifier, model_filename)

# # ------------------------------------------------------------
# # train Fuzzy KNN
# k_neighbors = 5
# fknn_classifier = FuzzyKNN(k=k_neighbors)
#
# # # Perform k-fold cross-validation and get the accuracy scores
# # k_fold = 5  # Number of folds for cross-validation
# # cv_results = cross_validate(fknn_classifier, X, y, cv=k_fold,
# #                         scoring=('accuracy', 'precision', 'recall', 'roc_auc'),
# #                         return_estimator=True)
# #                         # ,return_indices=True)
# #
# # # Calculate the mean accuracy and standard deviation of accuracy scores
# # mean_accuracy = cv_results['test_accuracy'].mean()
# # mean_precision = cv_results['test_precision'].mean()
# # mean_recall = cv_results['test_recall'].mean()
# # mean_rocauc = cv_results['test_roc_auc'].mean()
# # print(f"Mean accuracy: {mean_accuracy:.2f}")
# # print(f"Mean precision: {mean_precision:.2f}")
# # print(f"Mean recall: {mean_recall:.2f}")
# # print(f"Mean ROC_AUC: {mean_rocauc:.2f}")
#
# # # Save the trained model to a file
# # model_filename = './trained_model/fknn_model.joblib'
# # joblib.dump(fknn_classifier, model_filename)

# ------------------------------------------------------------
# scatter truth and predict on test data
# use doubling time classificaiton as truth
truth_label = np.load('Doublingtime_label_test.npy')
sol = np.squeeze(np.load('sol_test.npy'))
der = np.squeeze(np.load('der_test.npy'))
test_data = np.concatenate([sol, der], axis=0).T

# KNN train and predict
knn_classifier.fit(X, y)
knn_label = knn_classifier.predict(test_data)
accuracy_knn, precision_knn, recall_knn = evaluate_predictions(truth_label, knn_label)
roc_auc_score_knn = roc_auc_score(truth_label, 1 - knn_classifier.predict_proba(test_data)[:, 1], average='macro')

print("Accuracy KNN:", accuracy_knn)
print("Precision KNN:", precision_knn)
print("Recall KNN:", recall_knn)
print("ROC_AUC KNN:", roc_auc_score_knn)

# # FKNN train and predict
# fknn_classifier.fit(X, y)
# fknn_result = fknn_classifier.predict(test_data)  # result: [(label, votes),...]
# fknn_label = np.array([fknn_result[i][0] for i in range(len(fknn_result))])
#
# # # reverse label 0 - 1
# # fknn_label = - fknn_label + 1
# accuracy_fknn, precision_fknn, recall_fknn = evaluate_predictions(truth_label, fknn_label)
# roc_auc_score_fknn = roc_auc_score(truth_label, 1 - fknn_classifier.predict_proba(test_data)[:, 1], average='macro')
#
# print("Accuracy FKNN:", accuracy_fknn)
# print("Precision FKNN:", precision_fknn)
# print("Recall FKNN:", recall_fknn)
# print("ROC_AUC FKNN:", roc_auc_score_fknn)

fig1, axes1 = plt.subplots(1, 3, figsize=(9, 3), subplot_kw=dict(projection='3d'))
colors = ['b', 'r', 'g', 'orange']
ncenters = 2
x_data = sol[0, :]
y_data = sol[1, :]
z_data = sol[2, :]

for j in range(ncenters):
    axes1[0].scatter(x_data[truth_label == j], y_data[truth_label == j], z_data[truth_label == j],
               label=f'Regime {j+1}', s=.5, color=colors[j])
# j=1
# axes1[0].scatter(x_data[truth_label == j], y_data[truth_label == j], z_data[truth_label == j],
#                  label=f'Regime {j + 1}', s=.5, color=colors[j])
axes1[0].set_xlabel('X')
axes1[0].set_ylabel('Y')
axes1[0].set_zlabel('Z')
axes1[0].set_title('Truth')

for j in range(ncenters):
    axes1[1].scatter(x_data[knn_label == j], y_data[knn_label == j], z_data[knn_label == j],
               label=f'Regime {j+1}', s=.5, color=colors[j])
# axes1[1].scatter(x_data[knn_label == j], y_data[knn_label == j], z_data[knn_label == j],
#            label=f'Regime {j+1}', s=.5, color=colors[j])
axes1[1].set_xlabel('X')
axes1[1].set_ylabel('Y')
axes1[1].set_zlabel('Z')
axes1[1].set_title('KNN')

plt.subplots_adjust(wspace=0.25)  # Increased the width of intervals among subplots
plt.show()
