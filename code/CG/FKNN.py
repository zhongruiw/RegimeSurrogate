import operator
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score
from sklearn.base import BaseEstimator, ClassifierMixin
from numba import jit


# Define the distance calculation function using Numba's jit decorator
@jit(nopython=True)
def calculate_distances(X, x, n):
    distances = np.zeros(n)
    for i in range(n):
        distances[i] = np.linalg.norm(X[i] - x)
    return distances

class FuzzyKNN(BaseEstimator, ClassifierMixin):
    def __init__(self, k=3, plot=False):
        self.k = k
        self.plot = plot

    def fit(self, X, y=None):
        self._check_params(X, y)
        self.X = X
        self.y = y

        self.xdim = len(self.X[0])
        self.n = len(y)

        classes = list(set(y))
        classes.sort()
        self.classes = classes

        self.df = pd.DataFrame(self.X)
        self.df['y'] = self.y

        self.memberships = self._compute_memberships()

        self.df['membership'] = self.memberships

        self.fitted_ = True
        return self

    def predict(self, X):
        if self.fitted_ == None:
            raise Exception('predict() called before fit()')
        else:
            m = 2
            y_pred = []
            epsilon = 1e-10 # Small epsilon value to prevent division by zero

            for x in X:
                neighbors = self._find_k_nearest_neighbors(pd.DataFrame.copy(self.df), x)

                votes = {}
                for c in self.classes:
                    den = 0
                    for n in range(self.k):
                        dist = np.linalg.norm(x - neighbors.iloc[n, 0:self.xdim])
                        den += 1 / (dist ** (2 / (m - 1)) + epsilon)

                    neighbors_votes = []
                    for n in range(self.k):
                        dist = np.linalg.norm(x - neighbors.iloc[n, 0:self.xdim])
                        num = (neighbors.iloc[n].membership[c]) / (dist ** (2 / (m - 1))+ epsilon)

                        vote = num / den
                        neighbors_votes.append(vote)
                    votes[c] = np.sum(neighbors_votes)

                pred = max(votes.items(), key=operator.itemgetter(1))[0]
                y_pred.append((pred, votes))

            return y_pred

    def score(self, X, y):
        if self.fitted_ == None:
            raise Exception('score() called before fit()')
        else:
            predictions = self.predict(X)
            y_pred = [t[0] for t in predictions]
            confidences = [t[1] for t in predictions]

            return accuracy_score(y_pred=y_pred, y_true=y)

    def _find_k_nearest_neighbors(self, df, x):
        X = df.iloc[:, 0:self.xdim].values

        df['distances'] = calculate_distances(X, x, self.n)

        df.sort_values(by='distances', ascending=True, inplace=True)
        neighbors = df.iloc[0:self.k]

        return neighbors

    def _get_counts(self, neighbors):
        groups = neighbors.groupby('y')
        counts = {group[1]['y'].iloc[0]: group[1].count()[0] for group in groups}

        return counts

    def _compute_memberships(self):
        memberships = []
        for i in range(self.n):
            x = self.X[i]
            y = self.y[i]

            neighbors = self._find_k_nearest_neighbors(pd.DataFrame.copy(self.df), x)
            counts = self._get_counts(neighbors)

            membership = dict()
            for c in self.classes:
                try:
                    uci = 0.49 * (counts[c] / self.k)
                    if c == y:
                        uci += 0.51
                    membership[c] = uci
                except:
                    membership[c] = 0

            memberships.append(membership)
        return memberships

    def _check_params(self, X, y):
        if type(self.k) != int:
            raise Exception('"k" should have type int')
        elif self.k >= len(y):
            raise Exception('"k" should be less than no of feature sets')
        elif self.k % 2 == 0:
            raise Exception('"k" should be odd')

        if type(self.plot) != bool:
            raise Exception('"plot" should have type bool')