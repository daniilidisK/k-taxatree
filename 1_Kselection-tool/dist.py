from sklearn.metrics import pairwise_distances

def dist_python(X):

    return pairwise_distances(X,  metric='hamming', n_jobs = 4)
