from sklearn.metrics import pairwise_distances
#rom scipy.spatial import distance_matrix

def dist_python(X):
    
    return pairwise_distances(X,  metric='minkowski', p=1, n_jobs = 4)
