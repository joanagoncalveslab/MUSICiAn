import numpy as np
import pandas as pd
from numpy.linalg import svd
from skbio.stats.composition import ilr
from sklearn.preprocessing import StandardScaler
from sklearn.covariance import MinCovDet, EmpiricalCovariance
from scipy.stats import chi2

def impute_zero_values(profiles, detection_limit):
    profiles = profiles.replace({0: np.nan})
    # smallest_value = profiles.min().min()
    M = len(profiles.index)
    N = len(profiles.columns)
    ran = pd.DataFrame(np.random.uniform(low=0.1 * detection_limit, high=detection_limit, size=(M,N)), columns=profiles.columns, index=profiles.index)
    ran.update(profiles)
    ran.sum(axis=1)
    profiles_imp = ran
    print("{:.2f} of our data is zero".format((profiles_imp.isnull()).sum().sum()/(profiles_imp.shape[0] * profiles_imp.shape[1])))
    return profiles_imp

def outlier_detection(X, cov_estimator, transform_data=None):
    ss = StandardScaler().fit(X)
    X = ss.transform(X)
    cov = cov_estimator.fit(X)
    distances = cov.mahalanobis(X)
    pvalues = 1 - chi2.cdf(distances, X.shape[1] - 1)
    T = cov.location_
    C = cov.covariance_
    G_z, L_z, G_z_inv = svd(C) 
    pca_scores = X @ G_z
    loadings = G_z
    explained_variance_ratio = L_z/L_z.sum()

    if transform_data is not None:
        TD = transform_data
        TD = ss.transform(TD)
        transform_data = TD @ loadings
        transform_data_distances = cov.mahalanobis(TD)
        transform_data_pvalues = 1 - chi2.cdf(transform_data_distances, X.shape[1] - 1)
        return pca_scores, loadings, explained_variance_ratio, distances, pvalues, transform_data, transform_data_distances, transform_data_pvalues, cov

    return pca_scores, loadings, explained_variance_ratio, distances, pvalues, ss, cov


# Outlier detection
def logratio_outlier_detection(X, cov_estimator, transform_data=None):
    Z = ilr(X.values)

    cov = cov_estimator.fit(Z) 

    T = cov.location_
    C = cov.covariance_
    distances = cov.mahalanobis(Z)
    pvalues = 1 - chi2.cdf(distances, Z.shape[1] - 1)

    # For mathematical details on below, see Peter Filzmoser; Karel Hron; Clemens Reimann (2009). Principal component analysis for compositional data with outliers. 

    G_z, L_z, G_z_inv = svd(C) # (11)
    D = X.shape[1] 
    V = np.array([np.sqrt(i/(i+1)) * np.concatenate((np.repeat(1/i, i), np.array([-1]), np.repeat(0, D - i - 1))) for i in range(1, D)]).T # (4)

    C_y = V @ G_z * L_z @ G_z_inv @ V.T 
    G_y = V @ G_z # (13)

    Z_t = (Z - T.T) @ G_z # (10)
    Y_t = Z_t @ V.T # (12)

    pca_scores = Y_t
    loadings = G_y
    explained_variance_ratio = L_z/L_z.sum()

    if transform_data is not None:
        TD = transform_data
        Z2 = ilr(TD.values) 
        Z2_t = (Z2 - T.T) @ G_z 
        Y2_t = Z2_t @ V.T
        transform_data = Y2_t
        transform_data_distances = cov.mahalanobis(Z2)
        transform_data_pvalues = 1 - chi2.cdf(transform_data_distances, Z.shape[1] - 1)

        return pca_scores, loadings, explained_variance_ratio, distances, pvalues, transform_data, transform_data_distances, transform_data_pvalues, cov

    return pca_scores, loadings, explained_variance_ratio, distances, pvalues, cov

