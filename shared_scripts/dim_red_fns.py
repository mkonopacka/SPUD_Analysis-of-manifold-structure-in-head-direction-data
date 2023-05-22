'''November 26th 2018
Functions for dimensionality reduction.'''
import numpy as np
import numpy.linalg as la
from sklearn import decomposition, manifold

def run_dim_red(inp_data: np.ndarray, params: dict, method='iso', stabilize: bool = True) -> np.ndarray:
    '''Run dimensionality reduction on input data.

    Args:
        inp_data: Input data, shape (n_samples, n_features)
        params: Dictionary of parameters for dimensionality reduction. For isomap, the keys are `n_neighbors` and `target_dim`.
        method: Dimensionality reduction method. Currently only `iso` (isomap) is implemented.
        stabilize: Whether to apply variance stabilization to the input data.

    Returns:
        proj_data: Projected data, shape (n_samples, target_dim)
    '''
    # Variance stabilization option included, since we're usually working with Poisson-like data
    if stabilize:
        data_to_use = np.sqrt(inp_data)
    else:
        data_to_use = inp_data.copy()
    # Isomap
    if method == 'iso':
        iso_instance = manifold.Isomap(n_neighbors = params['n_neighbors'], n_components = params['target_dim'])
        proj_data = iso_instance.fit_transform(data_to_use)
    else:
        raise NotImplementedError(f"Method {method} not implemented.") # added by me 22.05.2023
    return proj_data