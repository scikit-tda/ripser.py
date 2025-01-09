import numpy as np
from ripser import estimate_ripser_memory

def test_memory_estimation():
    # Create a small point cloud
    n_points = 100
    n_features = 3
    X = np.random.rand(n_points, n_features)
    
    # Test memory estimation for different dimensions
    mem_dim0 = estimate_ripser_memory(X, maxdim=0)
    mem_dim1 = estimate_ripser_memory(X, maxdim=1)
    mem_dim2 = estimate_ripser_memory(X, maxdim=2)
    
    # Memory should increase with dimension
    assert mem_dim0 < mem_dim1 < mem_dim2
    
    # Test with distance matrix
    D = np.random.rand(n_points, n_points)
    D = (D + D.T) / 2  # Make it symmetric
    np.fill_diagonal(D, 0)  # Zero diagonal
    
    mem_dist = estimate_ripser_memory(D, maxdim=1, distance_matrix=True)
    assert mem_dist > 0

def test_memory_scaling():
    # Test how memory scales with number of points
    n_points_small = 50
    n_points_large = 100
    n_features = 3
    
    X_small = np.random.rand(n_points_small, n_features)
    X_large = np.random.rand(n_points_large, n_features)
    
    mem_small = estimate_ripser_memory(X_small, maxdim=1)
    mem_large = estimate_ripser_memory(X_large, maxdim=1)
    
    # Memory should scale super-linearly with number of points
    assert mem_large > 4 * mem_small  # Since complexity is O(n^3)