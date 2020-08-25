import pytest
import numpy as np

from ripser import ripser
from sklearn import datasets
from sklearn.metrics.pairwise import pairwise_distances
from scipy import sparse
import itertools


def makeSparseDM(X, thresh):
    """
    Helper function to make a sparse distance matrix
    """
    N = X.shape[0]
    D = pairwise_distances(X, metric="euclidean")
    [I, J] = np.meshgrid(np.arange(N), np.arange(N))
    I = I[D <= thresh]
    J = J[D <= thresh]
    V = D[D <= thresh]
    return sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr()


class TestLibrary:
    # Does the library install in scope? Are the objects in scope?
    def test_import(self):
        import ripser
        from ripser import ripser, Rips
        assert 1


class TestTransform:
    def test_input_warnings(self):
        data = np.random.random((3, 10))

        with pytest.warns(UserWarning, match="has more columns than rows") as w:
            ripser(data)

        data = np.random.random((3, 3))
        with pytest.warns(
            UserWarning, match="input matrix is square, but the distance_matrix"
        ) as w:
            ripser(data)

    def test_non_square_dist_matrix(self):
        data = np.random.random((3, 10))

        with pytest.raises(Exception, match="Distance matrix is not square"):
            ripser(data, distance_matrix=True)

    def test_too_many_perm(self):
        X = np.zeros((10, 2))

        with pytest.raises(Exception, match="points in greedy permutation is greater"):
            ripser(X, n_perm=11)

    def test_too_few_perm(self):
        X = np.zeros((10, 2))

        with pytest.raises(Exception, match="Should be a strictly positive number"):
            ripser(X, n_perm=-1)

    def test_sparse_greedyperm(self):
        D = np.zeros((4, 4))

        with pytest.raises(Exception):
            ripser(sparse.coo_matrix(D), distance_matrix=True, n_perm=2)


class TestParams:
    def test_defaults(self):
        data = np.random.random((100, 3))
        dgms = ripser(data)["dgms"]
        assert len(dgms) == 2

    def test_coeff(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        dgm3 = ripser(data, coeff=3)["dgms"]
        dgm2 = ripser(data)["dgms"]
        assert (
            dgm2 is not dgm3
        ), "This is a vacuous assertion, we only care that the above operations did not throw errors"

    def test_maxdim(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        # maxdim refers to the max H_p class, generate all less than
        dgms0 = ripser(data, maxdim=0)["dgms"]
        assert len(dgms0) == 1

        dgms1 = ripser(data)["dgms"]
        assert len(dgms1) == 2

        dgms2 = ripser(data, maxdim=2)["dgms"]
        assert len(dgms2) == 3

    def test_thresh(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        dgms0 = ripser(data, thresh=0.1)["dgms"]
        dgms1 = ripser(data)["dgms"]

        # Barcode of H_1 diagram will be smaller, right?
        assert len(dgms0[1]) < len(dgms1[1]), "Usually"

    def test_sparse(self):
        np.random.seed(10)
        thresh = 1.1

        # Do dense filtration with threshold
        data = (
            datasets.make_circles(n_samples=100)[0]
            + 5 * datasets.make_circles(n_samples=100)[0]
        )
        res0 = ripser(data, thresh=thresh)

        # Convert to sparse matrix first based on threshold,
        # then do full filtration
        D = makeSparseDM(data, thresh)
        res1 = ripser(D, distance_matrix=True)

        # The same number of edges should have been added
        assert res0["num_edges"] == res1["num_edges"]

        dgms0 = res0["dgms"]
        dgms1 = res1["dgms"]
        I10 = dgms0[1]
        I11 = dgms1[1]
        idx = np.argsort(I10[:, 0])
        I10 = I10[idx, :]
        idx = np.argsort(I11[:, 0])
        I11 = I11[idx, :]
        assert np.allclose(I10, I11)

    def test_sphere_sparse_H2(self):
        n = 3
        segment = [np.linspace(0, 1, 5)]
        endpoints = [np.linspace(0, 1, 2)]
        face = segment * (n - 1) + endpoints
        vertices = []
        for k in range(n):
            vertices.extend(itertools.product(*(face[k:] + face[:k])))
        coords = np.array(vertices)
        thresh = 1.5
        D = makeSparseDM(coords, thresh)
        rips = ripser(D, distance_matrix=True, maxdim=2, thresh=thresh)
        I2 = rips["dgms"][2]
        assert I2.shape[0] == 1
        assert np.allclose(1.0, I2[0, 1])

    def test_full_nonzerobirths(self):
        D = np.array([[1.0, 3.0], [3.0, 2.0]])
        h0 = ripser(D, distance_matrix=True, maxdim=0)["dgms"][0]
        h0 = h0[np.argsort(h0[:, 0]), :]
        assert h0[0, 0] == 1
        assert np.isinf(h0[0, 1])
        assert h0[1, 0] == 2
        assert h0[1, 1] == 3

    def test_greedyperm_dm_vs_pc(self):
        """
        Test that point cloud and distance matrix on point cloud
        give the same persistence diagrams and bottleneck bound
        """
        np.random.seed(100)
        X = np.random.randn(100, 3)
        D = pairwise_distances(X, metric="euclidean")
        dgms1 = ripser(X, n_perm=20)["dgms"]
        dgms2 = ripser(D, distance_matrix=True, n_perm=20)["dgms"]
        for I1, I2 in zip(dgms1, dgms2):
            I1 = I1[np.argsort(I1[:, 0] - I1[:, 1]), :]
            I2 = I2[np.argsort(I2[:, 0] - I2[:, 1]), :]
            assert np.allclose(I1, I2)

    def test_greedyperm_circlebottleneck(self):
        """
        Test a relationship between the bottleneck
        distance and the covering radius for a simple case
        where computing the bottleneck distance is trivial
        """
        N = 40
        np.random.seed(N)
        t = 2 * np.pi * np.random.rand(N)
        X = np.array([np.cos(t), np.sin(t)]).T
        res1 = ripser(X)
        res2 = ripser(X, n_perm=10)
        idx = res2["idx_perm"]
        h11 = res1["dgms"][1]
        h12 = res2["dgms"][1]
        assert res2["r_cover"] > 0
        assert np.max(np.abs(h11 - h12)) <= 2 * res2["r_cover"]
    
    def test_cocycle_indices_greedyperm(self):
        """
        Make sure the original vertex indices are retained in the
        cocycle simplex specifications, even when a greedy
        permutation is in use
        """
        N = 1000
        n_perm = 50
        np.random.seed(N)
        t = 2 * np.pi * np.random.rand(N)
        X = np.array([np.cos(t), np.sin(t)]).T
        res1 = ripser(X, n_perm=n_perm, do_cocycles=True, maxdim=1)
        cocycles1 = res1['cocycles'][1]
        idx_perm = res1['idx_perm']
        X = X[idx_perm, :]
        res2 = ripser(X, do_cocycles=True, maxdim=1)
        cocycles2 = res2['cocycles'][1]
        for cc1, cc2 in zip(cocycles1, cocycles2):
            assert cc1.shape[0] == cc2.shape[0]
            cc2[:, 0:-1] = idx_perm[cc2[:, 0:-1]]
            assert np.allclose(cc1, cc2)
    
    def test_sparse_format(self):
        """
        Test to make sure different formats for sparse matrices
        yield the same answer.  Test courtesy of Umberto Lupo (@ulupo)
        """
        data = np.array([6., 8., 2., 4., 5., 9., 10., 3., 1., 1.])
        row = np.array([0, 0, 0, 0, 1, 1, 1, 2, 2, 3])
        col = np.array([4, 1, 3, 2, 4, 3, 2, 3, 4, 4])
        dm = sparse.coo_matrix((data, (row, col)), shape=(5, 5))
        dgms1 = ripser(dm, distance_matrix=True)['dgms']
        dgms2 = ripser(dm.tocsr(), distance_matrix=True)['dgms']
        for dgm1k, dgm2k in zip(dgms1, dgms2):
            assert(np.allclose(dgm1k, dgm2k))
