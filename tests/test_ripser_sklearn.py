import pytest
import numpy as np

from ripser import Rips
from sklearn import datasets
from sklearn.metrics.pairwise import pairwise_distances
from scipy import sparse


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
        from ripser import Rips

        assert 1

    def test_instantiate(self):
        rip = Rips()
        assert rip is not None


class TestTransform:
    def test_input_warnings(self):

        rips = Rips()
        data = np.random.random((3, 10))

        with pytest.warns(UserWarning, match="has more columns than rows") as w:
            rips.transform(data)

        data = np.random.random((3, 3))
        with pytest.warns(
            UserWarning, match="input matrix is square, but the distance_matrix"
        ) as w:
            rips.transform(data)

    def test_non_square_dist_matrix(self):

        rips = Rips()
        data = np.random.random((3, 10))

        with pytest.raises(Exception):
            rips.transform(data, distance_matrix=True)


class TestParams:
    def test_verbose_true(self):
        data = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]])
        rips = Rips(verbose=True)
        dgm = rips.fit_transform(data)
        assert len(dgm) == 2
        assert dgm[0].shape == (4, 2)
        assert dgm[1].shape == (1, 2)

    def test_verbose_false(self):
        data = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]])
        rips = Rips(verbose=False)
        dgm = rips.fit_transform(data)
        assert dgm[0].shape == (4, 2)
        assert dgm[1].shape == (1, 2)

    def test_defaults(self):
        data = np.random.random((100, 3))

        rips = Rips()
        dgm = rips.fit_transform(data)

        assert len(dgm) == 2
        assert rips.coeff == 2

    def test_coeff(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        rips3 = Rips(coeff=3)
        dgm3 = rips3.fit_transform(data)

        rips2 = Rips(coeff=2)
        dgm2 = rips2.fit_transform(data)
        assert (
            dgm2 is not dgm3
        ), "This is a vacuous assertion, we only care that the above operations did not throw errors"

    def test_maxdim(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        # maxdim refers to the max H_p class, generate all less than

        rips0 = Rips(maxdim=0)
        dgm0 = rips0.fit_transform(data)
        assert len(dgm0) == 1

        rips1 = Rips(maxdim=1)
        dgm1 = rips1.fit_transform(data)
        assert len(dgm1) == 2

        rips2 = Rips(maxdim=2)
        dgm2 = rips2.fit_transform(data)
        assert len(dgm2) == 3

    def test_thresh(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        rips0 = Rips(thresh=0.1)
        dgm0 = rips0.fit_transform(data)

        rips1 = Rips(thresh=1)
        dgm1 = rips1.fit_transform(data)

        # Barcode of H_1 diagram will be smaller, right?
        assert len(dgm0[1]) < len(dgm1[1]), "Usually"

    def test_sparse(self):
        np.random.seed(10)
        thresh = 1.1

        # Do dense filtration with threshold
        data = (
            datasets.make_circles(n_samples=100)[0]
            + 5 * datasets.make_circles(n_samples=100)[0]
        )
        rips0 = Rips(thresh=thresh, maxdim=1)
        dgms0 = rips0.fit_transform(data)

        # Convert to sparse matrix first based on threshold,
        # then do full filtration
        rips1 = Rips(maxdim=1)
        D = makeSparseDM(data, thresh)
        dgms1 = rips1.fit_transform(D, distance_matrix=True)

        # The same number of edges should have been added
        assert rips0.num_edges_ == rips1.num_edges_

        I10 = dgms0[1]
        I11 = dgms1[1]
        idx = np.argsort(I10[:, 0])
        I10 = I10[idx, :]
        idx = np.argsort(I11[:, 0])
        I11 = I11[idx, :]
        assert np.allclose(I10, I11)

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
        rips1 = Rips(maxdim=1)
        rips2 = Rips(maxdim=1, n_perm=10)
        h11 = rips1.fit_transform(X)[1]
        h12 = rips2.fit_transform(X)[1]
        assert rips2.r_cover_ > 0
        assert np.max(np.abs(h11 - h12)) <= 2 * rips2.r_cover_
