"""
MIT License

Copyright (c) 2018 Christopher Tralie and Nathaniel Saul

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. 
"""

from itertools import cycle
import warnings

from scipy import sparse

import numpy as np
from sklearn.base import TransformerMixin
from sklearn.metrics.pairwise import pairwise_distances

from pyRipser import doRipsFiltrationDM as DRFDM
from pyRipser import doRipsFiltrationDMSparse as DRFDMSparse


def dpoint2pointcloud(X, i, metric):
    """
    Return the distance from the ith point in a Euclidean point cloud
    to the rest of the points
    Parameters
    ----------
    X: ndarray (n_samples, n_features)
        A numpy array of data 
    i: int
        The index of the point from which to return all distances
    metric: string or callable
        The metric to use when calculating distance between instances in a 
        feature array
    """
    ds = pairwise_distances(X, X[i, :][None, :], metric=metric).flatten()
    ds[i] = 0
    return ds


def get_greedy_perm(X, n_perm=None, distance_matrix=False, metric="euclidean"):
    """
    Compute a furthest point sampling permutation of a set of points
    Parameters
    ----------
    X: ndarray (n_samples, n_features)
        A numpy array of either data or distance matrix
    distance_matrix: bool
        Indicator that X is a distance matrix, if not we compute 
        distances in X using the chosen metric.
    n_perm: int
        Number of points to take in the permutation
    metric: string or callable
        The metric to use when calculating distance between instances in a 
        feature array
    Returns
    -------
    idx_perm: ndarray(n_perm)
        Indices of points in the greedy permutation
    lambdas: ndarray(n_perm)
        Covering radii at different points
    dperm2all: ndarray(n_perm, n_samples)
        Distances from points in the greedy permutation to points
        in the original point set
    """
    if not n_perm:
        n_perm = X.shape[0]
    # By default, takes the first point in the list to be the
    # first point in the permutation, but could be random
    idx_perm = np.zeros(n_perm, dtype=np.int64)
    lambdas = np.zeros(n_perm)
    if distance_matrix:
        dpoint2all = lambda i: X[i, :]
    else:
        dpoint2all = lambda i: dpoint2pointcloud(X, i, metric)
    ds = dpoint2all(0)
    dperm2all = [ds]
    for i in range(1, n_perm):
        idx = np.argmax(ds)
        idx_perm[i] = idx
        lambdas[i - 1] = ds[idx]
        dperm2all.append(dpoint2all(idx))
        ds = np.minimum(ds, dperm2all[-1])
    lambdas[-1] = np.max(ds)
    dperm2all = np.array(dperm2all)
    return (idx_perm, lambdas, dperm2all)


def ripser(
    X,
    maxdim=1,
    thresh=np.inf,
    coeff=2,
    distance_matrix=False,
    do_cocycles=False,
    metric="euclidean",
    n_perm=None,
):
    """Compute persistence diagrams for X.

    X can be a data set of points or a distance matrix. When using a data set
    as X it will be converted to a distance matrix using the metric specified.

    Parameters
    ----------

    X : ndarray (n_samples, n_features)
        A numpy array of either data or distance matrix (also pass `distance_matrix=True`). Can also be a sparse distance matrix of type scipy.sparse

    maxdim: int, optional, default 1
        Maximum homology dimension computed. Will compute all dimensions lower than and equal to this value.  For 1, H_0 and H_1 will be computed.

    thresh: float, default infinity
        Maximum distances considered when constructing filtration.  If infinity, compute the entire filtration.

    coeff: int prime, default 2
        Compute homology with coefficients in the prime field Z/pZ for p=coeff.

    distance_matrix: bool, optional, default False
        When True the input matrix X will be considered a distance matrix.

    do_cocycles: bool, optional, default False
        Computed cocycles will be available in the `cocycles` value
        of the return dictionary.

    metric: string or callable, optional, default "euclidean"
        Use this metric to compute distances between rows of X.

        "euclidean", "manhattan" and "cosine" are already provided metrics
        to choose from by using their name.

        You can provide a callable function and it will be used with two
        rows as arguments, it will be called once for each pair of rows in X.

        The computed distance will be available in the result dictionary under
        the key `dperm2all`.
    
    n_perm: int, optional, default None
        The number of points to subsample in a "greedy permutation,"
        or a furthest point sampling of the points.  These points
        will be used in lieu of the full point cloud for a faster
        computation, at the expense of some accuracy, which can 
        be bounded as a maximum bottleneck distance to all diagrams
        on the original point set

    Returns
    -------
    dict
        The result of the computation.

        .. note::
            Each list in `dgms` has a relative list in `cocycles`.

            >>> r = ripser(...)

            For each dimension ``d`` and index ``k`` then ``r['dgms'][d][k]``
            is the barcode associated to the representative cocycle
            ``r['cocycles'][d][k]``.

        The keys available in the dictionary are the:

            * ``dgms``: list (size maxdim) of ndarray (n_pairs, 2)
                For each dimension less than ``maxdim`` a list of persistence diagrams.
                Each persistent diagram is a pair (birth time, death time).
            * ``cocycles``: list (size maxdim) of list of ndarray
                For each dimension less than ``maxdim`` a list of representative cocycles.
                Each representative cocycle in dimension ``d`` is represented as a
                ndarray of ``(k,d+1)`` elements. Each non zero value of the cocycle
                is laid out in a row, first the ``d`` indices of the simplex and then
                the value of the cocycle on the simplex.  The indices of the simplex
                reference the original point cloud, even if a greedy permutation was used.
            * ``num_edges``: int
                The number of edges added during the computation
            * ``dperm2all``: ndarray(n_samples, n_samples) or ndarray (n_perm, n_samples) if n_perm
                The distance matrix used during the computation. When ``n_perm``
                is not None the distance matrix will only refers to the subsampled
                dataset.
            * ``idx_perm``: ndarray(n_perm) if ``n_perm`` > 0
                Index into the original point cloud of the points used
                as a subsample in the greedy permutation

                    >>> r = ripser(X, n_perm=k)
                    >>> subsampling = X[r['idx_perm']]

            * 'r_cover': float
                Covering radius of the subsampled points.
                If ``n_perm <= 0``, then the full point cloud was used and this is 0

    Examples
    --------

    .. code:: python

        from ripser import ripser, plot_dgms
        from sklearn import datasets
        from persim import plot_diagrams

        data = datasets.make_circles(n_samples=110)[0]
        dgms = ripser(data)['dgms']
        plot_diagrams(dgms, show = True)

    Raises
    ------

    ValueError
        If the distance matrix is not square.

    ValueError
        When using both a greedy permutation and a sparse distance matrix.

    ValueError
        When `n_perm` value is bigger than the number of rows in the matrix.

    ValueError
        When `n_perm` is non positive.

    Warns
    ----

        When using a square matrix without toggling `distance_matrix` to True.

        When there are more columns than rows (as each row is a different data point).

    """

    if distance_matrix:
        if not (X.shape[0] == X.shape[1]):
            raise ValueError("Distance matrix is not square")
    else:
        if X.shape[0] == X.shape[1]:
            warnings.warn(
                "The input matrix is square, but the distance_matrix "
                + "flag is off.  Did you mean to indicate that "
                + "this was a distance matrix?"
            )
        elif X.shape[0] < X.shape[1]:
            warnings.warn(
                "The input point cloud has more columns than rows; "
                + "did you mean to transpose?"
            )

    if n_perm and distance_matrix and sparse.issparse(X):
        raise ValueError(
            "Greedy permutation is not supported for sparse distance matrices"
        )
    if n_perm and n_perm > X.shape[0]:
        raise ValueError(
            "Number of points in greedy permutation is greater"
            + " than number of points in the point cloud"
        )
    if n_perm and n_perm < 0:
        raise ValueError(
            "Should be a strictly positive number of points in the greedy permutation"
        )

    idx_perm = np.arange(X.shape[0])
    r_cover = 0.0
    doing_permutation = False
    if n_perm and n_perm < X.shape[0]:
        doing_permutation = True
        idx_perm, lambdas, dperm2all = get_greedy_perm(
            X, n_perm=n_perm, distance_matrix=distance_matrix, metric=metric
        )
        r_cover = lambdas[-1]
        dm = dperm2all[:, idx_perm]
    else:
        if distance_matrix:
            dm = X
        else:
            dm = pairwise_distances(X, metric=metric)
        dperm2all = dm

    n_points = dm.shape[0]
    if not sparse.issparse(dm) and np.sum(np.abs(dm.diagonal()) > 0) > 0:
        # If any of the diagonal elements are nonzero,
        # convert to sparse format, because currently
        # that's the only format that handles nonzero
        # births
        dm = sparse.coo_matrix(dm)

    if sparse.issparse(dm):
        if sparse.isspmatrix_coo(dm):
            # If the matrix is already COO, we need to order the row and column indices
            # lexicographically to avoid errors. See issue #103
            row, col, data = dm.row, dm.col, dm.data
            lex_sort_idx = np.lexsort((col, row))
            row, col, data = row[lex_sort_idx], col[lex_sort_idx], data[lex_sort_idx]
        else:
            # Lexicographic ordering is performed by scipy upon conversion to COO
            coo = dm.tocoo()
            row, col, data = coo.row, coo.col, coo.data
        res = DRFDMSparse(
            row.astype(dtype=np.int32, order="C"),
            col.astype(dtype=np.int32, order="C"),
            np.array(data, dtype=np.float32, order="C"),
            n_points,
            maxdim,
            thresh,
            coeff,
            do_cocycles
        )
    else:
        I, J = np.meshgrid(np.arange(n_points), np.arange(n_points))
        DParam = np.array(dm[I > J], dtype=np.float32)
        res = DRFDM(DParam, maxdim, thresh, coeff, do_cocycles)

    # Unwrap persistence diagrams
    dgms = res["births_and_deaths_by_dim"]
    for dim in range(len(dgms)):
        N = int(len(dgms[dim]) / 2)
        dgms[dim] = np.reshape(np.array(dgms[dim]), [N, 2])

    # Unwrap cocycles
    cocycles = []
    for dim in range(len(res["cocycles_by_dim"])):
        cocycles.append([])
        for j in range(len(res["cocycles_by_dim"][dim])):
            ccl = res["cocycles_by_dim"][dim][j]
            n = int(len(ccl) / (dim + 2))
            ccl = np.reshape(np.array(ccl, dtype=np.int64), [n, dim + 2])
            ccl[:, -1] = np.mod(ccl[:, -1], coeff)
            if doing_permutation:
                # Retain original indices in the point cloud
                ccl[:, 0:-1] = idx_perm[ccl[:, 0:-1]]
            cocycles[dim].append(ccl)

    ret = {
        "dgms": dgms,
        "cocycles": cocycles,
        "num_edges": res["num_edges"],
        "dperm2all": dperm2all,
        "idx_perm": idx_perm,
        "r_cover": r_cover,
    }
    return ret


def lower_star_img(img):
    """
    Construct a lower star filtration on an image

    Parameters
    ----------
    img: ndarray (M, N)
        An array of single channel image data

    Returns
    -------
    I: ndarray (K, 2)
        A 0-dimensional persistence diagram corresponding to the sublevelset filtration
    """
    m, n = img.shape

    idxs = np.arange(m * n).reshape((m, n))

    I = idxs.flatten()
    J = idxs.flatten()
    V = img.flatten()

    # Connect 8 spatial neighbors
    tidxs = np.ones((m + 2, n + 2), dtype=np.int64) * np.nan
    tidxs[1:-1, 1:-1] = idxs

    tD = np.ones_like(tidxs) * np.nan
    tD[1:-1, 1:-1] = img

    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:

            if di == 0 and dj == 0:
                continue

            thisJ = np.roll(np.roll(tidxs, di, axis=0), dj, axis=1)
            thisD = np.roll(np.roll(tD, di, axis=0), dj, axis=1)
            thisD = np.maximum(thisD, tD)

            # Deal with boundaries
            boundary = ~np.isnan(thisD)
            thisI = tidxs[boundary]
            thisJ = thisJ[boundary]
            thisD = thisD[boundary]

            I = np.concatenate((I, thisI.flatten()))
            J = np.concatenate((J, thisJ.flatten()))
            V = np.concatenate((V, thisD.flatten()))

    sparseDM = sparse.coo_matrix((V, (I, J)), shape=(idxs.size, idxs.size))

    return ripser(sparseDM, distance_matrix=True, maxdim=0)["dgms"][0]


class Rips(TransformerMixin):
    """ sklearn style class interface for :code:`ripser` with :code:`fit` and :code:`transform` methods..

    Parameters
    ----------
    maxdim: int, optional, default 1
        Maximum homology dimension computed. Will compute all dimensions 
        lower than and equal to this value. 
        For 1, H_0 and H_1 will be computed.

    thresh: float, default infinity
        Maximum distances considered when constructing filtration. 
        If infinity, compute the entire filtration.

    coeff: int prime, default 2
        Compute homology with coefficients in the prime field Z/pZ for p=coeff.

    do_cocycles: bool
        Indicator of whether to compute cocycles, if so, we compute and store
        cocycles in the `cocycles_` dictionary Rips member variable

    n_perm: int
        The number of points to subsample in a "greedy permutation,"
        or a furthest point sampling of the points.  These points
        will be used in lieu of the full point cloud for a faster
        computation, at the expense of some accuracy, which can 
        be bounded as a maximum bottleneck distance to all diagrams
        on the original point set
    
    verbose: boolean
        Whether to print out information about this object
        as it is constructed

    Attributes
    ----------
    `dgm_`: list of ndarray, each shape (n_pairs, 2)
        After `transform`, `dgm_` contains computed persistence diagrams in
        each dimension
    
    cocycles_: list (size maxdim) of list of ndarray
        A list of representative cocycles in each dimension.  The list 
        in each dimension is parallel to the diagram in that dimension;
        that is, each entry of the list is a representative cocycle of
        the corresponding point expressed as an ndarray(K, d+1), where K is
        the number of nonzero values of the cocycle and d is the dimension
        of the cocycle.  The first d columns of each array index into
        the simplices of the (subsampled) point cloud, and the last column
        is the value of the cocycle at that simplex
    
    dperm2all_: ndarray(n_samples, n_samples) or ndarray (n_perm, n_samples) if n_perm
        The distance matrix used in the computation if n_perm is none.
        Otherwise, the distance from all points in the permutation to
        all points in the dataset
    
    metric_: string or callable
        The metric to use when calculating distance between instances in a 
        feature array. If metric is a string, it must be one of the options 
        specified in pairwise_distances, including "euclidean", "manhattan", 
        or "cosine". Alternatively, if metric is a callable function, it is 
        called on each pair of instances (rows) and the resulting value 
        recorded. The callable should take two arrays from X as input and 
        return a value indicating the distance between them.
    
    num_edges: int
        The number of edges added during the computation
    
    idx_perm: ndarray(n_perm) if n_perm > 0
        Index into the original point cloud of the points used
        as a subsample in the greedy permutation

    r_cover: float
        Covering radius of the subsampled points.  
        If n_perm <= 0, then the full point cloud was used and this is 0
    Examples
    --------
     .. code:: python

        from ripser import Rips
        import tadasets

        data = tadasets.dsphere(n=110, d=2)[0]
        rips = Rips()
        rips.transform(data)
        rips.plot()
    """

    def __init__(
        self,
        maxdim=1,
        thresh=np.inf,
        coeff=2,
        do_cocycles=False,
        n_perm=None,
        verbose=True,
    ):
        self.maxdim = maxdim
        self.thresh = thresh
        self.coeff = coeff
        self.do_cocycles = do_cocycles
        self.n_perm = n_perm
        self.verbose = verbose

        # Internal variables
        self.dgms_ = None
        self.cocycles_ = None
        self.dperm2all_ = None  # Distance matrix
        self.metric_ = None
        self.num_edges_ = None  # Number of edges added
        self.idx_perm_ = None
        self.r_cover_ = 0.0

        if self.verbose:
            print(
                "Rips(maxdim={}, thresh={}, coeff={}, do_cocycles={}, n_perm = {}, verbose={})".format(
                    maxdim, thresh, coeff, do_cocycles, n_perm, verbose
                )
            )

    def transform(self, X, distance_matrix=False, metric="euclidean"):
        result = ripser(
            X,
            maxdim=self.maxdim,
            thresh=self.thresh,
            coeff=self.coeff,
            do_cocycles=self.do_cocycles,
            distance_matrix=distance_matrix,
            metric=metric,
            n_perm=self.n_perm,
        )
        self.dgms_ = result["dgms"]
        self.num_edges_ = result["num_edges"]
        self.dperm2all_ = result["dperm2all"]
        self.idx_perm_ = result["idx_perm"]
        self.cocycles_ = result["cocycles"]
        self.r_cover_ = result["r_cover"]
        return self.dgms_

    def fit_transform(self, X, distance_matrix=False, metric="euclidean"):
        """
        Compute persistence diagrams for X data array and return the diagrams.

        Parameters
        ----------
        X: ndarray (n_samples, n_features)
            A numpy array of either data or distance matrix.

        distance_matrix: bool
            Indicator that X is a distance matrix, if not we compute a 
            distance matrix from X using the chosen metric.

        metric: string or callable
            The metric to use when calculating distance between instances in a 
            feature array. If metric is a string, it must be one of the options 
            specified in pairwise_distances, including "euclidean", "manhattan", 
            or "cosine". Alternatively, if metric is a callable function, it is 
            called on each pair of instances (rows) and the resulting value 
            recorded. The callable should take two arrays from X as input and 
            return a value indicating the distance between them.

        Returns
        -------
        dgms: list (size maxdim) of ndarray (n_pairs, 2)
            A list of persistence diagrams, one for each dimension less 
            than maxdim. Each diagram is an ndarray of size (n_pairs, 2) with 
            the first column representing the birth time and the second column 
            representing the death time of each pair.
        """
        self.transform(X, distance_matrix, metric)
        return self.dgms_

    def plot(
        self,
        diagrams=None,
        *args,
        **kwargs
    ):
        """A helper function to plot persistence diagrams. 

        Parameters
        ----------
        diagrams: ndarray (n_pairs, 2) or list of diagrams
            A diagram or list of diagrams as returned from self.fit. 
            If diagram is None, we use `self.dgm_` for plotting. 
            If diagram is a list of diagrams, then plot all on the same plot 
            using different colors.

        plot_only: list of numeric
            If specified, an array of only the diagrams that should be plotted.

        title: string, default is None
            If title is defined, add it as title of the plot.

        xy_range: list of numeric [xmin, xmax, ymin, ymax]
            User provided range of axes. This is useful for comparing 
            multiple persistence diagrams.

        labels: string or list of strings
            Legend labels for each diagram. 
            If none are specified, we use H_0, H_1, H_2,... by default.

        colormap: string, default is 'default'
            Any of matplotlib color palettes. 
            Some options are 'default', 'seaborn', 'sequential'. 
            See all available styles with
            
            .. code:: python

                import matplotlib as mpl
                print(mpl.styles.available)

        size: numeric, default is 20
            Pixel size of each point plotted.

        ax_color: any valid matplitlib color type. 
            See [https://matplotlib.org/api/colors_api.html](https://matplotlib.org/api/colors_api.html) for complete API.

        diagonal: bool, default is True
            Plot the diagonal x=y line.

        lifetime: bool, default is False. If True, diagonal is turned to False.
            Plot life time of each point instead of birth and death. 
            Essentially, visualize (x, y-x).

        legend: bool, default is True
            If true, show the legend.

        show: bool, default is True
            Call plt.show() after plotting. 
            If you are using self.plot() as part of a subplot, 
            set show=False and call plt.show() only once at the end.
        """
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import persim

        if diagrams is None:
            # Allow using transformed diagrams as default
            diagrams = self.dgms_

        persim.plot_diagrams(
            diagrams,
            *args,
            **kwargs
        )


__all__ = ["Rips", "ripser", "lower_star_img"]
