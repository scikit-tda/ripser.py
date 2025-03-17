import numpy as np

from ripser import lower_star_img

class TestLowerStar:
    def test_returns_dgm(self):
        img = np.random.random((100,400))
        dgm = lower_star_img(img)
        assert isinstance(dgm, np.ndarray)
        assert dgm.shape[1] == 2
        assert np.all(dgm[:,1] >= dgm[:,0])
    
    def test_single_point(self):
        img = np.ones((11,11))
        img[6,6] = 0

        dgm = lower_star_img(img)

        assert dgm.shape[0] == 1
        assert tuple(dgm[0]) == (0,np.inf)
    
    def test_zero_edge_bug(self):
        img = np.ones((5, 5))
        img[1:-1, 1:-1] = 0
        img[2, 2] = 1
        r1 = -lower_star_img(-img)
        r2 = 1-lower_star_img(1-img)
        assert(np.allclose(r1, r2))