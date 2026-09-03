import numpy as np

from dbm2d.exact_riemann import sod_exact


def test_sod_exact_initial_condition():
    x = np.array([0.25, 0.75])
    sol = sod_exact(x, 0.0, 0.5, (1.0, 0.0, 1.0), (0.125, 0.0, 0.1))
    assert np.allclose(sol["rho"], [1.0, 0.125])
    assert np.allclose(sol["p"], [1.0, 0.1])
