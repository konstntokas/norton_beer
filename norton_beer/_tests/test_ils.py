import numpy as np
import norton_beer.ils as ilsfun


def test_ils_analytical():
    """
    this function tests the the analytical ils generation against the numerical ils
    generation for Norton-Beer
    """
    ifglen = 1.5
    pars = np.round(np.arange(1.1, 2.01, 0.1), 1)
    ils_xs = [
        np.arange(-10, 10, 0.1),
        np.arange(-10, 10.01, 0.1)
    ]
    tags = ["even", "odd"]

    for (ils_x, tag) in zip(ils_xs, tags):
        for par in pars:
            ilsn, ilsn_x = ilsfun.generate_ils_numerical(ils_x, ifglen, par)
            ilsa = ilsfun.norton_beer(ils_x, ifglen, par)
            assert np.allclose(ilsn_x, ils_x, atol=1e-5), f"ilsx, {tag}, {par}"
            assert np.allclose(ilsn, ilsa, atol=1e-5), f"ils, {tag}, {par}"
