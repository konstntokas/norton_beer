"""
Description: This file contains the Norton-Beer apodization.

:copyright: 2023 (for authors see AUTHORS file). All rights reserved.
"""

import numpy as np
import logging


LOG = logging.getLogger(__name__)

"""
Predefined parameter for Norton-Beer; Ref: DOI: 10.1364/JOSAA.24.003644;
1.2, 1.4 and 1.6 correspond to the original weak, medium and strong modes
published by Norton and Beer 1976 (DOI: 10.1364/JOSA.66.000259)
"""
_NORTON_BEER_PARAMS = {
    1.0: np.array([1., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    1.1: np.array([0.701551, -0.639244, 0.937693, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    1.2: np.array([0.396430, -0.150902, 0.754472, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    1.3: np.array([0.237413, -0.065285, 0.827872, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    1.4: np.array([0.153945, -0.141765, 0.987820, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    1.5: np.array([0.077112, 0.0, 0.703371, 0.0, 0.219517, 0.0, 0.0, 0.0, 0.0]),
    1.6: np.array([0.039234, 0.0, 0.630268, 0.0, 0.234934, 0.0, 0.095563, 0.0, 0.0]),
    1.7: np.array([0.020078, 0.0, 0.480667, 0.0, 0.386409, 0.0, 0.112845, 0.0, 0.0]),
    1.8: np.array([0.010172, 0.0, 0.344429, 0.0, 0.451817, 0.0, 0.193580, 0.0, 0.0]),
    1.9: np.array([0.004773, 0.0, 0.232473, 0.0, 0.464562, 0.0, 0.298191, 0.0, 0.0]),
    2.0: np.array([0.002267, 0.0, 0.140412, 0.0, 0.487172, 0.0, 0.256200, 0.0, 0.113948])
}


def norton_beer(N, par):
    """ This function generates the Norton-Beer apodization

    Parameters
    ----------
    N : int
        number of samples
    par : float or 1darrray
        parameters of Norton-Beer apodization; if float is given,
        it must be one of [1.0, 1.1, ..., 2.0] to use the predefined
        parameters from <_NORTON_BEER_PARAMS>

    Returns
    -------
    apo : 1darray
        Norton-Beer apodization window

    Raises
    ------
    ValueError
        <par> for norton_beer apodization needs to be an 1darray with sum 1
        or one of the scalar in the following list [1.1, 1.2, ..., 2.0]

    Ref.:
    -----
    Norton and Beer 1976 https://doi.org/10.1364/JOSA.66.000259
    Naylor and Tahic 2007 https://doi.org/10.1364/JOSAA.24.003644
    """
    par = check_input(par)
    N_half = N // 2
    if N % 2 == 0:
        normv = 1.0 - ((np.arange(-N_half, N_half) + 0.5) / (N_half-0.5)) ** 2
    else:
        normv = 1.0 - (np.arange(-N_half, N_half+1) / N_half) ** 2
    apo = sum((a * (normv ** k) for (k, a) in enumerate(par)))
    return apo


def check_input(par):
    """ auxiliary function to assign FWHM to parameters, if given

    Parameters
    ----------
    par : float or 1darrray
        parameters of Norton-Beer apodization; if float is given,
        it must be one of [1.0, 1.1, ..., 2.0] to use the predefined
        parameters from <_NORTON_BEER_PARAMS>

    Returns
    -------
    par : 1darray
        parameters of Norton-Beer apodization

    Raises
    ------
    ValueError
        <par> for norton_beer apodization needs to one of
        the scalar in the following list [1.1, 1.2, ..., 2.0]
    ValueError
        "<par> for norton_beer apodization needs to be an 1darray with sum 1"
    """
    if np.isscalar(par):
        if np.any(np.isclose(par, list(_NORTON_BEER_PARAMS.keys()))):
            par = _NORTON_BEER_PARAMS[par]
        elif par == 1:
            par = [par]
        else:
            LOG.error("<par> for norton_beer apodization needs to one of "
                      "the scalar in the following list [1.1, 1.2, ..., 2.0]")
    else:
        if not np.isclose(np.sum(par), 1, atol=1e-3):
            LOG.error("<par> for norton_beer apodization needs to be an 1darray with sum 1")
    return par
