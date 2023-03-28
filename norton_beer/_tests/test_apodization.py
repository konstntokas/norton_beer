"""
Description: This file tests the function in apodization.py.

:copyright: 2023 (for authors see AUTHORS file). All rights reserved.
"""

import norton_beer.apodization as apo
import numpy as np
import pickle
from pathlib import Path
import os


def test_apodization():
    """
    Tests the function apodization
    """

    _DIR = Path(__file__).resolve().parent
    files = ["spatial_even.pkl", "spatial_odd.pkl"]
    for file in files:
        with open(os.path.join(_DIR, file), "rb") as filehandler:
            dat = pickle.load(filehandler)
        for (k, v) in dat.items():
            res = apo.norton_beer(len(v), k)
            assert np.allclose(res, v), f"{file} {k}"
