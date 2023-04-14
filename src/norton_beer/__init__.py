"""
norton_beer
===========
The 'norton_beer' library is a python tool box for
the Norton-Beer apodization functions.

Provides:
    1. Norton-Beer apodization window generation in spatial domain
    2. the analytical Fourier transform of the apodization window
       also known as instrument line shape (ILS)
    3. generation of new apodization window for a fixed spatial
       resolution of the ILS so that the sides lobes are minimized

Documentation:
    1. add paper when published
    2. add webpage of sphinx built documentation

Installation
------------

To use norton_beer, first install it using pip:

.. code-block:: console

   (.venv) $ pip install norton_beer


Available modules
-----------------
apodization
    generation of apodization window in spatial domain
ils
    generation of ILS in spectral domain
optimize
    generation of new apodizations for a fixed spatial
    resolution of the ILS so that the sides lobes are minimized

Further references
------------------
.. [1] Norton, R. H., & Beer, R. (1976). New apodizing functions
    for Fourier spectrometry. Journal of the Optical Society of America,
    66 (3), 259. https://doi.org/10.1364/JOSA.66.000259
.. [2] Norton, R. H., & Beer, R. (1977). Errata: New Apodizing
    Functions For Fourier Spectrometry. Journal of the Optical Society
    of America, 67 (3), 419. https://doi.org/10.1364/JOSA.67.000419
.. [3] Naylor, D. A., & Tahic, M. K. (2007). Apodizing functions
    for Fourier transform spectroscopy. Journal of the Optical Society
    of America A, 24 (11), 3644.260 https://doi.org/10.1364/JOSAA.24.003644
"""
import pkg_resources

__author__ = "Konstantin Ntokas"
__authors__ = ["Konstantin Ntokas", "Jörn Ungermann"]
__copyright__ = "Copyright 2023, Forschungszentrum Juelich GmbH"
__license__ = "GNU General Public License v3 (GPLv3)"
__version__ = pkg_resources.get_distribution(__name__).version
