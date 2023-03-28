import numpy as np
import math
import norton_beer.apodization as apo
from norton_beer.apodization import check_input


def norton_beer_old(k, ifglen, par):
    """ This function creates the analytical Fourier Transform
        of the Norton-Beer apodization

    Parameters
    ----------
    k : 1darray
        spectral axis
    ifglen : float
        signal length
    par : float or 1darray
        parameters of Norton-Beer apodization

    Returns
    -------
    result : 1darray
        Fourier transformation of Norton-Beer apodization window

    Ref.:
    -----
    Norton and Beer 1976 https://doi.org/10.1364/JOSA.66.000259
    Naylor and Tahic 2007 https://doi.org/10.1364/JOSAA.24.003644
    """

    def calc_q_fac(N):
        sinc_fac = np.empty((N, N))
        cos_fac = np.empty((N, N))

        # start value
        sinc_start = np.zeros(N)
        sinc_start[0] = 1
        sinc_fac[0, :] = sinc_start
        cos_fac[0, :] = np.zeros(N)
        # recursive calcuation of the factors
        for m in range(1, N):
            sinc_line = np.zeros(N)
            sinc_line[0] = 1
            sinc_line[1:] = 2*m * (2*m-1) * sinc_fac[m-1, :-1]
            sinc_fac[m, :] = sinc_line

            cos_line = np.zeros(N)
            cos_line[1] = 2*m
            cos_line[2:] = 2*m * (2*m-1) * cos_fac[m-1, 1:-1]
            cos_fac[m, :] = cos_line

        # account for negative values
        fac_neg = np.ones(N)
        for n in range(N):
            if n % 2 != 0:
                fac_neg[n] = -1
        fac_neg = fac_neg.reshape(1, -1)
        sinc_fac *= fac_neg
        cos_fac *= fac_neg * (-1)

        fac_sinc_mx = np.empty((N, N))
        fac_cos_mx = np.empty((N, N))
        bound = max(N, 5)
        fac_taylor_mx = np.empty((N, bound))
        apow_taylor = 2 * np.arange(bound)
        for n in range(N):
            # multiply rows with factors
            fac = np.array([math.comb(n, m) * (-1)**m for m in range(n+1)])
            fac = fac.reshape(-1, 1)
            sinc_fac_n = sinc_fac[:n+1, :] * fac
            cos_fac_n = cos_fac[:n+1, :] * fac
            sinc_fac_n = np.sum(sinc_fac_n, axis=0)
            cos_fac_n = np.sum(cos_fac_n, axis=0)
            fac_sinc_mx[n, :] = sinc_fac_n
            fac_cos_mx[n, :] = cos_fac_n

            # apply Taylor expansion
            K = 2 * N
            j_sum_idx = np.tile(np.arange(K), N)
            i_sum_idx = np.repeat(np.arange(N), K)

            fac_taylor = np.empty(len(i_sum_idx))
            fac_apow = np.empty(len(i_sum_idx), dtype=int)
            for i, (i_idx, j_idx) in enumerate(zip(i_sum_idx, j_sum_idx)):
                fac_taylor[i] = ((-1)**j_idx * (sinc_fac_n[i_idx] / math.factorial(2 * j_idx + 1)
                                 + cos_fac_n[i_idx] / math.factorial(2 * j_idx)))
                fac_apow[i] = int(2 * j_idx - 2 * i_idx)

            fac_taylor_sort = np.array([])
            for i in apow_taylor:
                fac_sum = np.sum(fac_taylor[fac_apow == i])
                fac_taylor_sort = np.append(fac_taylor_sort, fac_sum)
            fac_taylor_mx[n, :] = fac_taylor_sort

        return fac_sinc_mx, fac_cos_mx, fac_taylor_mx

    result = np.empty_like(k)
    a = np.pi * k * ifglen
    par = check_input(par)
    par = np.trim_zeros(par, trim="b")

    # get factors
    sinc_mx, cos_mx, taylor_mx = calc_q_fac(len(par))

    # case 1: if a < 1 -> analytical solution unstable -> expansion
    mask = abs(a) < 1
    apower = np.empty((taylor_mx.shape[1], len(a[mask])))
    apower[0, :] = 1
    apower[1, :] = a[mask]**2
    for i in range(2, apower.shape[0]):
        apower[i, :] = apower[i-1, :] * apower[1, :]
    qns = np.matmul(taylor_mx, apower)
    result[mask] = np.sum(par.reshape(-1, 1) * qns, axis=0)

    # case 2: if a >= 1 -> analytical solution
    sinca, cosa = np.sinc(a[~mask] / np.pi), np.cos(a[~mask])
    apower = np.empty((sinc_mx.shape[1], len(a[~mask])))
    apower[0, :] = np.ones(apower.shape[1])
    apower[1, :] = a[~mask]**2
    for i in range(2, apower.shape[0]):
        apower[i, :] = apower[i-1, :] * apower[1, :]
    apower_inv = 1 / apower
    qns = np.matmul(sinc_mx, apower_inv) * sinca + np.matmul(cos_mx, apower_inv) * cosa
    result[~mask] = np.sum(par.reshape(-1, 1) * qns, axis=0)

    return result


def norton_beer(k, ifglen, par):
    """ This function creates the analytical Fourier Transform
        of the Norton-Beer apodization

    Parameters
    ----------
    k : 1darray
        spectral axis
    ifglen : float
        signal length
    par : float or 1darray
        parameters of Norton-Beer apodization

    Returns
    -------
    result : 1darray
        Fourier transformation of Norton-Beer apodization window

    Ref.:
    -----
    Norton and Beer 1976 https://doi.org/10.1364/JOSA.66.000259
    Naylor and Tahic 2007 https://doi.org/10.1364/JOSAA.24.003644
    """

    def calc_q_fac(N):
        sinc_fac = np.zeros((N, N))
        sinc_fac[0, 0] = 1
        cos_fac = np.zeros((N, N))

        # recursive calcuation of the factors
        for m in range(1, N):
            sinc_fac[m, 0] = 1
            sinc_fac[m, 1:] = 2*m * (2*m-1) * sinc_fac[m-1, :-1]

            cos_fac[m, 1] = 2*m
            cos_fac[m, 2:] = 2*m * (2*m-1) * cos_fac[m-1, 1:-1]

        # account for negative values
        fac_neg = np.ones(N)
        fac_neg[::2] = -1
        sinc_fac *= -fac_neg[np.newaxis, :]
        cos_fac *= fac_neg[np.newaxis, :]

        # calculate outer sum with binomial coeff.
        fac_sinc_mx = np.zeros((N, N))
        fac_cos_mx = np.zeros((N, N))
        for n in range(N):
            fac = np.array([math.comb(n, m) * (-1)**m for m in range(n+1)])[:, np.newaxis]
            fac_sinc_mx[n, :] = np.sum(sinc_fac[:n+1, :] * fac, axis=0)
            fac_cos_mx[n, :] = np.sum(cos_fac[:n+1, :] * fac, axis=0)

        # calculate factor when taylor expansion is applied
        bound = max(N, 5)
        K = 2 * (bound - 1)
        idx0 = np.tile(np.arange(K), N)
        idx1 = np.repeat(np.arange(N), K)
        factorial0 = np.array([math.factorial(2 * idx) for idx in idx0]).astype(float)
        factorial1 = np.array([math.factorial(2 * idx + 1) for idx in idx0]).astype(float)
        fac_taylor = ((-1)**idx0[np.newaxis, :] * (fac_sinc_mx[:, idx1] / factorial1[np.newaxis, :]
                      + fac_cos_mx[:, idx1] / factorial0[np.newaxis, :]))
        fac_apow = (2 * idx0 - 2 * idx1).astype(int)

        fac_apow_uni = np.unique(fac_apow)
        for pow in fac_apow_uni[fac_apow_uni < 0]:
            assert np.allclose(np.sum(fac_taylor[:, fac_apow == pow], axis=1), 0.)
        fac_taylor_mx = np.zeros((N, bound))
        for i, pow in enumerate(2 * np.arange(bound)):
            fac_taylor_mx[:, i] = np.sum(fac_taylor[:, fac_apow == pow], axis=1)

        return fac_sinc_mx, fac_cos_mx, fac_taylor_mx

    ils = np.zeros_like(k)
    a = np.pi * k * ifglen
    par = check_input(par)
    par = np.trim_zeros(par, trim="b")

    # get factors
    sinc_mx, cos_mx, taylor_mx = calc_q_fac(len(par))

    # case 1: if a < 1 -> analytical solution unstable -> expansion
    mask = abs(a) < 1
    apower = np.zeros((taylor_mx.shape[1], len(a[mask])))
    apower[0, :] = 1
    apower[1, :] = a[mask]**2
    for i in range(2, apower.shape[0]):
        apower[i, :] = apower[i-1, :] * apower[1, :]
    qns = np.matmul(taylor_mx, apower)
    ils[mask] = np.sum(par[:, np.newaxis] * qns, axis=0)

    # case 2: if a >= 1 -> analytical solution
    sinca, cosa = np.sinc(a[~mask] / np.pi), np.cos(a[~mask])
    apower = np.zeros((sinc_mx.shape[1], len(a[~mask])))
    apower[0, :] = 1
    if apower.shape[0] > 1:
        apower[1, :] = a[~mask]**2
    if apower.shape[0] > 2:
        for i in range(2, apower.shape[0]):
            apower[i, :] = apower[i-1, :] * apower[1, :]
    apower_inv = 1 / apower
    qns = np.matmul(sinc_mx, apower_inv) * sinca + np.matmul(cos_mx, apower_inv) * cosa
    ils[~mask] = np.sum(par[:, np.newaxis] * qns, axis=0)

    return ils


def generate_ils_numerical(k, ifglen, par, nb_sample=100001):
    """ This function generates the ILS numerically via discrete Fourier Transform

    Parameters
    ----------
    k : 1darray
        spectral axis
    ifglen : float
        signal length
    par : float or 1darray
        parameters of Norton-Beer apodization
    nb_samples : int, optional
        number of samples, should be odd, by default 100001

    Returns
    -------
    ils : 1darray
        numerical Fourier transformation of
        Norton-Beer apodization window
    ils_x : 1darray
        corresponding spectral axis; should
        be close to <k>
    """

    par = apo.check_input(par)

    # many samples such that error between
    # numerical and analytical ILS decreases
    k_step = round(abs(k[0] - k[1]), 10)
    dx = ifglen / nb_sample
    n_pad = int(round(1 / (k_step * dx), 0))
    assert n_pad >= nb_sample

    # generate ILS via FFT
    apo_spat = apo.norton_beer(nb_sample, par)

    # zero padding
    if (n_pad % 2 != 0) & (nb_sample % 2 != 0):
        start = (n_pad - nb_sample) // 2
    elif (n_pad % 2 == 0) & (nb_sample % 2 != 0):
        start = (n_pad - nb_sample) // 2 + 1
    end = n_pad - start - nb_sample
    apo_spat = np.pad(apo_spat, (start, end))

    # get Fourier transform of apodization function
    apo_spat = np.fft.ifftshift(apo_spat, axes=-1)
    ils = np.fft.fft(apo_spat) / nb_sample
    ils_x = np.fft.fftfreq(n_pad, d=dx)
    ils = np.fft.fftshift(ils)
    ils_x = np.fft.fftshift(ils_x)
    assert not np.any(abs(np.imag(ils)) > 1e-10)
    ils = np.real(ils)

    # cut at range of original x-axis
    idx = (ils_x >= k[0] - k_step/2) & (ils_x <= k[-1] + k_step/2)
    ils = ils[idx]
    ils_x = ils_x[idx]

    return ils, ils_x
