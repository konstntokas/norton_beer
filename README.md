# norton_beer

The 'norton_beer' library is a python tool box for the Norton-Beer apodization functions. The Norton-Beer apodization is firstly descibed by Norton and Beer (1976), [doi:10.1364/JOSA.66.000259](https://www.osapublishing.org/abstract.cfm?URI=josa-66-3-259) and Norton and Beer (1977), [10.1364/JOSA.67.000419](https://opg.optica.org/abstract.cfm?URI=josa-67-3-419). It is further extended by Naylor and Tahic (2007), [10.1364/JOSAA.24.003644](https://www.osapublishing.org/abstract.cfm?URI=josaa-24-11-3644)

## Main features

- The apodization window in spatial domain can be generated for a given set of parameters
- The analytical Fourier transform of the apodization function can be generated for a given set of parameters
- New apodization functions can be generated, where the parameters are optimized so that for a given full with at half maximum (FWHM) the absolute maximal value of the side lobes is minimized
