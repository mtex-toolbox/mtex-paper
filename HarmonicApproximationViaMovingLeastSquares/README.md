# Fast and Stable Harmonic Approximation from Ill-distributed Data using Moving Least Squares — experiment code

R. Hielscher, T. Pöschl, E. Wünsche, TU Bergakademie Freiberg

This folder contains the MATLAB scripts that produce the figures, the table and
the numbers quoted in Section 3 of the paper (numerical experiments on the
sphere), together with their input data and the results shown in the paper.

## Requirements

- MATLAB; the scripts use `exportgraphics` and `writematrix`, which are
  available from R2020a on.
- [MTEX](https://mtex-toolbox.github.io) 7.1.0, which provides the HAMLS
  implementation (`S2FunMLS`, `S2FunHarmonic`) as well as the NFFT-based
  spherical Fourier transforms and the plotting routines used here. The
  experiments of the paper used FFTW 3.3.10 and NFFT 3.5.3. Earlier MTEX
  releases calibrate the regularization of the local MLS problems differently
  and give slightly different results.
- Only for `weather/prepareWeatherData.m`: an internet connection, `gunzip`, and
  the Mapping Toolbox (for the coast lines).

## Contents

```
weather/                    Section 3.1
  prepareWeatherData.m      downloads the GHCN data and writes data/
  weatherExample.m          HAMLS and LSQR on the weather data
  data/                     station data, coast lines
  results/                  outputs shown in the paper
synthetic/                  Section 3.2
  errorVsBandwidth.m        test function, error and runtime over the bandwidth
  equalTimeError.m          LSQR at multiples of the HAMLS runtime
  stepRuntimes.m            runtime of Steps 2 and 3 of Algorithm 1
  data/                     the test function
  results/errorVsBandwidth/
  results/equalTime/
  results/stepRuntimes/
```

Where each item of the paper comes from:

| Paper | Script | File(s) in `results/` |
| --- | --- | --- |
| Figure 1 | `weather/weatherExample.m` | `WeatherData.png`, `MLSApproximationWeatherData.png`, `HarmMLSApproximationWeatherData.png`, `LSQRApproximationWeatherData_{1,5,10,30}xT.png`, `ColorbarWeather.png` |
| Figure 2 | `weather/weatherExample.m` | `IterationStudy_reg_{1em07,1em08,1em09,1em10}.txt` |
| runtimes and iteration counts in Section 3.1 | `weather/weatherExample.m` | `LSQRPanels.txt` |
| Figure 3 | `synthetic/errorVsBandwidth.m` | `errorVsBandwidth/ToyFun.png`, `ToyData.png`, `ColorbarToy.eps` |
| Figure 4 | `synthetic/errorVsBandwidth.m` | `errorVsBandwidth/ApproximationError.txt` |
| runtimes in Section 3.2 | `synthetic/errorVsBandwidth.m` | `errorVsBandwidth/ApproximationTime.txt`, `LSQRError.mat` |
| Table 1 | `synthetic/stepRuntimes.m` | `stepRuntimes/RuntimeMLSQuadrature.txt` |
| Figure 5 | `synthetic/equalTimeError.m` | `equalTime/LSQRErrorTimeFactors.txt` |

Figures 2, 4 and 5 and Table 1 are drawn with pgfplots in the LaTeX source of
the paper, directly from these tab-separated tables; the column names are given
in the first line of each file. The `.mat` files hold the individual
measurements behind the tables and the intermediate results that later scripts
read.

## Running the scripts

The scripts are cell-mode scripts. Run their sections in order, since each
section uses the variables of the previous ones, and run them from within their
own folder, `weather/` or `synthetic/`, since all paths are relative to the
current working directory.

In `synthetic/`, run `errorVsBandwidth.m` first. It writes the data set and the
regularization parameters that the two other scripts depend on.

Each of the four experiment scripts has a summary section that prints the
numbers quoted in the text and collects them in a struct `paperNumbers`.

A run overwrites the files in the corresponding `results/` folder; make a copy
first to keep the ones of the paper. Most of the runtime is spent in LSQR. On
the Ryzen 7 5800X used in the paper, runtimes were approximately as follows:

- `weatherExample.m`: 2 hours (40 min for the iteration study of Figure 2, and
  1 hour for the LSQR panels)
- `errorVsBandwidth.m`: 3 hours for the LSQR, much longer for the search of the
  error-minimizing regularization parameter at every bandwidth
- `equalTimeError.m`: 2 days, almost all of it in the sweep over the time
  budgets. The script prints an estimate before this sweep starts.
- `stepRuntimes.m`: well under 1 hour

## Reproducibility

The random nodes of the synthetic example are drawn with a fixed seed, and
`stepRuntimes.m` checks that it redraws the data set saved by
`errorVsBandwidth.m`. The errors of Figure 4 therefore reproduce up to
rounding. For HAMLS this means agreement to about 13 digits. LSQR, however,
amplifies rounding over its thousands of iterations, and the rounding depends on
the number of threads; at the largest bandwidths its errors can therefore differ
from those of the paper in the second or third digit.

Runtimes depend on the machine; those of the paper were measured on the machine
described in Section 3. This affects not only the reported runtimes: the LSQR
iteration budgets of Figure 1d–g and of Figure 5 are derived from measured
runtimes, and so are the errors LSQR reaches within them. For the weather
example, setting `paperCounts = true` in `weatherExample.m` redraws the LSQR
panels at the iteration counts of the paper, listed in
`weather/results/LSQRPanels.txt`.

## Data

- `weather/data/Weather_Data.mat`: daily mean temperature (TAVG) of 21 March
  2018 at 6501 stations of the Global Historical Climatology Network – Daily,
  written by `prepareWeatherData.m`. Source: M. J. Menne et al., Global
  Historical Climatology Network – Daily (GHCN-Daily), Version 3, NOAA National
  Centers for Environmental Information, <https://doi.org/10.7289/V5D21VHZ>.
- `weather/data/coastLines.mat`: coast lines for the plots, from the Natural
  Earth 1:50m land polygons in `weather/data/land_contours/` (public domain,
  <https://www.naturalearthdata.com>).
- `synthetic/data/ToyExampleFun.mat`: the synthetic test function, an
  `S2FunHarmonic` of bandwidth 512. The scripts rescale it and add 300 to obtain
  the function f₀ of Section 3.2.
