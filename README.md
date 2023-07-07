<img src="TAM-Logo-white.png" align="right" height = "90" />

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/pypi/pyversions/event-bus-py2.svg)](https://www.python.org/download/releases/3.9/)

# MPyLab
Pronounce as Am-p-lab

Python toolkit for reading, plotting, and analyzing Micropulse Lidar data.

This toolkit took inspiration from peterkuma's project: https://github.com/peterkuma/mpl2nc which can read .mpl binary files and output netcdf files.

This project is currently developed and maintained by Bo Chen of Dr. Sarah Brooks Lab at Texas A&M University, Atmospheric Sciences Department.

## Development Plan:

- pympl.py: Reading, Processing, and plotting .mpl data (**In development**)
  - [x] Reading .MPL binary data files.
  - [x] Reading afterpulse, overlap, and deadtime correction .bin files.
  - [ ] Interpolate lidar data to even timesteps
  - [ ] Data selection based on start-time and end-time.
  - [ ] Method for making data and housekeeping data timeseries plot, single profile plot, and average profile plot.
  - [ ] Calculate Copol and Crosspol NRB, as well as Depolarization Ratio.
  - [ ] PBL, cloud, and aerosol layer detection based on the zero crossing algorithm.
  - [ ] Fernald method for backscatter coefficient and extinction coefficient retrieval
- mpl_aerosol.py: Aerosol retrieval based on MPL data (**Planned**)
- T-matrix calculaiton dataset for lidar aerosol retrieval (**Planned**)

## Documentation:

- [Wikipage](https://github.com/bochens/MPyLab/wiki)
- [`PyMPL` class](https://github.com/bochens/MPyLab/wiki/PyMPL-Class)
