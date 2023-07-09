<img src="assets/title.png" align="top" style="width: 100%;" />
<img src="assets/TAM-LogoBox.png" align="right" height = "90" />

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-3916/)

# MPyLab
Python toolkit for reading, plotting, and analyzing Micropulse Lidar data. Pronounce as Am-p-lab

This toolkit took inspiration from peterkuma's project: https://github.com/peterkuma/mpl2nc which can read .mpl binary files and output netcdf files.

This project is currently developed and maintained by Bo Chen of Dr. Sarah Brooks Lab at Texas A&M University, Atmospheric Sciences Department.

## Development Plan:
Note: the code is still in development and not well tested.

- pympl.py: Reading, Processing, and Interpolating .mpl data (**In development**)
  - [x] Reading .MPL binary data files.
  - [x] Reading afterpulse, overlap, and deadtime correction .bin files.
  - [x] Function for interpolating lidar data to even timesteps and deals with gaps in timestamps.
  - [x] Data selection based on start-time and end-time.
  - [x] Calculate SNR, R2 corrected, NRB, and Depol Ratio
  - [ ] PBL, cloud, and aerosol layer detection based on the zero crossing algorithm.
  - [ ] Fernald method for backscatter coefficient and extinction coefficient retrieval
    
- plotmpl.py: Ploting mpl data
  - [x] Function for 2d timeseries plot
  - [ ] Function for 1d timeseries plot
  - [ ] Function for single and average profile plot

- mpl_aerosol.py: Aerosol profile retrieval based on MPL data (**Planned**)
- T-matrix calculaiton dataset for lidar aerosol retrieval (**Planned**)

## Documentation:

- [Wikipage](https://github.com/bochens/MPyLab/wiki)

#### Micropulse Lidar data reading and processing
- [`PyMPL` class](https://github.com/bochens/MPyLab/wiki/pympl.PyMPL-Class)

#### Micropulse Lidar data visualization
- [`plotmpl` module](https://github.com/bochens/MPyLab/wiki/plotmpl)

