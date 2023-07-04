# pyMPL
Python toolkit for reading, plotting, and analyzing micro-pulse lidar data.

This toolkit took inspiration from peterkuma's project: https://github.com/peterkuma/mpl2nc which can read .mpl binary files and output netcdf files.

## pyMPL.py: Reading, Processing, and plotting .mpl data (**In development**)

- [x] Reading .MPL binary data files.
- [x] Reading afterpulse, overlap, and deadtime correction .bin files.
- [ ] Data selection based on start-time and end-time.
- [ ] Method for making data and housekeeping data timeseries plot, single profile plot, and average profile plot.
- [ ] Calculate Copol and Crosspol NRB, as well as Depolarization Ratio.
- [ ] PBL, cloud, and aerosol layer detection based on the zero crossing algorithm.
- [ ] Fernald method for backscatter coefficient and extinction coefficient retrieval

## MPL_aerosol.py: Aerosol retrieval based on MPL data (**Planned**)

## T-matrix calculaiton dataset for lidar aerosol retrieval (**Planned**)
