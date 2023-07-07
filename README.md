# MPyLab
Pronounce as Am-p-lab

Python toolkit for reading, plotting, and analyzing micro-pulse lidar data.\
This toolkit took inspiration from peterkuma's project: https://github.com/peterkuma/mpl2nc which can read .mpl binary files and output netcdf files.\
This project is currently developed and maintained by Bo Chen and Dr. Sarah Brooks Lab at Texas A&M University, Atmospheric Sciences Department.

## pympl.py: Reading, Processing, and plotting .mpl data (**In development**)

- [x] Reading .MPL binary data files.
- [x] Reading afterpulse, overlap, and deadtime correction .bin files.
- [ ] Interpolate lidar data to even timesteps
- [ ] Data selection based on start-time and end-time.
- [ ] Method for making data and housekeeping data timeseries plot, single profile plot, and average profile plot.
- [ ] Calculate Copol and Crosspol NRB, as well as Depolarization Ratio.
- [ ] PBL, cloud, and aerosol layer detection based on the zero crossing algorithm.
- [ ] Fernald method for backscatter coefficient and extinction coefficient retrieval

## mpl_aerosol.py: Aerosol retrieval based on MPL data (**Planned**)

## T-matrix calculaiton dataset for lidar aerosol retrieval (**Planned**)

# PyMPL Documentation

The `PyMPL` class is a Python class that provides functionalities for reading, processing, and interpolating data from MPL (Micropulse Lidar) instruments. It allows you to work with MPL data files, perform corrections, and calculate various lidar products.

## Class Attributes

- `_C`: Speed of light in meters per second (constant).
- `_record_entry`: A list of record entry names representing different fields in the MPL data.
- `_record_format`: A list of record entry formats specifying the data type of each entry.
- `_byte_nums`: A dictionary mapping data types to the corresponding number of bytes.
- `_byte_format`: A dictionary mapping data types to the corresponding byte format specifier.
- `_ap_header`: A list of afterpulse header entry names.
- `_ap_header_format`: A list of afterpulse header entry formats.

## Constructor

The `PyMPL` class constructor initializes an instance of the class and takes the following parameters:

- `data_input`: Path to an MPL data file or a dictionary containing the data.
- `ap_input`: Path to the afterpulse correction file or a dictionary containing the afterpulse data.
- `ov_input`: Path to the overlap correction file or a dictionary containing the overlap data.
- `dt_input`: Path to the deadtime correction file or a dictionary containing the deadtime data.

If file paths are provided, the constructor automatically reads the corresponding files. If dictionaries are provided, they are used directly.

## Methods

### `read_files(cls, data_path, ap_path, ov_path, dt_path)`

This class method reads the MPL data, afterpulse correction, overlap correction, and deadtime correction files.

- `data_path`: Path to an MPL data file or folder containing multiple MPL data files.
- `ap_path`: Path to the afterpulse correction file.
- `ov_path`: Path to the overlap correction file.
- `dt_path`: Path to the deadtime correction file.

### `read_mpl_single_file(cls, file_path)`

This class method reads a single MPL data file.

- `file_path`: Path to the MPL data file.

### `read_mpl_multiple_file(cls, file_path_list)`

This class method reads multiple MPL data files and combines their data.

- `file_path_list`: A list of paths to MPL data files.

### `read_mpl(cls, file_path)`

This class method reads either a single MPL data file or all MPL data files in a folder.

- `file_path`: Path to an MPL data file or folder containing multiple MPL data files.

### `read_afterpulse(cls, file_path, number_bins)`

This class method reads the afterpulse correction file.

- `file_path`: Path to the afterpulse correction file.
- `number_bins`: Number of bins in the MPL data.

### `read_overlap(cls, file_path)`

This class method reads the overlap correction file.

- `file_path`: Path to the overlap correction file.

### `read_deadtime(cls, file_path)`

This class method reads the deadtime correction file.

- `file_path`: Path to the deadtime correction file.

### `calculate_dtcf(self, raw_data)`

This method calculates the deadtime correction factor for the provided raw data.

- `raw_data`: Raw lidar data.

### `calculate_r2_corrected(self, raw_data, background)`

This method calculates the range-squared corrected data by subtracting the background from the raw data.

- `raw_data`: Raw lidar data.
- `background`: Background lidar data.

### `calculate_nrb(self, raw_data, background, ap_data, ap_background)`

This method calculates the normalized relative backscatter (NRB) by applying deadtime, background, and afterpulse corrections.

- `raw_data`: Raw lidar data.
- `background`: Background lidar data.
- `ap_data`: Afterpulse correction data.
- `ap_background`: Afterpulse background data.

### `select_time(self, start_time, end_time)`

This method returns a new `PyMPL` object with data within the specified time range.

- `start_time`: Start time for the selected data.
- `end_time`: End time for the selected data.

### `interpolate_single_data(self, new_time_array, data, time_resolution, gap_seconds=None)`

This method interpolates a single lidar data variable to a new time array.

- `new_time_array`: Array of new time values.
- `data`: Lidar data to be interpolated.
- `time_resolution`: Time resolution in seconds.
- `gap_seconds`: Maximum gap in timestamps for interpolation.

### `interpolate_data(self, start_time, end_time, time_resolution, gap_seconds=None)`

This method interpolates lidar data within a specified time range.

- `start_time`: Start time for the interpolation.
- `end_time`: End time for the interpolation.
- `time_resolution`: Time resolution in seconds.
- `gap_seconds`: Maximum gap in timestamps for interpolation.

## Attributes

The `PyMPL` class has the following attributes:

- `datetime`: Array of timestamps for the lidar data.
- `seconds_since_start`: Array of time differences in seconds from the start of data collection.
- `energy_monitor`: Array of energy monitor values.
- `temp_detector`: Array of detector temperature values.
- `temp_telescope`: Array of telescope temperature values.
- `temp_laser`: Array of laser temperature values.
- `number_bins`: Number of lidar bins.
- `bin_time`: Time duration per lidar bin.
- `range`: Array of lidar range values.
- `bin_resolution`: Resolution of lidar bins.
- `range_edges`: Array of lidar range values for the edges of bins.
- `raw_copol`: Array of raw copolarized lidar data.
- `raw_crosspol`: Array of raw cross-polarized lidar data.
- `r2_corrected_copol`: Array of range-squared corrected copolarized lidar data.
- `r2_corrected_crosspol`: Array of range-squared corrected cross-polarized lidar data.
- `nrb_copol`: Array of normalized relative backscatter (NRB) copolarized lidar data.
- `nrb_crosspol`: Array of NRB cross-polarized lidar data.

## Interpolation

The `PyMPL` class provides interpolation capabilities for lidar data. You can use the `interpolate_data` method to interpolate data within a specified time range. The interpolated data can be accessed using the following attributes:

- `interpolated_datetime`: Array of timestamps for the interpolated data.
- `interpolated_seconds_since_start`: Array of time differences in seconds from the start of the interpolated data.
- `interpolated_raw_copol`: Interpolated raw copolarized lidar data.
- `interpolated_raw_crosspol`: Interpolated raw cross-polarized lidar data.
- `interpolated_r2_corrected_copol`: Interpolated range-squared corrected copolarized lidar data.
- `interpolated_r2_corrected_crosspol`: Interpolated range-squared corrected cross-polarized lidar data.
- `interpolated_nrb_copol`: Interpolated normalized relative backscatter (NRB) copolarized lidar data.
- `interpolatednrb_crosspol`: Interpolated NRB cross-polarized lidar data.

Note: The interpolation is performed using linear interpolation, and any gaps in the data exceeding the specified `gap_seconds` parameter will not be interpolated.

