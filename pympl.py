import glob
import numpy as np
import os
import struct
import scipy
import copy
from netCDF4 import Dataset

class PyMPL:
    _C = 299792458 # m.s-1 speed of light

    _record_entry = ['unit_number','version','year','month','day','hours','minutes',
        'seconds','shots_sum','trigger_frequency','energy_monitor','temp_0',
        'temp_1','temp_2','temp_3','temp_4','background_average','background_std_dev',
        'number_channels','number_bins','bin_time','range_calibraiton','number_data_bins',
        'scan_scenario_flag','number_of_background_bins','azimuth_angle','elevation_angle',
        'campass_degrees','lidar_site','wavelength','gps_latitude','gps_longitude',
        'gps_altitude','ad_data_bad_flag','data_file_version','background_average_2',
        'background_std_dev_2','mcs_mode','first_data_bin','system_type',
        'sync_pulses_seen_per_second','first_background_bin','header_size',
        'ws_used','ws_inside_temperature','ws_outside_temperature',
        'ws_inside_humidity','ws_outside_humidity','ws_dew_point',
        'ws_wind_speed','ws_wind_direction','ws_baronmetric_pressure',
        'ws_rain_rate','channel_1_data','channel_2_data']
    
    _record_format = ['uint16','uint16','uint16','uint16','uint16','uint16',
                    'uint16','uint16','uint32','int32','uint32','uint32',
                    'uint32','uint32','uint32','uint32','float32','float32',
                    'uint16','uint32','float32','float32','uint16','uint16',
                    'uint16','float32','float32','float32','char*6','uint16',
                    'float32','float32','float32','byte','byte','float32',
                    'float32','byte','uint16','byte','uint16','uint16','uint16',
                    'byte','float32','float32','float32','float32','float32',
                    'float32','int16','float32','float32','float32 array',
                    'float32 array']

    _byte_nums = {'uint8':1,'uint16':2,'int16':2,'uint32':4,'int32':4,'float32':4,'float64':8,'byte':1,'char*6':6,'float32 array':4}
    _byte_format =  {'uint8':'<B','uint16':'<H','int16':'<h','uint32':'<I','int32':'<i','float32':'<f','float64':'<d','byte':'<B','char*6':'<cccccc',
                    'float32 array':'<'+'2000f'}
    
    _ap_header = ['ap_header','ap_file_version','ap_number_channels','ap_number_bins','ap_energy','ap_background_average_copol','ap_background_average_crosspol']
    _ap_header_format = ['uint32','uint16','uint8','uint32','float64','float64','float64']

    def __init__(self, data_input, ap_input, ov_input, dt_input, blind_range = 0.1):
        
        # make sure the instance constructor can take both dictionary and file path as input
        if all([isinstance(element, dict) for element in [data_input, ap_input, ov_input, dt_input]]):
            self.data_dict, self.ap_dict, self.ov_dict, self.dt_dict = data_input, ap_input, ov_input, dt_input
        else:
            self.data_dict, self.ap_dict, self.ov_dict, self.dt_dict = self.read_files(data_input, ap_input, ov_input, dt_input)
        
        self.blind_range = blind_range

        # prepare datetime array
        year_str    = np.char.zfill(np.array(self.data_dict['year']).astype(str), 4)
        month_str   = np.char.zfill(np.array(self.data_dict['month']).astype(str), 2)
        day_str     = np.char.zfill(np.array(self.data_dict['day']).astype(str), 2)
        hours_str   = np.char.zfill(np.array(self.data_dict['hours']).astype(str), 2)
        minutes_str = np.char.zfill(np.array(self.data_dict['minutes']).astype(str), 2)
        seconds_str = np.char.zfill(np.array(self.data_dict['seconds']).astype(str), 2)
        self.datetime = np.array([x+'-'+y+'-'+z+'T'+h+':'+m+':'+s \
                                  for x,y,z,h,m,s in zip(year_str,month_str,day_str,hours_str,minutes_str,seconds_str)], dtype='datetime64')
        self.seconds_since_start = (self.datetime-self.datetime[0]).astype('timedelta64[s]').astype(int) # Seconds since first data entry

        self.laser_energy   = np.array(self.data_dict['energy_monitor'])/1000   # micro Joules (uJ)
        self.temp_detector  = np.array(self.data_dict['temp_0'])/100            # degree C
        self.temp_telescope = np.array(self.data_dict['temp_2'])/100            # degree C
        self.temp_laser     = np.array(self.data_dict['temp_3'])/100            # degree C

        self.number_profile = self.datetime.size

        self.number_bins = self.data_dict['number_bins'][0]
        if np.all(np.array(self.data_dict['number_bins']) != self.number_bins):
            raise ValueError('number_bins of scans have to be the same')
        
        self.bin_time = self.data_dict['bin_time'][0] # bin width in seconds
        if np.all(np.array(self.data_dict['bin_time']) != self.bin_time):
            raise ValueError('bin_time of scans have to be the same')
        
        background_copol    = np.array(self.data_dict['background_average_2'])  # count per micro seconds (count us-1)
        background_crosspol = np.array(self.data_dict['background_average'])    # count per micro seconds (count us-1)

        #self.no_blind_range      = 0.5*self.bin_time*self._C*(np.arange(self.number_bins) + 1)*1e-3 #km
        # # a previous version uses +0.5 following Peter Kuma. However, using +1 generates the same r2corrected result as the sigmaMPL software
        self.no_blind_range      = 0.5*self.bin_time*self._C*(np.arange(self.number_bins) + 0.5)*1e-3 #km
        above_blind_range   = self.no_blind_range > self.blind_range

        self.range          = self.no_blind_range[above_blind_range] # do not include range below blind range
        self.bin_resolition = self.range[1]-self.range[0]
        self.range_edges    = np.append(self.range-self.bin_resolition/2, self.range[-1]+self.bin_resolition/2)

        self.raw_copol      = np.array(self.data_dict['channel_2_data'])[:,above_blind_range] # count per micro seconds (count us-1)
        self.raw_crosspol   = np.array(self.data_dict['channel_1_data'])[:,above_blind_range] # count per micro seconds (count us-1)

        self.snr_copol      = self.calculate_snr(self.raw_copol, background_copol, np.array(self.data_dict['background_std_dev_2']))     # unitless
        self.snr_crosspol   = self.calculate_snr(self.raw_crosspol, background_crosspol, np.array(self.data_dict['background_std_dev'])) # unitless
        ## Calculated product ##
        # Calculate Range Correction Product
        self.r2_corrected_copol    = self.calculate_r2_corrected(self.raw_copol, background_copol)        # counts km2 per micro seconds (counts km2 us-1)
        self.r2_corrected_crosspol = self.calculate_r2_corrected(self.raw_crosspol, background_crosspol)  # counts km2 per micro seconds (counts km2 us-1)

        # Calculate Normalized Relative Backscatter
        self.nrb_copol    = self.calculate_nrb(self.raw_copol,    background_copol    \
                                               ,self.ap_dict['ap_copol'],self.ap_dict['ap_background_average_copol'])       # counts km2 per micro seconds per micro Joules (counts km2 us-1 uJ-1)
        self.nrb_crosspol = self.calculate_nrb(self.raw_crosspol, background_crosspol \
                                               ,self.ap_dict['ap_crosspol'],self.ap_dict['ap_background_average_crosspol']) # counts km2 per micro seconds per micro Joules (counts km2 us-1 uJ-1)
        self.depol_ratio  = self.calculate_depol_ratio()    # unitless

        ## Interpolation product ##
        self.interpolation_flag                   = False   # True if the data has been interpolated
        self.interpolated_datetime                = None    # Use this as datetime for plotting interpolated data
        self.interpolated_laser_energy            = None    
        self.interpolated_temp_telescope          = None
        self.interpolated_temp_detector           = None
        self.interpolated_temp_laser              = None
        self.interpolated_seconds_since_start     = None
        self.interpolated_raw_copol               = None
        self.interpolated_raw_crosspol            = None
        self.interpolated_r2_corrected_copol      = None
        self.interpolated_r2_corrected_crosspol   = None
        self.interpolated_nrb_copol               = None
        self.interpolated_nrb_crosspol            = None
        self.interpolated_depol_ratio             = None
        self.interpolated_snr_copol               = None
        self.interpolated_snr_crosspol            = None

    def deepcopy(self):
        return copy.deepcopy(self)
        
    @classmethod
    def read_files(cls, data_path, ap_path, ov_path, dt_path):
        '''
        data_path: Path to mpl data or mpl data folder. Could also be a list of mpl data path.
        ap_path: Path to afterpulse correction file
        ov_path: Path to overlap correction file
        dt_path: Path to deadtime correction file
        '''
        # Read MPL data
        data_dict, number_bins = cls.read_mpl(data_path)            # Read MPL data
        ap_dict   = cls.read_afterpulse(ap_path, number_bins)       # Read Afterpulse data
        ov_dict   = cls.read_overlap(ov_path)                       # Read Overlap data
        dt_dict   = cls.read_deadtime(dt_path)                      # Read Deadtime data
        return data_dict, ap_dict, ov_dict, dt_dict

    @classmethod
    def read_mpl_single_file(cls, file_path):
        fin = open(file_path, "rb")
        data_dict = {}
        for i in range(len(cls._record_entry)):
            data_dict[cls._record_entry[i]] = []
        loop_bool = True
        number_bins = None
        while loop_bool:
            for i in range(len(cls._record_entry)):
                entry_type = cls._record_entry[i]
                entry_format = cls._record_format[i]
                if entry_type == 'channel_1_data' or entry_type == 'channel_2_data': 
                    byte_numbers_to_read = cls._byte_nums[entry_format] * number_bins
                else:
                    byte_numbers_to_read = cls._byte_nums[entry_format]
                binary_entry =  fin.read(byte_numbers_to_read)
                if binary_entry == b'': # check if reached the end of the file, if so then stop the loop and do not process
                    loop_bool = False
                else:
                    if entry_format in ['uint16','uint32','int16','int32','float32']:
                        unpacked_entry = struct.unpack(cls._byte_format[entry_format], binary_entry)
                        data_dict[entry_type].append(unpacked_entry[0])
                        if entry_type == 'number_bins':
                            number_bins = unpacked_entry[0]
                    elif entry_format == 'byte':
                        data_dict[entry_type].append(bin(int.from_bytes(binary_entry, byteorder='little'))[2:].zfill(8))
                    elif entry_format == 'char*6':
                        data_dict[entry_type].append(binary_entry.decode('latin-1'))
                    elif entry_format == 'float32 array':
                        unpacked_entry = struct.unpack(cls._byte_format[entry_format], binary_entry)
                        data_dict[entry_type].append(list(unpacked_entry))
        fin.close()
        return data_dict, number_bins

    @classmethod
    def read_mpl_multiple_file(cls, file_path_list):
        data_dict = {}
        for i in range(len(cls._record_entry)):
            data_dict[cls._record_entry[i]] = []
        for a_path in file_path_list:
            data_dict_to_merge, number_bins = cls.read_mpl_single_file(a_path)
            for i in range(len(cls._record_entry)):
                data_dict[cls._record_entry[i]].extend(data_dict_to_merge[cls._record_entry[i]])
        return data_dict, number_bins

    @classmethod
    def read_mpl(cls, file_path):
        if isinstance(file_path, str):
            if os.path.isdir(file_path):
                # read all mpl files contained in a folder
                filelist = glob.glob(os.path.join(file_path, '*.mpl'))
                filelist.sort(key=os.path.getmtime)
                return cls.read_mpl_multiple_file(filelist)
            elif os.path.isfile(file_path):
                # read a single mpl file
                return cls.read_mpl_single_file(file_path)
        elif isinstance(file_path, list):
            # read all mpl files from paths contianed in a list
            return cls.read_mpl_multiple_file(file_path)
    
    @classmethod
    def read_afterpulse(cls, file_path, number_bins):
        fin = open(file_path, "rb")
        ap_dict = {}
        for i in range(len(cls._ap_header)):
            entry_type = cls._ap_header[i]
            entry_format = cls._ap_header_format[i]
            byte_numbers_to_read = cls._byte_nums[entry_format]
            binary_entry =  fin.read(byte_numbers_to_read)
            ap_dict[entry_type] = struct.unpack(cls._byte_format[entry_format], binary_entry)[0]
        ap_number_bins = int(ap_dict['ap_number_bins'])
        for x in ['ap_range', 'ap_copol', 'ap_crosspol']:
            buf = fin.read(8*ap_number_bins)
            a = struct.unpack_from('<' + 'd'*ap_number_bins, buf)
            ap_dict[x] = np.array(a)
        fin.close()
        return ap_dict
        
    @classmethod
    def read_overlap(cls, file_path):
        fin = open(file_path, "rb")
        buf = fin.read()
        n = len(buf)//16
        fin.seek(0)
        ov_dict = {'ol_number_bins': n}
        for x in ['ol_range', 'ol_overlap']:
            binary_entry = fin.read(8*n)
            if len(binary_entry) < 8*n:
                raise IOError('Incomplete %s data' % x)
            ov_dict[x] = np.array(struct.unpack_from('<' + 'd'*n, binary_entry))
        fin.close()
        return ov_dict
    
    @classmethod
    def read_deadtime(cls, file_path):
        fin = open(file_path, "rb")
        buf = fin.read()
        n = len(buf)//4
        dt_dict = {'dt_number_coeff': n}
        dt_dict['dt_coeff'] = np.array(struct.unpack_from('<' + 'f'*n, buf))
        dt_dict['dt_coeff_degree'] = np.arange(n - 1, -1, -1)
        return dt_dict
    
    def calculate_snr(self, raw_data, background, background_std_dev):
        background_stacked = np.transpose(np.vstack([background] * raw_data.shape[1]))
        background_std_dev_stacked = np.transpose(np.vstack([background_std_dev] * raw_data.shape[1]))

        return (raw_data - background_stacked)/background_std_dev_stacked
    
    def calculate_dtcf(self, data): # calculate deadtime correction factor
        dt_coeff = self.dt_dict['dt_coeff']
        dt_coeff_degree = self.dt_dict['dt_coeff_degree']
        dt_number_coeff = self.dt_dict['dt_number_coeff']
        return np.sum([(data*1e3)**(dt_coeff_degree[i])*dt_coeff[i] for i in range(dt_number_coeff)], axis=0)
        # # Don't apply the 1000 factor
        #return np.sum([(data)**(dt_coeff_degree[i])*dt_coeff[i] for i in range(dt_number_coeff)], axis=0)

    def calculate_r2_corrected(self, raw_data, background):
        # the result is exactly the same as the software result
        background_stacked = np.transpose(np.vstack([background]*raw_data.shape[1]))
        return (raw_data - background_stacked)*np.square(self.range)

    def calculate_nrb(self, raw_data, background, ap_data, ap_background):
        ap_range = self.ap_dict['ap_range']
        ol_range = self.ov_dict['ol_range']
        ov_data  = self.ov_dict['ol_overlap']
        ap_data_interped = np.interp(self.range, ap_range, ap_data)
        ov_data_interped = np.interp(self.range, ol_range, ov_data)
        
        ap_data_stacked  = np.vstack([ap_data_interped] * raw_data.shape[0])
        ov_data_stacked  = np.vstack([ov_data_interped] * raw_data.shape[0])
        energy_stacked   = np.transpose(np.vstack([self.laser_energy] * raw_data.shape[1]))
        range_stacked    = np.vstack([self.range] * raw_data.shape[0])

        dt_corrected_raw = raw_data * self.calculate_dtcf(raw_data)
        background_correction = np.transpose(np.vstack([background * self.calculate_dtcf(background)] * raw_data.shape[1]))
        afterpulse_correction = (
            ap_data_stacked * self.calculate_dtcf(ap_data_stacked) - ap_background  * self.calculate_dtcf(ap_background)
        ) * energy_stacked  / self.ap_dict['ap_energy']

        # # Not applying 
        # dt_corrected_raw = raw_data
        # background_correction = np.transpose(np.vstack([background] * raw_data.shape[1]))
        # afterpulse_correction = (
        #     ap_data_stacked - ap_background 
        # ) * energy_stacked  / self.ap_dict['ap_energy']

        nrb = (dt_corrected_raw - background_correction - afterpulse_correction) * np.square(range_stacked) / (ov_data_stacked*energy_stacked)
        return nrb

    def calculate_depol_ratio(self):
        cross_over_co = self.nrb_crosspol/self.nrb_copol
        return cross_over_co/(cross_over_co+1)

    def select_time(self, start_time, end_time): # return a new PyMPL object with data within the start_time and end_time
        start_time = np.datetime64(start_time)  # include start
        end_time   = np.datetime64(end_time)    # include end
        selected_indicies = np.where(np.logical_and(self.datetime>=start_time, self.datetime<=end_time))[0]

        selected_data_dict = {}
        for key, value in self.data_dict.items():
            selected_data_dict[key] = [value[i] for i in selected_indicies]
        return PyMPL(selected_data_dict, self.ap_dict, self.ov_dict, self.dt_dict)

    @classmethod
    def select_snr(cls, data, snr_data, snr_limit): # return data with a minimum snr
        high_snr_data = data.copy()
        high_snr_data[snr_data<snr_limit] = np.nan
        return high_snr_data
    
    def interpolate_single_data(self, new_time_array, data, time_reoslution, gap_seconds = None):
        '''
        time_resolution: time resolution in seconds
        gap_seconds: determine when the gap in timestamps is too large to be interpolated
        '''
        data = np.array(data) # make sure its an np array
        if gap_seconds is None:
            gap_seconds = time_reoslution * 2

        new_seconds_passed = (new_time_array-new_time_array[0]).astype('timedelta64[s]').astype(int)
        original_seconds_passed  = (self.datetime-new_time_array[0]).astype('timedelta64[s]').astype(int)

        # this works because it finds the actural datetime of the gap so it does not matter where it started
        time_difference    = self.seconds_since_start[1:] - self.seconds_since_start[:-1] 
        gap_indicies       = np.where(time_difference>gap_seconds)[0]
        gap_starts_arr     = self.datetime[gap_indicies]
        gap_ends_arr       = self.datetime[gap_indicies+1]

        linear_interpolator = scipy.interpolate.interp1d(original_seconds_passed, data, kind='linear', axis= 0, bounds_error = False, fill_value = np.nan)
        interpolated_data   = linear_interpolator(new_seconds_passed) # is a numpy array

        for i in range(len(gap_indicies)):
            gap_starts = gap_starts_arr[i]
            gap_ends   = gap_ends_arr[i]
            new_gap_indicies  = np.where(np.logical_and(new_time_array>gap_starts, new_time_array<gap_ends))[0]
            interpolated_data[new_gap_indicies] = np.nan

        return interpolated_data
    
    @classmethod
    def make_time_array(cls, time_reoslution, start_time, end_time):
        start_time = np.datetime64(start_time)
        end_time   = np.datetime64(end_time) + np.timedelta64(time_reoslution,'s') # so the end time is also included

        return np.arange(start_time, end_time, np.timedelta64(time_reoslution, 's'))
    
    def interpolate_data(self, time_reoslution, start_time = None, end_time = None, gap_seconds = None):
        '''
        time_reoslution is seconds
        '''

        if start_time is None:
            start_time = self.datetime[0]

        if end_time is None:
            end_time = self.datetime[-1] 

        if gap_seconds is None:
            gap_seconds = time_reoslution * 2
        
        # return a new PyMPL object with interpolated data within the start_time and end_time
        
        new_time_array = self.make_time_array(time_reoslution, start_time, end_time)

        self.interpolation_flag                  = True
        self.interpolated_datetime               = new_time_array
        self.interpolated_seconds_since_start    = (new_time_array-new_time_array[0]).astype('timedelta64[s]').astype(int)
        self.interpolated_laser_energy           = self.interpolate_single_data(new_time_array, self.laser_energy, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_temp_detector          = self.interpolate_single_data(new_time_array, self.temp_detector, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_temp_telescope         = self.interpolate_single_data(new_time_array, self.temp_telescope, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_temp_laser             = self.interpolate_single_data(new_time_array, self.temp_laser, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_raw_copol              = self.interpolate_single_data(new_time_array, self.raw_copol, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_raw_crosspol           = self.interpolate_single_data(new_time_array, self.raw_crosspol, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_r2_corrected_copol     = self.interpolate_single_data(new_time_array, self.r2_corrected_copol, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_r2_corrected_crosspol  = self.interpolate_single_data(new_time_array, self.r2_corrected_crosspol, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_nrb_copol              = self.interpolate_single_data(new_time_array, self.nrb_copol, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_nrb_crosspol           = self.interpolate_single_data(new_time_array, self.nrb_crosspol, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_depol_ratio            = self.interpolate_single_data(new_time_array, self.depol_ratio, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_snr_copol              = self.interpolate_single_data(new_time_array, self.snr_copol, time_reoslution, gap_seconds=gap_seconds)
        self.interpolated_snr_crosspol           = self.interpolate_single_data(new_time_array, self.snr_crosspol, time_reoslution, gap_seconds=gap_seconds)

    def interpolation_reset(self):
        self.interpolation_flag                  = False
        self.interpolated_datetime               = None
        self.interpolated_seconds_since_start    = None
        self.interpolated_laser_energy           = None
        self.interpolated_temp_detector          = None
        self.interpolated_temp_telescope         = None
        self.interpolated_temp_laser             = None
        self.interpolated_raw_copol              = None
        self.interpolated_raw_crosspol           = None
        self.interpolated_r2_corrected_copol     = None
        self.interpolated_r2_corrected_crosspol  = None
        self.interpolated_nrb_copol              = None
        self.interpolated_nrb_crosspol           = None
        self.interpolated_depol_ratio            = None
        self.interpolated_snr_copol              = None
        self.interpolated_snr_crosspol           = None
    
    def write_mpl(self, output_dir, filename):
        """
        output the current data dict as a new binary mpl file.
        Useful when need to need to modify the .mpl binary file by changing self.data_dict
        """
        output_path = os.path.join(output_dir, filename+'.mpl')
        new_dict = self.data_dict.copy()
        
        new_file = open(output_path, 'wb')
        for i in range(self.number_profile):
            for j in range(len(self._record_entry)):
                entry_type = self._record_entry[j]
                entry_format = self._record_format[j]
                data_to_write = new_dict[entry_type][i]
                if entry_format in ['uint16','uint32','int16','int32','float32']:
                    binary_data = struct.pack(self._byte_format[entry_format], data_to_write)
                elif entry_format == 'byte':
                    integer_value = int(data_to_write, 2)
                    binary_data = struct.pack(self._byte_format[entry_format], integer_value)
                elif entry_format == 'char*6':
                    binary_data = data_to_write.encode('latin-1')
                elif entry_format == 'float32 array':
                    binary_data = struct.pack(self._byte_format[entry_format], *data_to_write)
                    
                new_file.write(binary_data)

    @classmethod
    def _movingaverage(cls, values, window):
        # helper function moving average with 1d array
        if window%2 == 0:
            raise ValueError('Use odd moving average window')
        weights = np.repeat(1.0, window)/window
        sma = np.convolve(values, weights, mode='same')
        return sma

    @classmethod
    def movingaverage(cls, values, window):
        # values could either be an 1d array or a 2d array.
        values = np.array(values)
        if values.ndim == 1:
            return cls._movingaverage(values, window)
        elif values.ndim == 2:
            result = []
            for i in range(values.shape[0]):
                result.append(cls._movingaverage(values[i], window))
            result = np.array(result)
            return result
        else:
            raise ValueError('values has to be an 1D array or a 2D array')
        

