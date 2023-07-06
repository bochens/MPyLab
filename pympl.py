import glob
import numpy as np
import os
import struct
from netCDF4 import Dataset

class PyMPL:
    _C = 299792458 # m.s-1 speed of light in the universe where this code is written

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

    def __init__(self, data_input, ap_input, ov_input, dt_input):
        
        # make sure the instance constructor can take both dictionary and file path as input
        if all([isinstance(element, dict) for element in [data_input, ap_input, ov_input, dt_input]]):
            self.data_dict, self.ap_dict, self.ov_dict, self.dt_dict = data_input, ap_input, ov_input, dt_input
        else:
            self.data_dict, self.ap_dict, self.ov_dict, self.dt_dict = self.read_files(data_input, ap_input, ov_input, dt_input)

        year_str    = np.char.zfill(np.array(self.data_dict['year']).astype(str), 4)
        month_str   = np.char.zfill(np.array(self.data_dict['month']).astype(str), 2)
        day_str     = np.char.zfill(np.array(self.data_dict['day']).astype(str), 2)
        hours_str   = np.char.zfill(np.array(self.data_dict['hours']).astype(str), 2)
        minutes_str = np.char.zfill(np.array(self.data_dict['minutes']).astype(str), 2)
        seconds_str = np.char.zfill(np.array(self.data_dict['seconds']).astype(str), 2)
        self.datetime = np.array([x+'-'+y+'-'+z+'T'+h+':'+m+':'+s \
                                  for x,y,z,h,m,s in zip(year_str,month_str,day_str,hours_str,minutes_str,seconds_str)], dtype='datetime64')

        self.energy_monitor    = np.array(self.data_dict['energy_monitor'])/1000
        self.temp_detector = np.array(self.data_dict['temp_0'])/100
        self.temp_telescope = np.array(self.data_dict['temp_2'])/100
        self.temp_laser = np.array(self.data_dict['temp_3'])/100

        self.number_bins = self.data_dict['number_bins'][0]
        if np.all(np.array(self.data_dict['number_bins']) != self.number_bins):
            raise ValueError('number_bins of scans have to be the same')
        
        self.bin_time = self.data_dict['bin_time'][0]
        if np.all(np.array(self.data_dict['bin_time']) != self.bin_time):
            raise ValueError('bin_time of scans have to be the same')

        self.raw_data_copol = np.array(self.data_dict['channel_1_data'])
        self.raw_data_crosspol = np.array(self.data_dict['channel_2_data'])
        background_copol = np.array(self.data_dict['background_average'])
        background_crosspol = np.array(self.data_dict['background_average_2'])
        self.range = 0.5*self.bin_time*self._C*(np.arange(self.number_bins) + 0.5)*1e-3 #km

        # Calculate Range Correction Product
        self.r2_corrected_copol    = self.calculate_r2_corrected(self.raw_data_copol, background_copol)
        self.r2_corrected_crosspol = self.calculate_r2_corrected(self.raw_data_crosspol, background_crosspol)

        # Calculate Normalized Relative Backscatter
        self.nrb_copol    = self.calculate_nrb(self.raw_data_copol   ,background_copol    \
                                               ,self.ap_dict['ap_copol']   ,self.ap_dict['ap_background_average_copol'])
        self.nrb_crosspol = self.calculate_nrb(self.raw_data_crosspol,background_crosspol \
                                               ,self.ap_dict['ap_crosspol'],self.ap_dict['ap_background_average_crosspol'])
    
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
        ol_dict = {'ol_number_bins': n}
        for x in ['ol_range', 'ol_overlap']:
            binary_entry = fin.read(8*n)
            if len(binary_entry) < 8*n:
                raise IOError('Incomplete %s data' % x)
            ol_dict[x] = np.array(struct.unpack_from('<' + 'd'*n, buf))
        fin.close()
        return ol_dict
    
    @classmethod
    def read_deadtime(cls, file_path):
        fin = open(file_path, "rb")
        buf = fin.read()
        n = len(buf)//4
        dt_dict = {'dt_number_coeff': n}
        dt_dict['dt_coeff'] = np.array(struct.unpack_from('<' + 'f'*n, buf))
        dt_dict['dt_coeff_degree'] = np.arange(n - 1, -1, -1)
        return dt_dict
    
    def calculate_dtcf(self, raw_data): # calculate deadtime correction factor
        dt_coeff = self.dt_dict['dt_coeff']
        dt_coeff_degree = self.dt_dict['dt_coeff_degree']
        dt_number_coeff = self.dt_dict['dt_number_coeff']
        return np.sum([(raw_data*1e3)**(dt_coeff_degree[i])*dt_coeff[i] for i in range(dt_number_coeff)], axis=0)

    def calculate_r2_corrected(self, raw_data, background):
        background_stacked = np.transpose(np.vstack([background]*raw_data.shape[1]))
        return (raw_data - background_stacked)*np.square(self.range)
    
    def calculate_nrb(self, raw_data, background, ap_data, ap_background):
        ap_range = self.ap_dict['ap_range']
        ol_range = self.ov_dict['ol_range']
        ov_data  = self.ov_dict['ol_overlap']
        ap_data_interped = np.interp(self.range, ap_range, ap_data) if self.ap_dict['ap_number_bins']!=self.number_bins else ap_data
        ov_data_interped = np.interp(self.range, ol_range, ov_data) if self.ov_dict['ol_number_bins']!=self.number_bins else ov_data
        
        ap_data_stacked  = np.vstack([ap_data_interped]*raw_data.shape[0])
        ov_data_stacked  = np.vstack([ov_data_interped]*raw_data.shape[0])
        energy_stacked   = np.transpose(np.vstack([self.energy_monitor]*raw_data.shape[1]))
        range_stacked    = np.vstack([self.range]*raw_data.shape[0])
        dt_corrected_raw = raw_data*self.calculate_dtcf(raw_data)
        background_correction = np.transpose(np.vstack([background*self.calculate_dtcf(background)]*raw_data.shape[1]))
        afterpulse_correction = (ap_data_stacked*self.calculate_dtcf(ap_data_stacked) \
                                 - ap_background*self.calculate_dtcf(ap_background))*energy_stacked/self.ap_dict['ap_energy']
        nrb = (dt_corrected_raw - background_correction - afterpulse_correction)*np.square(range_stacked)/(ov_data_stacked*energy_stacked)
        return nrb
    
    def get_closest_time_index(self, a_datetime):
        a_datetime = np.datetime64(a_datetime) # makesure the datetime type is numpy datetime65
        time_difference = (self.datetime - a_datetime)
        time_difference_seconds = time_difference.astype('timedelta64[s]').astype(int)
        nearest_index = np.argmin(np.abs(time_difference_seconds))
        return nearest_index

