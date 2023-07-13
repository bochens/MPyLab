import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import ListedColormap


def _make_lidar_colormap(cmap_str, under_color=np.array([0, 0, 0, 1]), over_color=np.array([1, 1, 1, 1])):
    cmap = cm.get_cmap(cmap_str, 256)
    newcolors = cmap(np.linspace(0, 1, 256))
    newcolors[:1, :] = under_color  # black
    newcolors[-1:,:] = over_color # white
    lidar_cmap = ListedColormap(newcolors)
    return lidar_cmap

##### lidar colormap #####
gnuplot2     = cm.get_cmap('gnuplot2', 256)
nipy_spectral= cm.get_cmap('nipy_spectral', 256)
lidar_jet    = _make_lidar_colormap('jet')
lidar_gist   = _make_lidar_colormap('gist_ncar')
lidar_turbo  = _make_lidar_colormap('turbo')

##### default parameters #####
default_figure_size = (7.2, 4.8)

def _get_fig_and_ax(fig, ax):
    if fig is None and ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=default_figure_size)
    elif fig is None and ax is not None:
        fig = ax.get_figure()
    elif fig is not None and ax is None:
        ax = fig.add_subplot(111)
    return fig, ax

def _format_datetime(datetime_arr):
    return_list = []
    formatted_str = ""

    date_and_time_format = "%H:%M:%S\n%m-%d-%Y"
    time_format = '%H:%M:%S'
    
    for idx, x in enumerate(datetime_arr):
        if idx == 0:
            formatted_str += x.item().strftime(date_and_time_format)
        else:
            previous_x = datetime_arr[idx-1]
            if previous_x.astype('datetime64[D]') == x.astype('datetime64[D]'):
                formatted_str += x.item().strftime(time_format)
            else:
                formatted_str += x.item().strftime(date_and_time_format)
        return_list.append(formatted_str)
        formatted_str = ""
    return np.array(return_list)

def plot_mpl_2d_timeseries(mpl_datetime, mpl_range, mpl_2d_data, fig=None, ax=None, range_max = 20, vmin=0, vmax=1, color_map = lidar_jet, colorbar_bool = True, colorbar_norm = 'linear', x_tick_number = None):
    '''
    range_max: max y range to plot in km
    '''
    fig, ax = _get_fig_and_ax(fig, ax)
    ax.set_facecolor('k')
    time_x  = np.arange(mpl_datetime.shape[0]+1)     # time dimension edge for pcolormesh
    y_resolution = mpl_range[1]-mpl_range[0]
    range_y = np.append(mpl_range-y_resolution/2, mpl_range[-1]+y_resolution/2)

    if vmin <= 0:
        log_vmin = np.nanmin(mpl_2d_data[mpl_2d_data>0])
    else:
        log_vmin = vmin

    if colorbar_norm == "linear":
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
    elif (colorbar_norm == "log") or (colorbar_norm == "Log"):
        norm = colors.LogNorm(vmin=log_vmin, vmax=vmax)

    handle = ax.pcolormesh(time_x, range_y, np.transpose(mpl_2d_data), cmap = color_map, shading='auto' , norm=norm)
    ax.set_ylim(0,range_max)

    if colorbar_bool:
        fig.colorbar(handle) # draw colorbar

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height

    # x ticking and label
    if x_tick_number is None:
        x_tick_number = int(width) # will automatically determine how many ticks and lebels to use depending on the figure width

    x_tick_positions_availiable = time_x[:-1] + 0.5
    x_tick_position_selector = np.linspace(x_tick_positions_availiable[0], x_tick_positions_availiable[-1], num = x_tick_number).astype(int)
    ax.set_xticks(x_tick_positions_availiable[x_tick_position_selector])
    formatted_time_label_string = _format_datetime(mpl_datetime[x_tick_position_selector])
    ax.set_xticklabels(formatted_time_label_string)

    ax.set_xlabel('time')
    ax.set_ylabel('range (km)')

def plot_mpl_1d_timeseries(mpl_datetime, mpl_1d_data, fig=None, ax=None, ylabel = None, tick_number = None, color = '#1f77b4'):
    fig, ax = _get_fig_and_ax(fig, ax)

    time_x = np.arange(mpl_datetime.size)
    handle = ax.plot(time_x, mpl_1d_data, color = color)

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height

    # x ticking and label
    if tick_number is None:
        tick_number = int(width) # will automatically determine how many ticks and lebels to use depending on the figure width

    x_tick_positions_availiable = time_x
    x_tick_position_selector = np.linspace(time_x[0], time_x[-1], num = tick_number).astype(int)
    ax.set_xticks(x_tick_positions_availiable[x_tick_position_selector])
    formatted_time_label_string = _format_datetime(mpl_datetime[x_tick_position_selector])
    ax.set_xticklabels(formatted_time_label_string)

    ax.set_xlabel('time')

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    
#def profile_plot