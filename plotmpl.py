import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import ListedColormap
from pympl import PyMPL

### lidar jet colormap ###
jet = cm.get_cmap('jet', 256)
newcolors = jet(np.linspace(0, 1, 256))
newcolors[:1, :] = np.array([0, 0, 0, 1])  # black
newcolors[-1:, :] = np.array([1, 1, 1, 1]) # white
lidar_cmap = ListedColormap(newcolors)

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

def plot_mpl_2d_timeseries(mpl_datetime, mpl_range_edges, mpl_2d_data, fig=None, ax=None, range_max = 20, vmin=0, vmax=1, colorbar_bool = True, tick_number = None):
    '''
    range_max: max y range to plot in km
    '''
    fig, ax = _get_fig_and_ax(fig, ax)
    ax.set_facecolor('k')
    number_of_timestamp = mpl_datetime.shape[0]   # number of mpl profiles
    time_x = np.arange(number_of_timestamp+1)     # time dimension edge for pcolormesh

    handle = ax.pcolormesh(time_x, mpl_range_edges, np.transpose(mpl_2d_data), cmap = lidar_cmap, shading='auto' , norm=colors.Normalize(vmin=0, vmax=1))
    ax.set_ylim(0,range_max)

    if colorbar_bool:
        fig.colorbar(handle) # draw colorbar

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height

    # x ticking and label
    if tick_number is None:
        tick_number = int(width) # will automatically determine how many ticks and lebels to use depending on the figure width
    x_tick_positions_availiable = time_x[:-1] + 0.5
    x_tick_position_selector = np.linspace(x_tick_positions_availiable[0], x_tick_positions_availiable[-1], num = tick_number).astype(int)
    ax.set_xticks(x_tick_positions_availiable[x_tick_position_selector])
    formatted_time_label_string = _format_datetime(mpl_datetime[x_tick_position_selector])
    ax.set_xticklabels(formatted_time_label_string)

    ax.set_xlabel('time')
    ax.set_ylabel('range (km)')
    