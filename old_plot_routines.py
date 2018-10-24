import qcodes as qc
from qcodes.plots.colors import color_cycle, colorscales
import numpy as np
from time import sleep

from data_utils import open_data

def format_plot(plot, left_axis):
    # First, let's resize the window
    plot.win.resize(600, 450)
    
    # Then, get the plotitem
    pl = plot.win.getItem(0, 0)
    
    # Set the Axes
    pl.setLabel("left", left_axis[0], left_axis[1])
    pl.setLabel("bottom", "Perpendicular Field", "T")
    
def plot_cooldown_update(date, num):
    plots = []
    data = open_data(date, num)
    
    for lockin in range(1, 2):
        for i, switch in enumerate(("A", "B", "C", "D", "E")):
            plot = qc.QtPlot()
            name = "{}_SR860_{}_X".format(switch, lockin)
            name_res = "{}_SR860_{}_resistance".format(switch, lockin)
            label = "V<sub>xy</sub>" if (lockin%2) == 1 else "V<sub>xx</sub>"
            label = (label, "V")
            tr = getattr(data, name_res, [])
            if not tr:
                tr = getattr(data, name, [])
            else:
                name = name_res
                label = ("Resistance", "Ohms")
            plot.add(tr, name=name, color=color_cycle[((lockin-1)*5 + i) % 10])
            plot.win.resize(600, 450)
            plots.append(plot)
        
    while True:
        data.read()
        for plot in plots:
            plot.update()
        sleep(10)

def plot_cooldown(datas):
    """
    Plot resistance of ohmics during cooldown
    Datas is the datasets that make up the cooldown
    """
    
    # Calculate how many data points are in each array
    lengths = []
    for data in datas:
        time = data.time
        lengths.append(np.where(np.isnan(time))[0][0] - 1)
    
    # Make a new DataSet
    new_data = qc.DataSet(location="data/Cooldown")
    new_length = np.sum(lengths)
    # And add new dataarrays for each dataarray in the original data
    for d_arr in datas[0].arrays.values():
        data_array = qc.DataArray(
                name=d_arr.name,
                full_name=d_arr.full_name, 
                label=d_arr.full_name,
                array_id=d_arr.array_id,
                unit=d_arr.unit,
                is_setpoint=d_arr.is_setpoint,
                shape=(new_length, *d_arr.shape[1:]))
        data_array.init_data()
        new_data.add_array(data_array)
    
    # Then, update each of the set arrays
    for key, d_arr in new_data.arrays.items():
        d_arr.set_arrays = tuple(new_data.arrays[s_arr.name + "_set"] for s_arr in datas[0].arrays[key].set_arrays)
        
    # Then, fill in each item
    cumsum = 0
    for data, length in zip(datas, lengths):
        for key, d_arr in data.arrays.items():
            new_data.arrays[key][cumsum:cumsum+length] = d_arr.ndarray[0:length]
        cumsum += length
    
    # We also need to make time keep counting upwards
    cumsum = 0
    offs = 0
    for i, l in enumerate(lengths[:-1]):
        cumsum += l
        offs += new_data.time[l-1]
        new_data.time[cumsum:cumsum+lengths[i+1]] += offs
    
    return new_data

def plot_all_field_sweeps_comb(data):
    """ 
    Plot field sweeps under the condition that all sweeps are included
    in a single data file
    """
    plots = []
    for sw in ("B", "C", "D", "E"):
        for lockin in range(1, 4):
            plot = qc.QtPlot()
            name = "{}_SR860_{}_X".format(sw, lockin)                          
            name_res = "{}_SR860_{}_resistance".format(sw, lockin)
            label = "V<sub>xy</sub>" if (lockin%2) == 1 else "V<sub>xx</sub>"
            label = (label, "V")
            tr = getattr(data, name_res, [])
            if not tr:
                tr = getattr(data, name, [])
            else:
                name = name_res
                label = ("Resistance", "Ohms")
            plot.add(tr, name=name, color=color_cycle[lockin-1])
            format_plot(plot, label)
            plots.append(plot)
            plot.save(filename="{}.png".format(name))
    return plots

def plot_all_field_sweeps(date, start_num, switches=("B", "C", "D", "E")):
    """
    Plot field sweeps where data is split into multiple consecutive datasets.
    This is the case when one DataSet is taken for each switch configuration
    """
    plots = []
    for i, sw in enumerate(switches):
        data = open_data(date, start_num+i)
        for lockin in range(1, 5):
            plot = qc.QtPlot()
            name = "SR860_{}_X_preamp".format(lockin)
            name_res = "Voltage_{}".format(lockin)                                 # fyi I changed 
            label = "V<sub>xx</sub>" if (lockin%2) == 1 else "V<sub>xy</sub>"   # some of these
            label = (label, "V")
            tr = getattr(data, name_res, [])
            if not tr:
                tr = getattr(data, name, [])
            else:
                name = name_res
                label = ("Resistance", "Ohms")
            plot.add(tr, name=name, color=color_cycle[lockin-1])
            format_plot(plot, label)
            plots.append(plot)
            plot.save(filename="{}_{}.png".format(sw, name))
    return plots

def plot_update_currents(date, num):
    plots = []
    data = open_data(date, num)
    
    for lockin in range(1, 4):
        plot = qc.QtPlot()
        name_ithaco = "SR860_{}_X_ithaco".format(lockin)
        name_volt = "SR860_{}_X_preamp".format(lockin)
        tr = getattr(data, name_ithaco, [])
        if not tr:
            tr = getattr(data, name_volt, [])
            label = ("Voltage", "V")
            name = name_volt
        else:
            label = ("Current", "A")
            name = name_ithaco
        plot.add(tr, name=name, color=color_cycle[lockin-1])
        format_plot(plot, label)
        plots.append(plot)
        
    while True:
        data.read()
        for plot in plots:
            plot.update()
        sleep(10)

def plot_update(date, num):
    plots = []
    data = open_data(date, num)
    
    for lockin in range(1, 5):
        plot = qc.QtPlot()
        name = "SR860_{}_X".format(lockin)
        name_res = "SR860_{}_resistance".format(lockin)
        label = "V<sub>xy</sub>" if (lockin%2) == 1 else "V<sub>xx</sub>"
        label = (label, "V")
        tr = getattr(data, name_res, [])
        if not tr:
            tr = getattr(data, name, [])
        else:
            name = name_res
            label = ("Resistance", "Ohms")
        plot.add(tr, name=name, color=color_cycle[lockin-1])
        format_plot(plot, label)
        plots.append(plot)
        
    while True:
        data.read()
        for plot in plots:
            plot.update()
        sleep(10)
        
def add_label(plot, text, posn=(100, 30)):
    qtplot = plot.win._handler._import("pyqtgraph")
    l = qtplot.LegendItem(offset=posn)
    l.setParentItem(plot.win.getItem(0, 0))
    txt = qtplot.LabelItem()
    txt.setAttr("color", "000000")
    txt.setText(text)
    l.layout.addItem(txt, 0, 0)
 