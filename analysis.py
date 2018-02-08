# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:02:06 2018

@author: spauka
"""

import qcodes as qc
import numpy as np
from qcodes.plots.colors import color_cycle, colorscales
import os, re
from time import sleep

def find_data(date, num):
    date = str(date)
    data_dir = os.path.join("data", date)
    files = os.listdir(data_dir)
    for file in files:
        if re.match("#{:03}".format(num), file):
            return os.path.join(data_dir, file)
        
def open_data(date, num):
    loc = find_data(date, num)
    if loc is None:
        return None
    data = qc.DataSet(location=loc)
    data.read()
    data.read_metadata()
    return data
        
def format_plot(plot, left_axis):
    # First, let's resize the window
    plot.win.resize(600, 450)
    
    # Then, get the plotitem
    pl = plot.win.getItem(0, 0)
    
    # Set the Axes
    pl.setLabel("left", left_axis[0], left_axis[1])
    pl.setLabel("bottom", "Perpendicular Field", "T")

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

def plot_all_field_sweeps(date, start_num):
    """
    Plot field sweeps where data is split into multiple consecutive datasets.
    This is the case when one DataSet is taken for each switch configuration
    """
    plots = []
    for i, sw in enumerate(switches):
        data = qc.DataSet(location=find_data(date, start_num+i))
        data.read()
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
            plot.save(filename="{}_{}.png".format(sw, name))
    return plots

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
