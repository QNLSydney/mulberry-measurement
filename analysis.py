# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:02:06 2018

@author: spauka
"""

import qcodes as qc
import numpy as np
import scipy.constants as const
import pyqtgraph as qtplot
from qcodes.plots.colors import color_cycle, colorscales
import os, re
from time import sleep
from pathlib import Path

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
    
def plot_cooldown_update(date, num):
    plots = []
    data = open_data(date, num)
    
    for lockin in range(1, 4):
        plot = qc.QtPlot()
        for i, switch in enumerate(("A", "B", "C", "D", "E")):
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
            if switch == "C" and lockin == 3:
                continue
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

def plot_update_currents(date, num):
    plots = []
    data = open_data(date, num)
    
    for lockin in range(1, 4):
        plot = qc.QtPlot()
        name_ithaco = "SR860_{}_X_ithaco".format(lockin)
        name_curr = "SR860_{}_X".format(lockin)
        label = "Current"
        label = (label, "A")
        tr = getattr(data, name_ithaco, [])
        if not tr:
            tr = getattr(data, name_curr, [])
        else:
            name = name_curr
        plot.add(tr, name=name_ithaco, color=color_cycle[lockin-1])
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
    
def calc_rho(field, data_xx, data_xy, curr=1e-9, width=50e-6, length=100e-6, 
             name="Analyzed_Field_Sweep"):
    path = "data" / Path(name)
    
    np_field = field.ndarray.copy()
    np_xx = data_xx.ndarray.copy()
    np_xy = data_xy.ndarray.copy()
    
    # Calculate resistances
    rho_xy = np_xy/curr
    # rho_xx is scaled by the width+length of the hall bar to get a sheet resistivity
    rho_xx = (np_xx/curr) * (width/length)
    
    # Calculate density
    # We want to do a fit between 1, -1 tesla as there is some extra structure
    # above these values
    min_ind = np.where(np.abs(np_field + 1) < 0.005)[0][0]
    max_ind = np.where(np.abs(np_field - 1) < 0.005)[0][0]
    min_ind, max_ind = np.min((min_ind, max_ind)), np.max((min_ind, max_ind))
    res = np.polyfit(np_field[min_ind:max_ind], rho_xy[min_ind:max_ind], 1)
    poly = np.poly1d(res)
    
    # Then the density is given by 1/(|e| dp/dB), in cm^2
    density = 1/(const.e * res[0])
    density *= 1e-4
    
    print("Density is: {:.2e} cm^-2".format(density))
    
    # And the calculated mobility is given by mu=1/(rho_xx |e| ns)
    mu = 1/(rho_xx * const.e * density)
    
    # And let's quote the density slightly offset from 0, at 0.1T
    mob_ind = np.where(np.abs(np_field - 0.1) < 0.005)[0][0]
    mobility = mu[mob_ind]
    
    print("Mobility is: {:.2e} cm^2/V s".format(mobility))
    
    # And finally, let's create a new dataset. Leave location unset for now...
    dataset = qc.DataSet(location=str(path))
    
    da_field = qc.DataArray(
            array_id="Field", label="Field", unit="T", is_setpoint=True,
            preset_data=np_field)
    da_reduced_field = qc.DataArray(
            array_id="Reduced_Field", label="Field", unit="T", is_setpoint=True,
            preset_data=np_field[min_ind:max_ind])
    da_poly = qc.DataArray(
            array_id="Poly_Deg", label="Polynomial Degree", is_setpoint=True,
            preset_data=list(range(2))[::-1])
    
    da_rho_xy = qc.DataArray(
            array_id="Rho_XY", label="Rho XY", unit="Ohms", 
            set_arrays=(da_field,), preset_data=rho_xy)
    da_rho_xy_fit = qc.DataArray(
            array_id="fit_Rho_XY", label="Rho XY", unit="Ohms",
            set_arrays=(da_reduced_field,), preset_data=poly(np_field[min_ind:max_ind]))
    da_rho_xy_coef = qc.DataArray(
            array_id="coef_Rho_XY", label="Poly Coefficients",
            set_arrays=(da_poly,), preset_data=res)
    da_rho_xx = qc.DataArray(
            array_id="Rho_XX", label="Rho XX", unit="Ohms",
            set_arrays=(da_field, ), preset_data=rho_xx)
    da_mu = qc.DataArray(
            array_id="mu", label="Mobility", unit="cm<sup>2</sup> (V s)<sup>-1</sup>",
            set_arrays=(da_field, ), preset_data=mu)
    
    dataset.add_array(da_field)
    dataset.add_array(da_reduced_field)
    dataset.add_array(da_poly)
    dataset.add_array(da_rho_xy)
    dataset.add_array(da_rho_xy_fit)
    dataset.add_array(da_rho_xy_coef)
    dataset.add_array(da_rho_xx)
    dataset.add_array(da_mu)
    
    # Save the data
    dataset.finalize()
    
    # Make some nice plots
    # Showing Density (Rho_XY) analysis
    plot_rho_xy = qc.QtPlot()
    plot_rho_xy.add(dataset.Rho_XY)
    plot_rho_xy.add(dataset.fit_Rho_XY)
    add_label(plot_rho_xy, "Using {} current<br>".format(qtplot.siFormat(curr, suffix="A")) +
                      "From a linear fit:<br>" +
                      "dρ/dB = {}<br>".format(qtplot.siFormat(res[0], suffix="Ω")) +
                      "n<sub>s</sub> = 1/(|e| dρ/dB) = {:e} cm<sup>-2</sup>".format(density),
            posn=(100, 30))
    plot_rho_xy.save(filename=str(path/"rho_xy.png"))
    # Showing Mobility analysis
    plot_mob = qc.QtPlot()
    plot_mob.add(dataset.mu, color=color_cycle[5])
    add_label(plot_mob, "Mobility extracted from:<br>" +
                        "μ = 1/ρ<sub>xx</sub> |e| n<sub>s</sub>, with n<sub>s</sub>= {:.2e} cm<sup>-2</sup><br>".format(density) +
                        "And using W = {}, L = {}".format(qtplot.siFormat(width, suffix="m"),
                                       qtplot.siFormat(length, suffix="m")),
                                       posn=(-30, -60))
    plot_mob.save(filename=str(path/"mobility.png"))
    
    return dataset
