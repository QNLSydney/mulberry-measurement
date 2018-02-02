# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:02:06 2018

@author: spauka
"""

import qcodes as qc
from qcodes.plots.colors import color_cycle, colorscales
import os, re

def find_data(date, num):
    date = str(date)
    data_dir = os.path.join("data", date)
    files = os.listdir(data_dir)
    for file in files:
        if re.match("#{:03}".format(num), file):
            return os.path.join(data_dir, file)
        
def format_plot(plot, left_axis):
    # First, let's resize the window
    plot.win.resize(600, 450)
    
    # Then, get the plotitem
    pl = plot.win.getItem(0, 0)
    
    # Set the Axes
    pl.setLabel("left", left_axis[0], left_axis[1])
    pl.setLabel("bottom", "Perpendicular Field", "T")
    
def plot_all(data):
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

def plot_all_2(date, start_num):
    plots = []
    for i, sw in enumerate(("B", "C", "D", "E")):
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