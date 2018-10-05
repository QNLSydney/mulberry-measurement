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
from scipy.optimize import curve_fit
from scipy.special import digamma
from scipy.integrate import quad
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.fftpack
import datetime

def find_data(date, num):
    date = str(date)
    data_dir = os.path.join("data", date)
    files = os.listdir(data_dir)
    for file in files:
        if re.match("#{:03}".format(num), file):
            return os.path.join(data_dir, file)
    return None
        
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
        data = qc.DataSet(location=find_data(date, start_num+i))
        data.read()
        for lockin in range(1, 5):
            plot = qc.QtPlot()
            name = "SR860_{}_X_preamp".format(lockin)
            name_res = "Voltage".format(lockin)                                 # fyi I changed 
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
    
def calc_rho(field, data_xx, data_xy, curr=1e-9, width=50e-6, length=100e-6, 
             field_center=1, name="Analyzed_Field_Sweep"):
    path = "data" / Path(name)
    
    np_field = field.ndarray.copy()
    np_xx = data_xx.ndarray.copy()
    np_xy = data_xy.ndarray.copy()
    
    curr = np.average(-1*(data.SR860_3_X_preamp.ndarray.copy())/1e6)
    print(curr)
    
    # Calculate resistances
    rho_xy = np_xy/curr
    # rho_xx is scaled by the width+length of the hall bar to get a sheet resistivity
    rho_xx = (np_xx/curr) * (width/length)
    
    # Calculate density
    # We want to do a fit between (field_center,-field_center) tesla as there is some extra structure
    # above these values
    min_ind = np.where(np.abs(np_field + field_center) < 0.01)[0][0]               
    max_ind = np.where(np.abs(np_field - field_center) < 0.01)[0][0]    
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




def leak_test_plot(date, num, width=50e-6):
    """
    
    """
    data = open_data(date, num)   
                        
    V_xx = (data.SR860_1_X_preamp) 
    #V_xy = (data.SR860_2_X_preamp)
    curr = (data.SR860_3_X_preamp)
    #Q = np.arctan(data.SR860_3_Y_preamp/data.SR860_3_X_preamp)*180/3.141
    V_g = (data.yoko_t_voltage_set)
      
    # Calculate resistances
    #rho_xy = V_xy/curr
    # rho_xx is scaled by the width+length of the hall bar to get a sheet resistivity
    #rho_xx = (V_xx/curr) * (width/100e-6)
    curr = curr

    fig, ax1 = plt.subplots()
    ax1.plot(V_g,V_xx,'rx', label = '{} um device'.format(1e6*width))
    ax1.set_xlabel('top gate voltage (V)')
    ax1.set_ylabel('V_xx (V)', color='r')
    ax1.tick_params('y', colors='r')
    ax1.set_ylim([0.0,2e-5])
    
    ax2 = ax1.twinx()
    ax2.plot(V_g,np.abs(curr),'b-')
    ax2.set_xlabel('top gate voltage (V)')
    ax2.set_ylabel('current (uA)', color='b')
    ax2.tick_params('y', colors='b')  
    ax2.set_ylim([0, 0.0042])
    
    fig.tight_layout()
    plt.legend()
    plt.show()
    name = "leak_test_{}_{}_".format(date,num)
    plt.savefig(name)
      
    
    

def Hikami_fit(date, num, current=False, width=50e-6, length=100e-6, direction=1, l_mfp=3.51e-7):
    """
    Fits weak (anti-)localization with the Hikami-Larkin-Nagaoka formula 
    (or the modification in the strong spin orbit coupling regime).
    Use with "top_gate_step_B_field_sweep" function.
    """
    # Define the fitting functions
    e = 1.60217662e-19
    h_bar = 1.0545718e-34
    pi = 3.141
    c = e**2/(pi**2*h_bar)
    alpha = +0.5     # is -1 for weak localization and +1/2 for weak anti-localization
    
    B_e_mfp = h_bar/(4.*e*l_mfp**2)
        
    def func(B, B_phi, B_so, B_e=B_e_mfp):                   
        return (1/2)*c*(np.log(B_phi/B)-digamma((1/2)+(B_phi/B)))+c*(np.log((B_so+B_e)/B)-digamma((1/2)+((B_so+B_e)/B)))-(3/2)*c*(np.log(((4/3)*B_so+B_phi)/B)-digamma((1/2)+((4/3)*B_so+B_phi)/B))
    
    #(1/2)*c*(digamma((1/2)+(B_phi/B)+(B_so/B)) + (1/2)*digamma((1/2)+(B_phi/B)+(2*B_so/B)) - (1/2)*digamma((1/2)+(B_phi/B)) - np.log((B_so+B_phi)/B) - (1/2)*np.log((B_so+B_phi)/B) + (1/2)*np.log((B_phi)/B))
    #(1/2)*c*(np.log(B_phi/B)-digamma((1/2)+(B_phi/B)))+c*(np.log((B_so+B_e)/B)-digamma((1/2)+((B_so+B_e)/B)))-(3/2)*c*(np.log(((4/3)*B_so+B_phi)/B)-digamma((1/2)+((4/3)*B_so+B_phi)/B))
    
    def func_strong_so(B, B_phi):
        return alpha*(1/2)*c*(np.log(B_phi/B)-digamma((1/2)+(B_phi/B)))
    
    # Import the data
    data = open_data(date, num)
    if data == None:
        return 'change'
    xdata = (data.ami_field_set.ndarray.copy()*direction) 
    ydata = data.SR860_1_X_preamp.ndarray.copy()     
    if current == True:
        curr = data.SR860_3_X_preamp.ndarray.copy()*1e-6
    else:
        curr = 1e-9
    
    # Calculate conductivity
    # rho_xx is scaled by the width and length of the hall bar to get a sheet resistivity
    rho_xx = (ydata/curr) * (width/length)
    sig_xx = 1/rho_xx
    
    # Find Zero values (average)
    a1 = np.where(np.abs(xdata) < 0.00005)    # select low field ange
    print(a1)
    b1 = np.argmax(sig_xx[a1])               # find sigma max in this rage
    c1 = a1[0]
    sig_xx_0 = sig_xx[c1[0]+b1]
    xdata = xdata - xdata[c1[0]+b1]          # corrects slight field offset

    # Calculate Magnetoconductivity
    D_sig_xx = sig_xx - sig_xx_0
    
    #plt.plot(xdata, D_sig_xx, 'b-', label='data')
    #plt.xlim(-0.005,0.005)
    plt.xlabel('B (T)')
    plt.ylabel('Delta conductance (S)')
    #plt.ylim(-0.00005,0.000003)
        
    # Specify range of data to be fit (> 0)
    fit_data_x = xdata[np.where(xdata > 0.0)]
    fit_data_y = D_sig_xx[np.where(xdata > 0.0)]
    
    B_max = B_e_mfp
    
    fit_data_y = fit_data_y[np.where(fit_data_x < B_max)]
    fit_data_x = fit_data_x[np.where(fit_data_x < B_max)]    
       
    trial_values = np.array([1.1e-8,1.8e-8])#,1.5e-3])
         
    popt, pcov = curve_fit(func, fit_data_x, fit_data_y, p0=trial_values, bounds=([0.0,0.0],[1.1,1.1]))
    plt.plot(fit_data_x, fit_data_y, 'b-', label='data')
    plt.plot(xdata, D_sig_xx, 'g-', label='data')
    plt.plot(fit_data_x, func(fit_data_x, *popt), 'r-')
    #print(pcov,np.diag(pcov))
    plt.show()    
    #perr_phi = (np.sqrt(pcov[0,0]))
    #perr_so = (np.sqrt(pcov[1,1]))
    #perr_e = (np.sqrt(pcov[2,2]))
       
    l_phi = np.sqrt(h_bar/(4*e*popt[0]))/1e-9
    l_so = np.sqrt(h_bar/(4*e*popt[1]))/1e-9
    l_e = np.sqrt(h_bar/(4*e*B_e_mfp))/1e-9
    
    B_max = B_e_mfp
    print('l_phi = ' + str(l_phi)  + ' nm')
    print('l_so = ' + str(l_so) + ' nm')
    print('l_e = ' + str(l_e) + ' nm')
    print('B_max = ' + str(B_max) + ' T')
    return (l_phi, l_so, l_e)
    
    
    
    
def Hikami_fit_looper(date,start_num,le_array):
    """
    Loops over the Hikami fit function for multiple data sets
    Use with "top_gate_step_B_field_sweep" function.
    """
    tg_values =  np.linspace(-0.7,0.05,50)  #  np.linspace(-0.25,0.42241379,27)  # np.linspace(-0.47857143,0.12142857,29) #np.linspace(-0.653125,0.190625,19) # np.linspace(-0.653125,-0.23125,10)   # [-0.5,-0.4,-0.25,-0.15,0.05,0.0,0.05,0.15,0.25] #np.linspace(0.14285714,0.25,6) #np.linspace(-0.5,0.12142857,30)   #np.linspace(-0.5,0.25,36)
    L_phi = np.zeros(len(tg_values))
    L_so = np.zeros(len(tg_values))
    L_e = np.zeros(len(tg_values))
    L_phi_rev = np.zeros(len(tg_values))
    L_so_rev = np.zeros(len(tg_values))
    L_e_rev = np.zeros(len(tg_values))
    for i, item in enumerate(tg_values):  # Need to update coeffs depending on l_e
        coeffs = np.array(le_array)
        l_mfp = np.poly1d(coeffs)
        if Hikami_fit(date, start_num+i, current=True, direction=1, l_mfp=l_mfp(item)) == 'change':
            day = str(int(date[-2:])+1)  
            if len(day) < 2:
                day = '0' + day
            date = date[:-2] + day
            start_num = 1-i  
        print(start_num+i)
        l_phi, l_so, l_e = Hikami_fit(date, start_num+i, current=True, direction=1, l_mfp=l_mfp(item))
        L_phi[i] = l_phi
        L_so[i] = l_so
        L_e[i] = l_e
        l_phi_rev, l_so_rev, l_e_rev = Hikami_fit(date, start_num+i, current=True, direction=-1, l_mfp=l_mfp(item))
        L_phi_rev[i] = l_phi_rev
        L_so_rev[i] = l_so_rev
        L_e_rev[i] = l_e_rev
    
    
    plt.plot(tg_values[0:],L_phi[0:],'ro')
    plt.plot(tg_values,L_phi_rev,'rx')
    plt.ylabel('l_phi (nm)', color='r')
    plt.xlabel('V_tg (V)')
    #plt.xlim(0.138,0.255)
    #plt.ylim(0,2500)
    l_phi_fit = np.polyfit(tg_values,L_phi, 5)
    l_phi_poly = np.poly1d(l_phi_fit)
    #plt.plot(tg_values,l_phi_poly(tg_values),'b--')
    #l_phi_fit_rev = np.polyfit(tg_values,L_phi_rev, 5)
    #l_phi_poly_rev = np.poly1d(l_phi_fit_rev)
    #plt.plot(tg_values,l_phi_poly_rev(tg_values),'b-.')
    plt.show()
        
    plt.plot(tg_values[0:],L_so[0:],'bo')
    plt.plot(tg_values,L_so_rev,'bx')
    plt.ylabel('l_so (nm)', color='b')
    plt.xlabel('V_tg (V)')
    #plt.ylim(0,800)
    l_so_fit = np.polyfit(tg_values,L_so, 5)
    l_so_poly = np.poly1d(l_so_fit)
    #plt.plot(tg_values,l_so_poly(tg_values),'r--')
    #l_so_fit_rev = np.polyfit(tg_values,L_so_rev, 5)
    #l_so_poly_rev = np.poly1d(l_so_fit_rev)
    #plt.plot(tg_values,l_so_poly_rev(tg_values),'r-.')
    plt.show()
        
    plt.plot(tg_values[0:],L_e[0:],'go-')
    plt.plot(tg_values,L_e_rev,'gx-')
    plt.ylabel('l_e (nm)', color='g')
    plt.xlabel('V_tg (V)')
    #plt.ylim(100,500)
    plt.show()
    
    #Av_L_e = L_e 
    #Av_L_so = (L_so + L_so_rev)/2
    
    plt.plot(L_e[0:],L_so[0:],'ro', label='V < -0.42')
    #l_e_so_fit = np.polyfit(L_e[3:25],L_so[3:25], 1)
    #print(l_e_so_fit)
    #l_e_so_poly = np.poly1d(l_e_so_fit)
    #plt.plot(L_e[:29],l_e_so_poly(L_e[:29]),'b-')
    
    #plt.plot(L_e[25:35],L_so[25:35],'mx')
    
    #plt.plot(L_e[35:],L_so[35:],'go', label='V > -0.38')
    #hi_l_e_so_fit = np.polyfit(L_e[35:],L_so[35:], 1)
    #print(hi_l_e_so_fit)
    #hi_l_e_so_poly = np.poly1d(hi_l_e_so_fit)
    #plt.plot(L_e[35:],hi_l_e_so_poly(L_e[35:]),'b-')
    
    #plt.plot(L_e_rev[0:8],L_so_rev[0:8],'ro-.', label='V < -0.08 (density ~ 0.6e12 cm-2)')
    #plt.plot(L_e_rev[8:17],L_so_rev[8:17],'mo', label='mfp too long for good fit')
    #plt.plot(L_e_rev[17:],L_so_rev[17:],'go-.', label='V > 0.15 (density ~ 1.0e12 cm-2)')
    plt.xlabel('l_e (nm)')
    plt.ylabel('l_so (nm)')
    #plt.legend()
    #plt.ylim(0,1100)
    #plt.xlim(0,350)
    plt.show()
    
    plt.plot(L_e,L_phi,'ro-')
    plt.show()
  
    
    
    
    
def ILP_fit(date, num, current=False, width=50e-6, length=100e-6, direction=1, l_mfp=3.51e-7):
    """
    Fits weak (anti-)localization with the ILP formula 
    Use with "top_gate_step_B_field_sweep" function.
    """
    # Define the fitting functions
    e = 1.60217662e-19
    h_bar = 1.0545718e-34
    pi = 3.141
     
    B_e_mfp = h_bar/(4.*e*l_mfp**2)
        
    def func(B, B_phi, B_so, B_e=B_e_mfp):                   
        return 
    
    # Import the data
    data = open_data(date, num)
    if data == None:
        return 'change'
    xdata = (data.ami_field_set.ndarray.copy()*direction)+0.0001  # corrects slight offset
    ydata = data.SR860_1_X_preamp.ndarray.copy()     
    if current == True:
        curr = data.SR860_3_X_preamp.ndarray.copy()*-1e-6
    else:
        curr = 1e-9
    
    # Calculate conductivity
    # rho_xx is scaled by the width+length of the hall bar to get a sheet resistivity
    rho_xx = (ydata/curr) * (width/length)
    sig_xx = 1/rho_xx
    
    # Find Zero values (average)
    a = np.where(np.abs(xdata) < 0.0011)
    sig_xx_0 = sum(sig_xx[a])/len(sig_xx[a])
    
    # Calculate Magnetoconductivity
    D_sig_xx = sig_xx - sig_xx_0
    
    #plt.plot(xdata, D_sig_xx, 'b-', label='data')
    #plt.xlim(-0.005,0.005)
    plt.xlabel('B (T)')
    plt.ylabel('Delta conductance (S)')
    #plt.ylim(-0.00005,0.000003)
        
    # Specify range of data to be fit (> 0)
    fit_data_x = xdata[np.where(xdata > 0.0)]
    fit_data_y = D_sig_xx[np.where(xdata > 0.0)]
    
    B_max = B_e_mfp
    
    fit_data_x = fit_data_x[np.where(fit_data_x < B_max)]
    fit_data_y = fit_data_y[np.where(fit_data_x < B_max)]
       
    trial_values = np.array([1.1e-6,1.8e-6])#,1.5e-3])    
       
    popt, pcov = curve_fit(func, fit_data_x, fit_data_y, p0=trial_values, bounds=([0.0,0.0],[1.1e-1,1.1e-1]))
    plt.plot(fit_data_x, fit_data_y, 'b-', label='data')
    plt.plot(fit_data_x, func(fit_data_x, *popt), 'r-')
    #print(pcov,np.diag(pcov))
    plt.show()    
    #perr_phi = (np.sqrt(pcov[0,0]))
    #perr_so = (np.sqrt(pcov[1,1]))
    #perr_e = (np.sqrt(pcov[2,2]))
       
    l_phi = np.sqrt(h_bar/(4*e*popt[0]))/1e-9
    l_so = np.sqrt(h_bar/(4*e*popt[1]))/1e-9
    l_e = np.sqrt(h_bar/(4*e*B_e_mfp))/1e-9
    
    B_max = B_e_mfp
    print('l_phi = ' + str(l_phi)  + ' nm')
    print('l_so = ' + str(l_so) + ' nm')
    print('l_e = ' + str(l_e) + ' nm')
    print('B_max = ' + str(B_max) + ' T')
    return (l_phi, l_so, l_e)
    
    
    
    

def N_e_vs_mu(date, start_num, tomorrow_num, current=False, width=25e-6, lockin_set='top'):
    """
    Calculate and plot the mobility (and rho_xx) vs. the electron density 
    for the "top_gate_E_field_sweep_B_field_stepped" function. Also plots 
    current (amplitude, phase), rho_xx, n_e and mu as a function of applied 
    top-gate voltage
    """
    curr=1e-9
    np_field = np.array([-1.0,-0.5,-0.25,-0.05,0.25,0.5,1.0])          # work on this !!!! to specify fields automatically
    data_points = 200    
    pair_num = 5           # number of devices under test                                   # work on this !!!! to specify data points automatically
    
    
    Density_out = np.zeros(data_points)    
    Mobility_out = np.zeros(data_points)
    Res_out = np.zeros(data_points)
    I_out = np.zeros(data_points)
    Q_out = np.zeros(data_points)
    Vg_out = np.zeros(data_points)
    
    i = 0
    while i < data_points:   
        
        V_xx = np.zeros(np_field.shape[0])
        V_xy = np.zeros(np_field.shape[0])
        I = np.zeros(np_field.shape[0])
        Q = np.zeros(np_field.shape[0])
        V_g = np.zeros(np_field.shape[0])
        
        reset = False
    
        for j, field in enumerate(np_field):
            # Import the data
            if reset == True:
                mid_num = tomorrow_num - minus
                date1 = date2
            else:
                mid_num = start_num
                date1 = date
            num = (j*pair_num)+mid_num                                       # work on this !!!! to specify which PAIR
            data = open_data(date1, num)  
            print(date1,num, field)
            if data == None:
                day = str(int(date[-2:])+1)   
                print(day)
                if len(day) < 2:
                    day = '0' + day
                tomorrow_date = date[:-2] + day
                if int(day) > 31:
                    tomorrow_date = date[:-5] + '0' + str(int(date[-5:-3])+1) + '-01'
                date2 = tomorrow_date
                num = tomorrow_num  
                data = open_data(date2, num) 
                print(date2,num, field)
                reset = True
                minus = (j*pair_num)
            
            #print(i,j,field,data)
            
            if lockin_set == 'top':            
                V_xx[j] = (data.SR860_1_X_preamp[i]) 
                V_xy[j] = (data.SR860_2_X_preamp[i])
            else:
                V_xx[j] = (data.SR860_3_X_preamp[i]) 
                V_xy[j] = (data.SR860_4_X_preamp[i])
                    
            if current==True:
                I[j] = (1*data.SR860_3_X_preamp[i]/1e6)
                Q[j] = np.arctan(data.SR860_3_Y_preamp[i]/data.SR860_3_X_preamp[i])*180/3.141
                curr = I   
                
            V_g[j] = (data.yoko_t_voltage_set[i])
        
        # Calculate resistances
        rho_xy = V_xy/curr
        # rho_xx is scaled by the width+length of the hall bar to get a sheet resistivity
        rho_xx = (V_xx/curr) * (width/100e-6)
                        
        # Calculate density
        res = np.polyfit(np_field, rho_xy, 1)
        poly = np.poly1d(res)
        
        # Then the density is given by 1/(|e| dp/dB), in cm^2
        density = 1/(const.e * res[0])
        density *= 1e-4
                     
        # And the calculated mobility is given by mu=1/(rho_xx |e| ns)
        mu = 1/(rho_xx * const.e * density)
                
        # And let's quote the density slightly offset from 0, at 0.05T
        mob_ind = np.where(np.abs(np_field + 0.05) < 0.01)[0][0]
        mobility = mu[mob_ind]
                
        Density_out[i] = density
        Mobility_out[i] = mobility
        res_ind = np.where(np.abs(np_field + 0.05) < 0.01)[0][0]
        Res_out[i] = rho_xx[res_ind]
        I_out[i] = curr[res_ind]
        Q_out[i] = Q[res_ind]
        Vg_out[i] = np.average(V_g)                
        
        
        i=i+1
    
    
    
    if current==True:
        fig, ax1 = plt.subplots()
        ax1.plot(np.abs(Density_out),np.abs(Mobility_out),'rx', label = '{} um device'.format(1e6*width))
        ax1.set_ylabel('mobility (cm2/Vs)', color='r')
        ax1.set_xlabel('electron density (cm-2)')
        ax1.tick_params('y', colors='r')
        ax1.set_xlim(0,3.6e+12)
        ax1.set_ylim([0.0,4.6e+4])
        ax2 = ax1.twinx()
        ax2.plot(np.abs(Density_out),Res_out,'bx')
        ax2.set_xlabel('electron density (cm-2)')
        ax2.set_ylabel('rho_xx (Ohms/sq)', color='b')
        ax2.tick_params('y', colors='b')  
        ax2.set_ylim([0, 4000])
        fig.tight_layout()
        plt.legend()
        plt.show()
        name = "Mob_vs_Dens_{}_{}_".format(date,start_num)
        plt.savefig(name)
        
        fig2 = plt.plot(Vg_out, np.abs(Res_out), 'bx', label = '{} um device'.format(1e6*width))
        plt.ylabel('Sheet resistance (Ohms/sq)')
        plt.xlabel('applied top gate voltage (V)')
        #plt.xlim(-0.9,0.9)
        plt.ylim([100,4000000])
        plt.yscale('log')
        plt.show()
        
        fig3 = plt.plot(Vg_out, Q_out, 'gx', label = '{} um device'.format(1e6*width))
        plt.ylabel('Current phase (degs)')
        plt.xlabel('applied top gate voltage (V)')
        #plt.xlim(-0.9,0.9)
        plt.yscale('linear')
        plt.show()
        
        fig4, ax3 = plt.subplots()
        ax3.plot(Vg_out, np.abs(Density_out), 'mx')
        ax3.set_ylabel('electron density (cm-2)', color='m')
        ax3.set_xlabel('applied top gate voltage (V)')
        ax3.tick_params('y', colors='m')
        #ax3.set_xlim(-0.9,1.0)
        ax3.set_ylim(0,3.6e+12)
        ax4 = ax3.twinx()
        ax4.plot(Vg_out, np.abs(Mobility_out), 'cx', label = '{} um device'.format(1e6*width))
        ax4.set_ylabel('mobility (cm2/Vs)', color='c')
        ax4.tick_params('y', colors='c')  
        ax4.set_ylim(0,4.5e+4)
        fig4.tight_layout()
        plt.legend()
        plt.show()
        name = "Mob_vs_Dens_{}_{}_".format(date,start_num)
        plt.savefig(name)
        
        #i=0
        #while i < len(Vg_out):
        #    print(Vg_out[i], np.abs(Density_out[i]))
        #    i=i+1
            
        fig5 = plt.plot(Vg_out, I_out, 'bx', label = '{} um device'.format(1e6*width))
        plt.ylabel('Current (A)')
        plt.xlabel('applied top gate voltage (V)')
        #plt.xlim(-0.9,1.0)
        plt.show()
        
        e = 1.60217662e-19
        h_bar = 1.0545718e-34
        pi = math.pi
        
        def l_e(mu, n_e):
            return (h_bar/e)*mu*1.0e-4*np.sqrt(2.0*pi*n_e*1.0e+4)
        
        fig6 = plt.plot(Vg_out, l_e(np.abs(Mobility_out),np.abs(Density_out)), 'b-')
        
        Vg_out_fit = Vg_out[np.where(np.abs(Vg_out) < 0.7)]
        l_e_fit = l_e(np.abs(Mobility_out),np.abs(Density_out))[np.where(np.abs(Vg_out) < 0.7)]
        mfp_fit = np.polyfit(Vg_out_fit, l_e_fit, 17)
        poly = np.poly1d(mfp_fit)
        plt.plot(Vg_out_fit,poly(Vg_out_fit),'r--')
        plt.ylabel('l_e (m)')
        plt.xlabel('applied top gate voltage (V)')
        #plt.ylim(0,1e-5)
        plt.show()
        
        print(list(mfp_fit))
        
        Vg_values = Vg_out[np.where(np.abs(I_out)>0.9*I_out[np.argmax(np.abs(I_out))])]
        Mobility_out_values = Mobility_out[np.where(np.abs(I_out)>0.9*I_out[np.argmax(np.abs(I_out))])]
        Density_out_values = Density_out[np.where(np.abs(I_out)>0.9*I_out[np.argmax(np.abs(I_out))])]
        
        Vg_0 = np.argmin(np.abs(Vg_values))        
        mu_peak_index = np.argmax(np.abs(Mobility_out_values))
        
        mu_peak = "%7.4g" % Mobility_out_values[mu_peak_index]
        mu_peak_dens = "%7.4g" % Density_out_values[mu_peak_index]
        mu_peak_bias = "%7.4g" % Vg_values[mu_peak_index] 
        mu_0 = "%7.4g" % Mobility_out_values[Vg_0]
        Dens_0 = "%7.4g" % Density_out_values[Vg_0]
        
        #print('Gate bias at 0 (V) = ' + str(Vg_values[Vg_0]))
        
        print('zero bias mobility = ' + str(mu_0))
        print('zero bias density = ' + str(Dens_0))
        
        print('peak mobility = ' + str(mu_peak))
        print('bias at peak mobility = ' + str(mu_peak_bias))        
        print('denisty at peak mobility = ' + str(mu_peak_dens))
        
    else:
        plt.plot(np.abs(Density_out),np.abs(Mobility_out),'rx', label = '{} um device'.format(1e6*width))
        plt.xlabel('electron density (cm-2)')
        plt.ylabel('mobility (cm2/Vs)')
        plt.xlim(0,3.6e+12)
        plt.ylim(0.0,3.3e+4)
        plt.legend()
        plt.show()
        name = "Mob_vs_Dens_{}_{}_".format(date,start_num)
        plt.savefig(name)
        
    return #Density_out,Mobility_out
   
    
    
    
def invert_field(date, num):
    """
    Plot the data as a function of inverse applied magnetic field
    and calculate and plot the FFT of the oscillations. 
    Use with "inverse_B_field_sweep" function
    """
    # Define the fitting functions
    e = 1.60217662e-19
    h_bar = 1.0545718e-34
    pi = math.pi
          
    # Import the data
    data = open_data(date, num)
    xdata = data.ami_field_set.ndarray.copy()
    ydata = data.SR860_1_X_preamp.ndarray.copy()   
    zdata = data.SR860_2_X_preamp.ndarray.copy()
    curr = data.SR860_3_X_preamp.ndarray.copy()
    
   
    # Calculate rho_xx
    rho_xx = (ydata/(curr*-1e-6)) * (1/4)
    
    # Calculate rho_xy
    r_H = 1/((zdata/(curr*-1e-6))/(h_bar*pi/e**2)) # Factor of 2 !!!
    print((h_bar*2*pi/e**2))
         
    # Invert B
    inv_B = (1/np.abs(xdata))
    
    plt.plot(inv_B, rho_xx)
    
    yf = scipy.fftpack.fft(ydata)
    yf_rho = scipy.fftpack.fft(rho_xx)
    N = len(inv_B)
    T = np.abs((inv_B[-1]-inv_B[0])/(N-1))             
    xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    print(N,T,xdata[-1],xdata[0])

    fig, ax = plt.subplots()
    ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
    ax.set_ylim([0.0,1.0e-7])
    ax.set_xlim([5,50])
    ax.grid(color='g', linestyle='--', linewidth=0.1)
    
    #ax2 = ax.twinx()
    #ax2.plot(xf, 2.0/N * np.abs(yf_rho[:N//2]))
    #ax2.set_ylim([0.0,10.0])
    plt.show()
    
    plt.plot(xdata, r_H)
    plt.grid(color='g', linestyle='--', linewidth=0.1)
    plt.xlim(1.4,2.0)
    plt.ylim(-19,-13)
    plt.show()
    
    
    
def lockin_IV_plot(date, num, series_R=1.0e7):
    data = open_data(date, num)
    V_out = data.SR860_1_sine_outdc_set.ndarray.copy()
    V_in = data.SR860_1_X_preamp.ndarray.copy()   
    I_in = data.SR860_3_X_preamp.ndarray.copy()*-1.0e6
    
    dVdI = (V_in/I_in)
    I = V_out/series_R
    
    plt.plot(I,dVdI,'r-', label = '')
    plt.xlabel('current (A)')
    plt.ylabel('dV/dI (Ohms)')
    #plt.xlim(0,3.6e+12)
    plt.ylim(1.0e-9,1.075e-8)
    plt.legend()
    plt.show()    
    
    plt.plot(I,V,'ro', label = '')
    plt.xlabel('current (A)')
    plt.ylabel('V (Ohms)')
    #plt.xlim(0,3.6e+12)
    #plt.ylim(-2e-16,2e-16)
    plt.legend()
    plt.show()
    
       
    
    
def Dingle(date,num):
    return 
    
    
def multi_plotter():
    B_x, B_y = N_e_vs_mu('2018-05-16',30,1,current=True,width=25e-6,lockin_set='top')
    C_x, C_y = N_e_vs_mu('2018-05-16',31,2,current=True,width=25e-6,lockin_set='top')
    D_x, D_y = N_e_vs_mu('2018-05-16',32,3,current=True,width=25e-6,lockin_set='top')
    E_x, E_y = N_e_vs_mu('2018-05-16',33,4,current=True,width=25e-6,lockin_set='top')
    plt.plot(np.abs(B_x), np.abs(B_y), 'b-', label=' 7')
    plt.plot(np.abs(C_x), np.abs(C_y), 'g-', label=' 13')
    plt.plot(np.abs(D_x), np.abs(D_y), 'r-', label=' 18')
    plt.plot(np.abs(E_x), np.abs(E_y), 'm-', label=' 19')
    plt.xlim(0.3e+12,2.5e+12)
    plt.ylim([0.0,3.0e+4])
    plt.xlabel('electron density (cm-2)')
    plt.ylabel('mobility (cm2/Vs)')
    plt.legend()
    plt.show()
    
def rxx_rxy_plot(date, num,  sw, loc, width=50e-6):
    data = open_data(date, num)
    field = np.array(data.ami_field_set)
    I = np.array(data.SR860_3_X_preamp) * 1e-6
    Vxx = np.array(data.SR860_1_X_preamp)
    Vxy = np.array(data.SR860_2_X_preamp)
    
    I_mean = I.mean()
    Rxx = Vxx/I
    Rxy = Vxy/I
    
    # Calculate density
    res = np.polyfit(field, Rxy, 1)
    poly = np.poly1d(res)
    
    # Then the density is given by 1/(|e| dp/dB), in cm^2
    density = 1/(const.e * res[0])
    density *= 1e-4
                 
    # And the calculated mobility is given by mu=1/(rho_xx |e| ns)
    mu = 1/(Rxx  * (width/100e-6) * const.e * density)
    # And let's quote the density slightly offset from 0, at 0.05T
    mob_ind = np.where(np.abs(field + 0.05) < 0.01)[0][0]
    mobility = mu[mob_ind]
    print(f"Density: {density:.3e}")
    print(f"Mobility: {mobility:.3e}")
    
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Field (T)")
    
    ax1.plot(field, Rxx, 'b')
    ax1.set_ylabel("R_xx (Ohms)")
    ax1.tick_params('y', colors='b')
    
    ax2 = ax1.twinx()
    ax2.plot(field, Rxy, 'r')
    ax2.plot(field, poly(field), 'k-')
    ax2.set_ylabel("R_xy (Ohms)")
    ax2.tick_params('y', colors='r')
    
    fig.tight_layout()
    plt.title(f"Field Sweep, MB: {sw}, Loc: {loc}")
    plt.show()
    #plt.savefig(f"field_sweep_{date}_{num}.png")
    
def wl_scan_2d_plot(date, start, num_sweeps):
    data = open_data(date, start)
    
    yoko_array = np.zeros((num_sweeps,))
    data_arrays = []
    field_arrays = []
    field_arrays.append(np.array(data.ami_field_set))
    data_arrays.append(np.array(data.SR860_1_X_preamp))
    yoko_array[0] = data.metadata['station']['instruments']['yoko_t']['parameters']['voltage']['value']
    
    offs = 0
    vmin = float('inf')
    vmax = float('-inf')
    for i in range(1, num_sweeps):
        print(date, start+i-offs)
        data = open_data(date, start+i-offs)
        
        if data is None:
            date = datetime.datetime.strptime(date, "%Y-%m-%d")
            date += datetime.timedelta(days=1)
            date = date.strftime("%Y-%m-%d")
            offs = start+i-1
            data = open_data(date, start+i-offs)
            if data is None:
                print(f"Error opening data at: {date}, {start+i-offs}")
                break
        
        data_array = np.array(data.SR860_1_X_preamp)/(np.array(data.SR860_3_X_preamp).mean()*1e-6)
        data_arrays.append(data_array)
        field_arrays.append(np.array(data.ami_field_set))
        yoko_array[i] = data.metadata['station']['instruments']['yoko_t']['parameters']['voltage']['value']
        vmin = min(vmin, data.SR860_1_X_preamp.min())
        vmax = max(vmax, data.SR860_1_X_preamp.max())
    
    fig, ax = plt.subplots()
    ax.set_xlabel("TG Voltage (V)")
    ax.set_ylabel("Field")
    
    for tg_volt, field, val in zip(yoko_array, field_arrays, data_arrays):
        tg_volt_arr = np.full((field.shape[0],), tg_volt)
        ax.scatter(tg_volt_arr, field, c=val, s=1, cmap='plasma')
    plt.show()
    
    return (yoko_array, field_arrays, data_arrays)