# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:02:06 2018

@author: spauka
"""

import math
from pathlib import Path

import numpy as np
from numpy import ma

import scipy.constants as const
from scipy.optimize import curve_fit
from scipy.special import digamma
import scipy.fftpack

import pyqtgraph as qtplot
import matplotlib.pyplot as plt

import qcodes as qc
from qcodes.plots.colors import color_cycle

from data_utils import *
from plot_tools import *

plt.ion()

def calc_rho(field, data_xx, data_xy, curr=1e-9, width=50e-6, length=100e-6,
             field_center=1, name="Analyzed_Field_Sweep"):
    """
    Calculate rho_xx, rho_xy from given data sweeps:
    Parameters:
     - field: field setpoints
     - data_xx: V_xx DataArray
     - data_xy: V_xy DataArray
     - curr: Either the current as a number or as a DataArray. If a DataArray is given,
             the trace will be averaged to extract a mean current.
     - width/length: Width of the hall bar/distance between V_xx probes to extract
                     a sheet resistance
     - field_center:
     - name: Name to save output data
    """
    path = "data" / Path(name)

    np_field = field.ndarray.copy()
    np_xx = data_xx.ndarray.copy()
    np_xy = data_xy.ndarray.copy()

    if isinstance(curr, qc.DataArray):
        curr = np.average(-1*(curr.ndarray.copy())/1e6)
        print(curr)

    # Calculate resistances
    rho_xy = np_xy/curr
    # rho_xx is scaled by the width+length of the hall bar to get a sheet resistivity
    rho_xx = (np_xx/curr) * (width/length)

    # Calculate density
    # We want to do a fit between (field_center, -field_center) tesla as there is some extra structure
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
        set_arrays=(da_field, ), preset_data=rho_xy)
    da_rho_xy_fit = qc.DataArray(
        array_id="fit_Rho_XY", label="Rho XY", unit="Ohms",
        set_arrays=(da_reduced_field, ), preset_data=poly(np_field[min_ind:max_ind]))
    da_rho_xy_coef = qc.DataArray(
        array_id="coef_Rho_XY", label="Poly Coefficients",
        set_arrays=(da_poly, ), preset_data=res)
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
    fig, ax = plt.subplots()
    ax.plot(dataset.Field, dataset.Rho_XY, 'r')
    ax.plot(dataset.Field, dataset.fit_Rho_XY, 'k')
    ax.set_xlabel(format_ax(dataset.Field))
    ax.set_ylabel(format_ax(dataset.Rho_XY))
    ax.text(0.05, 0.95, "Using {} current<br>".format(qtplot.siFormat(curr, suffix="A")) +
                        "From a linear fit:<br>" +
                        "dρ/dB = {}<br>".format(qtplot.siFormat(res[0], suffix="Ω")) +
                        "n<sub>s</sub> = 1/(|e| dρ/dB) = {:e} cm<sup>-2</sup>".format(density),
            transform=ax.transAxes)
    fig.savefig(path/"rho_xy.png")

    # Showing Mobility analysis
    fig, ax = plt.subplots()
    ax.plot(dataset.Field, dataset.mu, c=color_cycle[5])
    ax.set_xlabel(format_ax(dataset.Field))
    ax.set_ylabel(format_ax(dataset.mu))
    ax.text(0.05, 0.95, "Mobility extracted from:<br>" +
                        "μ = 1/ρ<sub>xx</sub> |e| n<sub>s</sub>, with n<sub>s</sub>= {:.2e} cm<sup>-2</sup><br>".format(density) +
                        "And using W = {}, L = {}".format(qtplot.siFormat(width, suffix="m"),
                                                          qtplot.siFormat(length, suffix="m")),
            transform=ax.transAxes)
    fig.savefig(path/"mobility.png")

    return dataset


def leak_test_plot(date, num, width=50e-6):
    """
    Plot gate leakage sweep
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
    ax1.plot(V_g, V_xx, 'rx', label='{} um device'.format(1e6*width))
    ax1.set_xlabel('top gate voltage (V)')
    ax1.set_ylabel('V_xx (V)', color='r')
    ax1.tick_params('y', colors='r')
    ax1.set_ylim([0.0, 2e-5])

    ax2 = ax1.twinx()
    ax2.plot(V_g, np.abs(curr), 'b-')
    ax2.set_xlabel('top gate voltage (V)')
    ax2.set_ylabel('current (uA)', color='b')
    ax2.tick_params('y', colors='b')
    ax2.set_ylim([0, 0.0042])

    fig.tight_layout()
    plt.legend()
    plt.show()
    name = "leak_test_{}_{}_".format(date, num)
    plt.savefig(name)

def Hikami_fit(date, num, current=False, width=50e-6, length=100e-6, direction=1, l_mfp=1.1310514919623452e-07, B_max=0.015):
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
    if data is None:
        return 'change'
    xdata = (data.yoko_mag_current_set.ndarray.copy()*direction)/15.600624025
    ydata = data.SR860_1_X_preamp.ndarray.copy()*1e-3
    if current:
        curr = data.SR860_2_X_preamp.ndarray.copy()*1e-6
    else:
        curr = 1e-9

    # Calculate conductivity
    # rho_xx is scaled by the width and length of the hall bar to get a sheet resistivity
    rho_xx = (ydata/curr) * (width/length)
    sig_xx = 1/rho_xx

    # Find Zero values (average)
    a1 = np.where(np.abs(xdata) < 0.0001)    # select low field range
    b1 = np.argmax(sig_xx[a1])               # find sigma max in this rage
    c1 = a1[0]

    sig_xx_0 = sig_xx[c1[0]+b1]
    xdata = xdata - xdata[c1[0]+b1]          # corrects slight field offset

    # Calculate Magnetoconductivity
    D_sig_xx = sig_xx - sig_xx_0

    #plt.plot(xdata, D_sig_xx, 'b-', label='data')
    #plt.xlim(-0.005, 0.005)
    plt.xlabel('B (T)')
    plt.ylabel('Delta conductance (S)')
    #plt.ylim(-0.00002, 0.000015)

    # Specify range of data to be fit (> 0)
    fit_data_x = xdata[np.where(xdata > 0.0)]
    fit_data_y = D_sig_xx[np.where(xdata > 0.0)]

    #B_max = 0.015#B_e_mfp*0.75
    #if B_max > 0.003:
    #    B_max = 0.003

    fit_data_y = fit_data_y[np.where(fit_data_x < B_max)]
    fit_data_x = fit_data_x[np.where(fit_data_x < B_max)]

    trial_values = np.array([1.2e-6, 7.2e-7])#, 1.5e-3])

    popt, pcov = curve_fit(func, fit_data_x, fit_data_y, p0=trial_values, bounds=([0.0, 0.0], [1.1, 1.1]))
    plt.plot(fit_data_x, fit_data_y, 'b-', label='data')
    plt.plot(xdata, D_sig_xx, 'g.', label='data')
    plt.plot(fit_data_x, func(fit_data_x, *popt), 'r-')

    #### Trial functions
    #l_phi_trial = 2.2e-6
    #l_so_trial = 1.2e-6

    #B_phi_trial = h_bar/(4.*e*l_phi_trial**2)
    #B_so_trial = h_bar/(4.*e*l_so_trial**2)

    #plt.plot(fit_data_x, func(fit_data_x, B_phi_trial, B_so_trial), 'c-')
    ###

    B_phi_sd = np.sqrt(pcov[0, 0])
    B_so_sd = np.sqrt(pcov[1, 1])


    #print(pcov, np.diag(pcov))
    plt.show()

    ##is the fit any good?
    decision = input('is the fit good? y/n ')
    if decision == 'n':
        return 'try_again'

    #perr_phi = (np.sqrt(pcov[0, 0]))
    #perr_so = (np.sqrt(pcov[1, 1]))
    #perr_e = (np.sqrt(pcov[2, 2]))

    l_phi = np.sqrt(h_bar/(4*e*popt[0]))/1e-9
    l_so = np.sqrt(h_bar/(4*e*popt[1]))/1e-9
    l_e = np.sqrt(h_bar/(4*e*B_e_mfp))/1e-9

    l_phi_sd = (B_phi_sd*l_phi)/popt[0]
    l_so_sd = (B_so_sd*l_so)/popt[1]

    B_max = B_e_mfp
    print('l_phi = ' + str(l_phi)  + ' nm')
    print('l_so = ' + str(l_so) + ' nm')
    print('l_e = ' + str(l_e) + ' nm')
    print('B_max = ' + str(B_max) + ' T')

    alpha = (h_bar**2/(0.026*9.1e-31*l_so*1e-9))/e

    print('alpha = ' + str(alpha) + ' eV.m')

    return l_phi, l_so, l_e, l_phi_sd, l_so_sd




def Hikami_fit_looper(date, start_num, le_array):
    """
    Loops over the Hikami fit function for multiple data sets
    Use with "top_gate_step_B_field_sweep" function.
    """
    tg_values = np.linspace(-0.25, 0.5, 50)  #  np.linspace(-0.25, 0.42241379, 27)  # np.linspace(-0.47857143, 0.12142857, 29) #np.linspace(-0.653125, 0.190625, 19) # np.linspace(-0.653125, -0.23125, 10)   # [-0.5, -0.4, -0.25, -0.15, 0.05, 0.0, 0.05, 0.15, 0.25] #np.linspace(0.14285714, 0.25, 6) #np.linspace(-0.5, 0.12142857, 30)   #np.linspace(-0.5, 0.25, 36)
    L_phi = np.zeros(len(tg_values))
    L_so = np.zeros(len(tg_values))
    L_e = np.zeros(len(tg_values))
    L_phi_sd = np.zeros(len(tg_values))
    L_so_sd = np.zeros(len(tg_values))


    L_phi_rev = np.zeros(len(tg_values))
    L_so_rev = np.zeros(len(tg_values))
    L_e_rev = np.zeros(len(tg_values))
    L_phi_rev_sd = np.zeros(len(tg_values))
    L_so_rev_sd = np.zeros(len(tg_values))



    for i, item in enumerate(tg_values):  # Need to update coeffs depending on l_e
        B_max = 0.015
        coeffs = np.array(le_array)
        l_mfp = np.poly1d(coeffs)
        if Hikami_fit(date, start_num+i, current=True, direction=1, l_mfp=l_mfp(item), B_max=B_max) == 'change':
            day = str(int(date[-2:])+1)
            if len(day) < 2:
                day = '0' + day
            date = date[:-2] + day
            start_num = 1-i
        print(start_num+i)
        if Hikami_fit(date, start_num+i, current=True, direction=1, l_mfp=l_mfp(item), B_max=B_max) == 'try_again':
            B_max = float(input('maximum field to fit (T) = '))
            if Hikami_fit(date, start_num+i, current=True, direction=1, l_mfp=l_mfp(item), B_max=B_max) == 'try_again':
                skip = input('skip this plot? y/n ')
                if skip == 'y':
                    l_phi, l_so, l_e, l_phi_sd, l_so_sd = ('nan', )*5
                    l_phi_rev, l_so_rev, l_e_rev, l_phi_rev_sd, l_so_rev_sd = ('nan', )*5
                if skip == 'n':
                    l_phi, l_so, l_e, l_phi_sd, l_so_sd = Hikami_fit(date, start_num+i, current=True, direction=1, l_mfp=l_mfp(item), B_max=B_max)
        else:
            l_phi, l_so, l_e, l_phi_sd, l_so_sd = Hikami_fit(date, start_num+i, current=True, direction=1, l_mfp=l_mfp(item), B_max=B_max)
            l_phi_rev, l_so_rev, l_e_rev, l_phi_rev_sd, l_so_rev_sd = Hikami_fit(date, start_num+i, current=True, direction=-1, l_mfp=l_mfp(item), B_max=B_max)

        L_phi[i] = l_phi
        L_so[i] = l_so
        L_e[i] = l_e
        L_phi_sd[i] = l_phi_sd
        L_so_sd[i] = l_so_sd

        L_phi_rev[i] = l_phi_rev
        L_so_rev[i] = l_so_rev
        L_e_rev[i] = l_e_rev
        L_phi_rev_sd[i] = l_phi_rev_sd
        L_so_rev_sd[i] = l_so_rev_sd

    #plt.subplot(2, 2, 1)
    #plt.errorbar(tg_values[0:], L_phi[0:], yerr=L_phi_sd[0:], fmt='o', elinewidth=1)
    #plt.errorbar(tg_values, L_phi_rev, yerr=L_phi_rev_sd[0:], fmt='x', elinewidth=1)
    L_phi_av = (L_phi + L_phi_rev)/2
    L_phi_av_sd = np.sqrt(L_phi_sd**2+L_phi_rev_sd**2)
    plt.errorbar(tg_values[0:], L_phi_av[0:], yerr=L_phi_av_sd[0:], fmt='o', elinewidth=1)
    plt.ylabel('l_phi (nm)', color='r')
    plt.xlabel('V_tg (V)')
    plt.show()

    #plt.subplot(2, 2, 2)
    #plt.errorbar(tg_values[0:], L_so[0:], yerr=L_so_sd[0:], fmt='o', elinewidth=1)
    #plt.errorbar(tg_values, L_so_rev, yerr=L_so_rev_sd[0:], fmt='x', elinewidth=1)
    L_so_av = (L_so + L_so_rev)/2
    L_so_av_sd = np.sqrt(L_so_sd**2+L_so_rev_sd**2)
    plt.errorbar(tg_values, L_so_av, yerr=L_so_av_sd[0:], fmt='o', elinewidth=1)
    plt.ylabel('l_so (nm)', color='b')
    plt.xlabel('V_tg (V)')
    plt.ylim(0, 2300)
    plt.show()

    #plt.subplot(2, 2, 3)
    plt.plot(tg_values[0:], L_e[0:], 'go-')
    plt.ylabel('l_e (nm)', color='g')
    plt.xlabel('V_tg (V)')
    plt.show()


    #plt.subplot(2, 2, 4)
    plt.errorbar(L_e[4:15], L_so_av[4:15], yerr=L_so_av_sd[4:15], fmt='o', elinewidth=1, label='whole range')
    l_e_so_fit = np.polyfit(L_e[4:15], L_so_av[4:15], 1)
    l_e_so_poly = np.poly1d(l_e_so_fit)
    plt.plot(L_e[4:15], l_e_so_poly(L_e[4:15]), 'b-')

    plt.errorbar(L_e[15:22], L_so_av[15:22], yerr=L_so_av_sd[15:22], fmt='o', elinewidth=1, label='whole range')

    plt.errorbar(L_e[22:50], L_so_av[22:50], yerr=L_so_av_sd[22:50], fmt='o', elinewidth=1, label='whole range')
    l_e_so_fit = np.polyfit(L_e[22:50], L_so_av[22:50], 1)
    l_e_so_poly = np.poly1d(l_e_so_fit)
    plt.plot(L_e[22:50], l_e_so_poly(L_e[22:50]), 'm-')

    plt.xlabel('l_e (nm)')
    plt.ylabel('l_so (nm)')
    #plt.legend()
    plt.show()

    idx = np.isfinite(L_e) & np.isfinite(L_phi_av)
    plt.errorbar(L_e[idx], L_phi_av[idx], yerr=L_phi_av_sd[idx], fmt='o', elinewidth=1)
    plt.xlabel('l_e (nm)')
    plt.ylabel('l_phi (nm)')
    l_e_phi_fit = np.polyfit(L_e[idx], L_phi_av[idx], 2)
    l_e_phi_poly = np.poly1d(l_e_phi_fit)
    plt.plot(L_e[idx], l_e_phi_poly(L_e[idx]), label='tentative quadratic fit...')
    plt.legend()
    plt.show()





def ILP_fit(date, num, current=False, width=50e-6, length=100e-6, direction=1, l_mfp=3.51e-7):
    """
    Fits weak (anti-)localization with the ILP formula
    Use with "top_gate_step_B_field_sweep" function.
    """
    # Define the fitting functions
    e = scipy.constants.e
    h_bar = scipy.constants.hbar

    B_e_mfp = h_bar/(4.*e*l_mfp**2)

    def func(B, B_phi, B_so, B_e=B_e_mfp):
        return

    # Import the data
    data = open_data(date, num)
    if data is None:
        return 'change'
    xdata = (data.ami_field_set.ndarray.copy()*direction)+0.0001  # corrects slight offset
    ydata = data.SR860_1_X_preamp.ndarray.copy()
    if current:
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
    #plt.xlim(-0.005, 0.005)
    plt.xlabel('B (T)')
    plt.ylabel('Delta conductance (S)')
    #plt.ylim(-0.00005, 0.000003)

    # Specify range of data to be fit (> 0)
    fit_data_x = xdata[np.where(xdata > 0.0)]
    fit_data_y = D_sig_xx[np.where(xdata > 0.0)]

    B_max = B_e_mfp

    fit_data_x = fit_data_x[np.where(fit_data_x < B_max)]
    fit_data_y = fit_data_y[np.where(fit_data_x < B_max)]

    trial_values = np.array([1.1e-6, 1.8e-6])#, 1.5e-3])

    popt, pcov = curve_fit(func, fit_data_x, fit_data_y, p0=trial_values, bounds=([0.0, 0.0], [1.1e-1, 1.1e-1]))
    plt.plot(fit_data_x, fit_data_y, 'b-', label='data')
    plt.plot(fit_data_x, func(fit_data_x, *popt), 'r-')
    #print(pcov, np.diag(pcov))
    plt.show()
    #perr_phi = (np.sqrt(pcov[0, 0]))
    #perr_so = (np.sqrt(pcov[1, 1]))
    #perr_e = (np.sqrt(pcov[2, 2]))

    l_phi = np.sqrt(h_bar/(4*e*popt[0]))/1e-9
    l_so = np.sqrt(h_bar/(4*e*popt[1]))/1e-9
    l_e = np.sqrt(h_bar/(4*e*B_e_mfp))/1e-9

    B_max = B_e_mfp
    print('l_phi = ' + str(l_phi)  + ' nm')
    print('l_so = ' + str(l_so) + ' nm')
    print('l_e = ' + str(l_e) + ' nm')
    print('B_max = ' + str(B_max) + ' T')
    return (l_phi, l_so, l_e)

def Exponent(date, start_num, tomorrow_num, current=False, width=50e-6, lockin_set='top'):
    """
    Calculate and plot the the exponent relating mu to n as a funtion of n
    """
    curr = 1e-9
    np_field = np.array([-1.0, -0.5, -0.25, -0.05, 0.25, 0.5, 1.0])          # work on this !!!! to specify fields automatically
    data_points = 200
    pair_num = 5          # number of devices under test                                   # work on this !!!! to specify data points automatically


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

        for j in range(np_field.size):
            # Import the data
            if reset:
                mid_num = tomorrow_num - minus
                date1 = date2
            else:
                mid_num = start_num
                date1 = date
            num = (j*pair_num)+mid_num                                       # work on this !!!! to specify which PAIR
            data = open_data(date1, num)
            #print(date1, num, field)
            if data is None:
                day = str(int(date[-2:])+1)
                #print(day)
                if len(day) < 2:
                    day = '0' + day
                tomorrow_date = date[:-2] + day
                if int(day) > 31:
                    tomorrow_date = date[:-5] + '0' + str(int(date[-5:-3])+1) + '-01'
                date2 = tomorrow_date
                num = tomorrow_num
                data = open_data(date2, num)
                #print(date2, num, field)
                reset = True
                minus = (j*pair_num)

            #print(i, j, field, data)

            if lockin_set == 'top':
                V_xx[j] = (data.SR860_1_X_preamp[i])
                V_xy[j] = (data.SR860_2_X_preamp[i])
            else:
                V_xx[j] = (data.SR860_3_X_preamp[i])
                V_xy[j] = (data.SR860_4_X_preamp[i])

            if current:
                I[j] = (-1*data.SR860_3_X_preamp[i]/1e6)
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


        i = i+1



    if current:

        Mobility_out_exps = Mobility_out[np.where(np.abs(I_out) > 0.9*I_out[np.argmax(np.abs(I_out))])]
        Density_out_exps = Density_out[np.where(np.abs(I_out) > 0.9*I_out[np.argmax(np.abs(I_out))])]


        #Mobility_out_exps_new = np.zeros(Mobility_out_exps.shape)
        ###remove outliers
        #for i, item in enumerate(Mobility_out_exps):
        #    j = i
        #    if j < 2:
        #        j=2
        #    if j > Mobility_out_exps.shape[0] - 3:
        #        j=Mobility_out_exps.shape[0] - 3
        #    section = np.concatenate((Mobility_out_exps[j-2:j-1], Mobility_out_exps[j+1:j+3]))
        #    if i == 14:
        #        print(section)
        #        print(i, item, np.abs(item - (np.mean(section))), np.std(section))
        #    if np.abs(item - (np.mean(section)) > np.std(section)):
        #        print('point ' + str(i) + ' expunged at n = ' + str(Density_out_exps[i]))
        #        Mobility_out_exps_new[i] = 'nan'
        #    else:
        #        Mobility_out_exps_new[i] = item

        #print(np.log(np.abs(Mobility_out_exps_new)))

        dens_start = 0.3e12
        dens_end = 1.0e12
        sample_size = 0.15e12

        alpha = []
        alpha_sd = []
        n = []

        start_idx = np.abs((np.abs(Density_out_exps) - dens_start)).argmin()
        end_idx = np.abs((np.abs(Density_out_exps) - dens_end)).argmin()
        sample_idx = 15# np.abs((np.abs(Density_out_exps) - (dens_start + sample_size))).argmin() - start_idx

        print(start_idx, end_idx, sample_idx)

        value = start_idx
        while value < end_idx - sample_idx:

            mu_fit_data_x2 = Density_out_exps[value:value+sample_idx]
            mu_fit_data_y2 = Mobility_out_exps[value:value+sample_idx]

            mu_fit_data_y2[mu_fit_data_y2 == 0] = 'nan'

            log_x = np.log(np.abs(mu_fit_data_x2))
            log_y = np.log(np.abs(mu_fit_data_y2))

            idx = np.isfinite(log_x) & np.isfinite(log_y)
            grad, cov = np.polyfit(log_x[idx], log_y[idx], 1, cov=True)

            plt.plot(np.abs(mu_fit_data_x2), mu_fit_data_y2)
            plt.show()

            plt.plot(log_x, log_y, 'b.')
            fit = np.poly1d(grad)
            plt.plot(log_x, fit(log_x), 'r-')
            plt.show()
            alpha_val = grad[0]
            alpha_sd_val = np.sqrt(cov[0, 0])
            alpha.append(alpha_val)
            alpha_sd.append(alpha_sd_val)
            n.append(np.abs(Density_out_exps[value]))
            value = value + 1

        alpha_grads = []
        for i, value in enumerate(alpha):
            j = i
            if j == 0:
                j = 1
            if j == len(alpha)-1:
                j = len(alpha)-2
            alpha_grad = (alpha[j+1]-alpha[j-1])/(n[j+1]-n[j-1])
            alpha_grads.append(alpha_grad)


        plt.axhline(y=0, color='k', linestyle='-', linewidth=5.0)
        plt.plot(n, alpha_grads, 'g-.')
        plt.ylim(-0.8e-11, 0.1e-11)
        plt.xlabel('electron density (cm-2)')
        plt.ylabel('dalpha/dn')
        plt.show()


        plt.errorbar(n, alpha, yerr=alpha_sd, fmt='o', elinewidth=1)
        plt.legend()
        plt.axhline(y=0.5, color='r', linestyle='--')
        plt.axhline(y=1.0, color='r', linestyle='--')
        plt.axhline(y=1.5, color='r', linestyle='--')
        plt.axhline(y=1.7, color='r', linestyle='--')

        #plt.ylim(0.4, 2.6)
        plt.xlabel('electron density (cm-2)')
        plt.ylabel('alpha')
        plt.show()


    else:
        print('you are probably using the wrong function!')









def invert_field(date, num):
    """
    Plot the data as a function of inverse applied magnetic field
    and calculate and plot the FFT of the oscillations.
    Use with "inverse_B_field_sweep" function
    """
    # Define the fitting functions
    e = const.e
    h_bar = const.hbar
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

    print(N, T, xdata[-1], xdata[0])

    fig, ax = plt.subplots()
    ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
    ax.set_ylim([0.0, 1.0e-7])
    ax.set_xlim([5, 50])
    ax.grid(color='g', linestyle='--', linewidth=0.1)

    #ax2 = ax.twinx()
    #ax2.plot(xf, 2.0/N * np.abs(yf_rho[:N//2]))
    #ax2.set_ylim([0.0, 10.0])
    plt.show()

    plt.plot(xdata, r_H)
    plt.grid(color='g', linestyle='--', linewidth=0.1)
    plt.xlim(1.4, 2.0)
    plt.ylim(-19, -13)
    plt.show()



def lockin_IV_plot(date, num, series_R=1.0e7):
    data = open_data(date, num)
    V_out = data.SR860_1_sine_outdc_set.ndarray.copy()
    V_in = data.SR860_1_X_preamp.ndarray.copy()
    I_in = data.SR860_3_X_preamp.ndarray.copy()*-1.0e6

    dVdI = (V_in/I_in)
    I = V_out/series_R

    plt.plot(I, dVdI, 'r-', label='')
    plt.xlabel('current (A)')
    plt.ylabel('dV/dI (Ohms)')
    #plt.xlim(0, 3.6e+12)
    plt.ylim(1.0e-9, 1.075e-8)
    plt.legend()
    plt.show()

    plt.plot(I, V_in, 'ro', label='')
    plt.xlabel('current (A)')
    plt.ylabel('V (Ohms)')
    #plt.xlim(0, 3.6e+12)
    #plt.ylim(-2e-16, 2e-16)
    plt.legend()
    plt.show()


def rxx_rxy_plot(date, num, loc="UNK", width=50e-6):
    data = open_data(date, num)
    field = np.array(data.ami_field_set)
    I = np.array(data.SR860_3_X_preamp) * 1e-6
    Vxx = np.array(data.SR860_1_X_preamp)
    Vxy = np.array(data.SR860_2_X_preamp)
    sw = instr_param_val(data, 'md', 'select')

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
    """
    Plot weak localization scans

    Note: this version of the function requires a fixed sweep range, split
    between coarse and fine ranges. This was true for sweeps run after 2018-10-05
    (66bca5a)
    """
    # Create arrays to hold data
    yoko_array = np.zeros((num_sweeps, ))
    data_arrays = []
    field_arrays = []

    # Load data arrays
    for i, data in enumerate(open_data_sequence(date, start, num_sweeps)):
        data_array = np.array(data.SR860_1_X_preamp)/(np.array(data.SR860_2_X_preamp)*1e-6)
        data_arrays.append(data_array)
        field_arrays.append(np.array(data.yoko_mag_current_set)/15.600624025)
        yoko_array[i] = instr_param_val(data, 'yoko_t', 'voltage')

    fine_field_array = np.linspace(0.0015, -0.0015, 601)
    fine_data_array = np.zeros((601, num_sweeps))

    coarse_field_array = np.linspace(0.011, -0.011, 221)
    coarse_data_array = np.full((221, num_sweeps), float('NaN'))

    for i, (field, val) in enumerate(zip(field_arrays, data_arrays)):
        numpoints = field.size
        midindex = (numpoints-600)//2
        startind = np.where(np.isclose(field[0], coarse_field_array))[0][0]

        fine_data_array[:, i] = val[midindex:midindex + 601]
        coarse_data_array[startind:startind+midindex, i] = val[:midindex]
        coarse_data_array[125:125+midindex, i] = val[midindex + 601:]

        assert all(np.isclose(fine_field_array, field[midindex:midindex + 601]))
        assert all(np.isclose(coarse_field_array[startind:startind+midindex], field[:midindex]))
        assert all(np.isclose(coarse_field_array[125:125+midindex], field[midindex + 601:]))

        dmin = min(fine_data_array[:, i].min(), np.nanmin(coarse_data_array[:, i]))
        dmax = max(fine_data_array[:, i].max(), np.nanmax(coarse_data_array[:, i]))
        fine_data_array[:, i] = (fine_data_array[:, i] - dmin)/(dmax-dmin)
        coarse_data_array[:, i] = (coarse_data_array[:, i] - dmin)/(dmax-dmin)

        #ax.scatter(tg_volt_arr, field, c=val, s=5, cmap='plasma')
    masked_coarse_data_array = ma.masked_invalid(coarse_data_array)

    # Figure out corrected field/voltage arrays
    x_step = yoko_array[1] - yoko_array[0]
    p_voltage_array = np.arange(yoko_array[0] - (x_step/2),
                                yoko_array[-1] + 2*(x_step/2),
                                x_step)
    pc_field_array = np.arange(0.011 + 0.0001/2,
                               -0.011 - 2*0.0001/2,
                               -0.0001)
    pf_field_array = np.arange(0.0015 + 0.000005/2,
                               -0.0015 - 2*0.000005/2,
                               -0.000005)
    assert p_voltage_array.size == num_sweeps+1

    fig, ax = plt.subplots()
    ax.set_xlabel("TG Voltage (V)")
    ax.set_ylabel("Field")
    ax.pcolormesh(p_voltage_array, pc_field_array, masked_coarse_data_array, cmap='plasma')
    ax.pcolormesh(p_voltage_array, pf_field_array, fine_data_array, cmap='plasma')
    #ax.set_xlim(min(yoko_array), max(yoko_array))
    ax.set_ylim(-0.011, 0.011)
    fig.tight_layout()
    plt.show()

    return (yoko_array, field_arrays, data_arrays)
