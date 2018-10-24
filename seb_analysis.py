# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:02:06 2018

@author: spauka
"""

import math
from time import sleep
from pathlib import Path
import datetime
from collections import namedtuple

import numpy as np
from numpy import ma

import scipy.constants as const
from scipy.optimize import curve_fit
from scipy.special import digamma
from scipy.integrate import quad
import scipy.fftpack

import pyqtgraph as qtplot
import matplotlib.pyplot as plt

import qcodes as qc
from qcodes.plots.colors import color_cycle

from data_utils import *
from plot_tools import *

plt.ion()

AnalyzedDevice = namedtuple("AnalyzedDevice", ['vg', 'density', 'mobility', 'res', 'curr', 'phase'])

def N_e_vs_mu_all(date, start_num, stop_num, width=50e-6, current=True, lockin_set='top', plot=True, subset=None):
    """
    Calculate and plot the mobility vs. electron density for all devices
    run with "top_gate_E_field_sweep_B_field_stepped".

    As the title suggests, this function assumes that we cycle through
    devices (i.e. A,B,C,D,E) with constant field, with TG swept.

    If we only want to analyze a subset of the data, set subset to the devices
    that should be analyzed (i.e. ('B', 'C'))
    """
    data_arrays = open_data_sequence(date, start_num, stop_num)

    num_devices, devices = detect_cycle([inst_param_val(data, 'md', 'select') for data in data_arrays])

    # And fill in the fields
    fields = [inst_param_val(data, 'ami', 'field') for data in data_arrays[::num_devices]]
    # Let's also double check that each sweep has the correct field set
    for i in range(1, num_devices):
        fields_check = [inst_param_val(data, 'ami', 'field') for data in data_arrays[i::num_devices]]
        if fields != fields_check:
            raise ValueError("The fields between devices are not equal. "
                             f"Unequal fields were: {fields}, {fields_check}")

    print(devices, fields)

    analyzed_devices = {}
    for i, device in enumerate(devices):
        if subset is not None and device not in subset:
            continue
        analyzed_devices[device] = N_e_vs_mu(data_arrays[i::num_devices], fields, device,
                                             lockin_set=lockin_set, width=width, current=current, plot=plot)
    return analyzed_devices

def N_e_vs_mu(data_arrays, fields, sw, current=False, width=50e-6, lockin_set='top', plot=True):
    """
    Calculate and plot the mobility (and rho_xx) vs. the electron density
    for the "top_gate_E_field_sweep_B_field_stepped" function. Also plots
    current (amplitude, phase), rho_xx, n_e and mu as a function of applied
    top-gate voltage
    """
    curr = 1e-9
    if current and lockin_set != 'top':
        raise RuntimeError("Lockin set must be top if current is being read (third lockin must be current...)")
    np_field = np.array(fields)

    data_points = data_arrays[0].yoko_t_voltage_set.size

    Density_out = np.zeros(data_points)
    Mobility_out = np.zeros(data_points)
    Res_out = np.zeros(data_points)
    Curr_out = np.zeros(data_points)
    Phase_out = np.zeros(data_points)
    Vg_out = np.zeros(data_points)

    for i in range(data_points):
        V_xx = np.zeros(np_field.size)
        V_xy = np.zeros(np_field.size)
        Curr = np.zeros(np_field.size)
        Phase = np.zeros(np_field.size)
        V_g = np.zeros(np_field.size)

        for j, data in enumerate(data_arrays):
            if lockin_set == 'top':
                V_xx[j] = (data.SR860_1_X_preamp[i])
                V_xy[j] = (data.SR860_2_X_preamp[i])
            else:
                V_xx[j] = (data.SR860_3_X_preamp[i])
                V_xy[j] = (data.SR860_4_X_preamp[i])

            if current:
                Curr[j] = (-1*data.SR860_3_X_preamp[i]/1e6)
                Phase[j] = np.arctan(data.SR860_3_Y_preamp[i]/data.SR860_3_X_preamp[i])*180/3.141
                curr = Curr

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
        mob_ind = np.where(np.abs(np_field) - 0.05 < 0.005)[0][0]
        mobility = mu[mob_ind]

        Density_out[i] = np.abs(density)
        Mobility_out[i] = np.abs(mobility)
        res_ind = np.where(np.abs(np_field) - 0.05 < 0.01)[0][0]
        Res_out[i] = np.abs(rho_xx[res_ind])
        if current:
            Curr_out[i] = np.abs(Curr[res_ind])
            Phase_out[i] = Phase[res_ind]
        else:
            Curr_out[i] = curr
            Phase_out[i] = np.nan
        Vg_out[i] = np.average(V_g)

    if current:
        if plot:
            density_mobility_plot(Density_out, Mobility_out, Res_out, title=f"Position {sw}")
            sheet_resistance_plot(Vg_out, Res_out, title=f"Position {sw}")
            phase_plot(Vg_out, Phase_out, title=f"Position {sw}")
            density_mobility_voltage_plot(Vg_out, Density_out, Mobility_out, title=f"Position {sw}")
            current_vg_plot(Vg_out, Curr_out, title=f"Position {sw}")
            mfp_plot(Vg_out, Mobility_out, Density_out, title=f"Position {sw}")

        Vg_values = Vg_out[np.where(np.abs(Curr_out) > 0.9*Curr_out[np.argmax(np.abs(Curr_out))])]
        Mobility_out_values = Mobility_out[np.where(np.abs(Curr_out) > 0.9*Curr_out[np.argmax(np.abs(Curr_out))])]
        Density_out_values = Density_out[np.where(np.abs(Curr_out) > 0.9*Curr_out[np.argmax(np.abs(Curr_out))])]

        Vg_0 = np.argmin(np.abs(Vg_values))
        mu_peak_index = np.argmax(np.abs(Mobility_out_values))

        mu_peak = "%7.4g" % Mobility_out_values[mu_peak_index]
        mu_peak_dens = "%7.4g" % Density_out_values[mu_peak_index]
        mu_peak_bias = "%7.4g" % Vg_values[mu_peak_index]
        mu_0 = "%7.4g" % Mobility_out_values[Vg_0]
        Dens_0 = "%7.4g" % Density_out_values[Vg_0]

        print('Gate bias at 0 (V) = ' + str(Vg_values[Vg_0]))

        print('zero bias mobility = ' + str(mu_0))
        print('zero bias density = ' + str(Dens_0))
        print('peak mobility = ' + str(mu_peak))
        print('bias at peak mobility = ' + str(mu_peak_bias))
        print('denisty at peak mobility = ' + str(mu_peak_dens))

    else:
        if plot:
            density_mobility_plot(Density_out, Mobility_out, Res_out, title=f"Position {sw}")
            sheet_resistance_plot(Vg_out, Res_out, title=f"Position {sw}")
            density_mobility_voltage_plot(Vg_out, Density_out, Mobility_out, title=f"Position {sw}")
            mfp_plot(Vg_out, Mobility_out, Density_out, title=f"Position {sw}")

    # Package device data into a named tuple
    device_data = AnalyzedDevice(Vg_out, Density_out, Mobility_out, Res_out, Curr_out, Phase_out)
    return device_data

def extract_peak_mu(Vg, density, mobility, current, dev="", pos=""):
    """
    Extract peak mobility by using a polynomial fit, and taking an average of values
    around the turning point at the peak of the polynomial
    """
    #valid_range = np.intersect1d(np.where(current > 0.75e-9), np.where(0.5e12 < density < 3e12))
    valid_range = np.where(np.all((current > 0.75e-9, density > 0.5e12, density < 3e12, mobility > 5000, Vg > -0.7), axis=0))
    Vg = Vg[valid_range]
    if Vg.size < 10:
        print("Not enough data for this dataset")
        return
    density = density[valid_range]
    mobility = mobility[valid_range]

    vg_mob_fit = np.polyfit(Vg, mobility, 12)
    vg_mob_poly = np.poly1d(vg_mob_fit)

    den_sort_ind = np.argsort(density)
    m_density = density[den_sort_ind]
    m_mobility = mobility[den_sort_ind]
    den_mob_fit = np.polyfit(m_density, m_mobility, 12)
    den_mob_poly = np.poly1d(den_mob_fit)

    vg_max_mob = np.argmax(vg_mob_poly(Vg))
    m1 = ""
    m1 += f"Maximum mobility at Vg={Vg[vg_max_mob]:.3}: {mobility[vg_max_mob-1:vg_max_mob+2].mean():.0f}\n"
    m1 += f"Density at the above is: {density[vg_max_mob-1:vg_max_mob+2].mean():.3e}"
    avg_mob = mobility[vg_max_mob-1:vg_max_mob+2].mean()
    avg_den = density[vg_max_mob-1:vg_max_mob+2].mean()
    avg_vg = Vg[vg_max_mob]
    print(m1)

    den_max_mob = np.argmax(den_mob_poly(m_density))
    m2 = ""
    m2 += f"Maximum mobility at density={m_density[den_max_mob-1:den_max_mob+2].mean():.3e}: "
    m2 += f"{m_mobility[den_max_mob-1:den_max_mob+2].mean():.0f}\n"
    m2 += f"Vg at the above is: {Vg[den_sort_ind[den_max_mob]]:.3}"
    avg_mob += m_mobility[den_max_mob-1:den_max_mob+2].mean()
    avg_den += m_density[den_max_mob-1:den_max_mob+2].mean()
    avg_vg += Vg[den_sort_ind[den_max_mob]]
    print(m2)

    print("---")
    print(f"Avg Mobility: {avg_mob/2:.0f}")
    print(f"Avg Density: {avg_den/2:.3e}")
    print(f"Avg Vg: {avg_vg/2:.3}")
    print("---")

    fig, ax = plt.subplots()
    ax.plot(Vg, mobility, 'cx')
    ax.plot(Vg, vg_mob_poly(Vg), 'k--')
    ax.set_xlabel("Top Gate Voltage (V)")
    ax.set_ylabel("Mobility (cm^2/V s)")
    ax.set_title(f"Device {dev}, Pos: {pos}")
    ax.text(0.95, 0.95, m1, transform=ax.transAxes, verticalalignment="top", horizontalalignment="right")
    fig.savefig(f"dev_{dev}_vg_mob", bbox_inches="tight")

    fig, ax2 = plt.subplots()
    ax2.plot(m_density, m_mobility, 'mx')
    ax2.plot(m_density, den_mob_poly(m_density), 'k--')
    ax2.set_title(f"Device {dev}, Pos: {pos}")
    ax2.set_xlabel("Density (cm^-2)")
    ax2.set_ylabel("Mobility (cm^2/V s)")
    ax2.text(0.95, 0.95, m2, transform=ax2.transAxes, verticalalignment="top", horizontalalignment="right")
    fig.savefig(f"dev_{dev}_den_mob", bbox_inches="tight")
