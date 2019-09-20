from collections import defaultdict
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colors
import scipy.constants as const
from scipy.optimize import curve_fit
import scipy.stats

from qcodes.dataset.data_set import load_by_id
from data_utils import detect_cycle, open_data_sequence

# Initialize matplotlib
plt.rc('text', usetex=False)
plt.rc('text.latex', preamble=r"\usepackage{siunitx}")
#plt.rcParams['font.sans-serif'] = ['Helvetica', 'sans-serif']
plt.ion()

def density_mobility_plot(electron_density, mobility, rho_xx=None, width=50e-6, title=""):
    fig, ax1 = plt.subplots()

    # Plot mobility
    ax1.plot(np.abs(electron_density), np.abs(mobility), 'ro', alpha=1, label='{} um device'.format(1e6*width))
    ax1.set_title(title)
    ax1.set_ylabel(r"Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")#, color='r')
    ax1.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    #ax1.tick_params('y', colors='r')
    ax1.set_xlim(0.3e12, 2.5e+12)
    ax1.set_ylim([0, 4.7e+4])

    # Plot rho_xx, if given
    if rho_xx is not None:
        ax2 = ax1.twinx()
        ax2.plot(np.abs(electron_density), np.abs(rho_xx), 'bx')
        ax2.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
        ax2.set_ylabel(r'rhoxx $\left(\Omega^{-1}\right)$', color='b')
        ax2.tick_params('y', colors='b')
        ax2.set_ylim([0, 4000])
    fig.tight_layout()
    return fig

def density_mobility_log_plot(electron_density, mobility, width=50e-6, title=""):
    # Plot mobility
    fig = density_mobility_plot(electron_density, mobility, width=width, title=title)
    ax = fig.axes[0]

    # Set log scale
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1e11, 5e12)
    ax.set_ylim(1e3, 1e5)

    # fit power law
    filt_idx = np.where(np.all((electron_density < 7e11, electron_density > 3e11, mobility > 8000), axis=0))
    electron_density = electron_density[filt_idx]
    mobility = mobility[filt_idx]

    # Plot the part which will be used
    ax.plot(electron_density, mobility, 'ob', alpha=0.5)

    def power_law(x, A, alpha):
        return A*(x**alpha)
    popt, _ = curve_fit(power_law, electron_density, mobility, p0=(1e-21, 2.5), maxfev=1800)
    print(popt)
    log_d = np.log10(electron_density)
    log_m = np.log10(mobility)
    p = np.polyfit(log_d, log_m, 1)
    print(p)
    pf = np.poly1d(p)

    ax.plot(electron_density, 10**pf(log_d), 'k-')
    ax.plot(electron_density, power_law(electron_density, *popt), 'y-')

    fig, ax = plt.subplots()
    ax.plot(log_d, log_m, 'xr', alpha=0.6)
    ax.plot(log_d, pf(log_d), 'k.')
    ax.plot(log_d, np.log10(power_law(electron_density, *popt)), 'y-')

def mobility_voltage_plot(Vg, mobility, title=""):
    fig, ax = plt.subplots()

    # Plot mobility
    ax.plot(Vg, mobility, 'mx')
    ax.set_ylabel(r"Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$", color='r')
    ax.set_xlabel('applied top gate voltage (V)')
    ax.tick_params('y', colors='c')
    ax.set_ylim(0, 4.5e4)
    ax.set_title(title)
    fig.tight_layout()
    return fig

def find_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

def density_mobility_voltage_plot(Vg, electron_density, mobility=None, width=50e-6, title="", highlight=None):
    fig, ax1 = plt.subplots()
    # Plot electron density
    ax1.plot(Vg, np.abs(electron_density), 'mx')
    ax1.set_ylabel(r"Density $\left(\si{\per\square\centi\meter}\right)$", color='m')
    ax1.set_xlabel('applied top gate voltage (V)')
    ax1.tick_params('y', colors='m')
    ax1.set_ylim(0, 3.0e+12)

    # Plot mobility
    if mobility is not None:
        # Filter mobility by density
        filt = np.where(electron_density > 0.2e12)

        ax2 = ax1.twinx()
        ax2.plot(Vg[filt], np.abs(mobility[filt]), 'co', label='{} um device'.format(1e6*width))
        ax2.set_ylabel(r"Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$", color='c')
        ax2.tick_params('y', colors='c')
        ax2.set_ylim(0, 4.7e+4)

    if highlight is not None:
        points = np.array(tuple(find_nearest(Vg, h_vg) for h_vg in highlight))
        print(points)
        ax1.plot(Vg[points], np.abs(electron_density[points]), 'o')
        ax2.plot(Vg[points], np.abs(mobility[points]), 'o')

    ax1.set_xlim(-0.7, 0.5)
    ax1.set_title(title)
    fig.tight_layout()
    return fig

def sheet_resistance_plot(Vg, rho_xx, width=50e-6, title=""):
    fig, ax = plt.subplots()
    ax.plot(Vg, np.abs(rho_xx), 'bx', label='{} um device'.format(1e6*width))
    ax.set_title(title)
    ax.set_ylabel(r'Sheet Resistance $\left(\Omega^{-1}\right)$', color='b')
    ax.set_xlabel('applied top gate voltage (V)')
    #ax.set_xlim(-0.9, 0.9)
    ax.set_ylim([100, 4000000])
    ax.set_yscale('log')
    return fig

def phase_plot(Vg, phase, width=50e-6, title=""):
    fig, ax = plt.subplots()
    ax.plot(Vg, phase, 'gx', label='{} um device'.format(1e6*width))
    ax.set_title(title)
    ax.set_ylabel('Current phase (degs)')
    ax.set_xlabel('applied top gate voltage (V)')
    #ax.set_xlim(-0.9, 0.9)
    ax.set_yscale('linear')
    return fig

def current_vg_plot(Vg, current, title=""):
    fig, ax = plt.subplots()
    ax.plot(Vg, current, 'bx')
    ax.set_title(title)
    ax.set_ylabel('Current (A)')
    ax.set_xlabel('applied top gate voltage (V)')
    #ax.set_xlim(-0.9, 1.0)
    return fig

def mfp_plot(Vg, mobility, density, title=""):
    def l_e(mu, n_e):
        e = const.e
        h_bar = const.hbar
        pi = const.pi
        return (h_bar/e)*mu*1.0e-4*np.sqrt(2.0*pi*n_e*1.0e+4)

    fig, ax = plt.subplots()
    ax.plot(Vg, l_e(np.abs(mobility), np.abs(density)), 'b-')
    Vg_out_fit = Vg[np.where(np.abs(Vg) < 0.7)]
    l_e_fit = l_e(np.abs(mobility), np.abs(density))[np.where(np.abs(Vg) < 0.7)]
    mfp_fit = np.polyfit(Vg_out_fit, l_e_fit, 17)
    poly = np.poly1d(mfp_fit)
    ax.plot(Vg_out_fit, poly(Vg_out_fit), 'r--')
    ax.set_ylabel('l_e (m)')
    ax.set_xlabel('applied top gate voltage (V)')
    ax.set_title(title)
    print(list(mfp_fit))
    #ax.set_ylim(0, 1e-5)
    return fig

def format_ax(data):
    return f"{data.label} ({data.unit})"

def plot_all(analyzed_devices, plot_func, extract, recolor=True):
    """
    Plot all the analyzed devices with the plotting function plot_func.
    Extract the parameters listed in extract.

    For example:
    >>> plot_all(analyzed_devices, density_mobility_plot, ("density", "mobility", "res"))
    """
    keys = tuple(analyzed_devices.keys())
    data = defaultdict(list)
    for key in keys:
        for param in extract:
            data[param].append(getattr(analyzed_devices[key], param))

    # Format the data in the correct format
    for param in extract:
        data[param] = np.array(data[param]).T

    # Run the plot function
    fig = plot_func(*data.values(), title="All Devices")
    ax = fig.axes[0]

    # If we need to recolor, then do that now
    if recolor:
        for color, line in zip(colors.TABLEAU_COLORS, ax.get_lines()):
            line.set_color(color)

    # Add a legend
    ax.legend(keys)

def plot_den_mob_analysis():
    data = np.genfromtxt("DenMobAnalysis.txt", delimiter=",", skip_header=True, converters={1: str, 2: str}, encoding="utf-8")

    fig, ax = plt.subplots()

    c = {'Near': 'g', 'Far': 'r'}
    m = {'Near': 'o', 'Far': '^'}
    treatments = {"A": ("No Pretreatment\nTMA/H2O", 3),
                  "B": ("TMA Prepulse\nTMA/H2O", 1),
                  "C": ("TMA Prepulse\nTMA/O3", 2),
                  "D": ("ArH Plasma\nTMA/H2O", 4),
                  "E": ("ArH Plasma\nTMA/O3", 0)
                  }
    # treatments = {"A": ("Treatment A", 3),
    #               "B": ("Treatment B", 1),
    #               "C": ("Treatment C", 2),
    #               "D": ("Treatment D", 4),
    #               "E": ("Treatment E", 0)
    #               }

    # Plot points
    for loc in c:
        x_vals = [treatments[d[1]][1] for d in data if d[1] in treatments and d[2] == loc]
        y_vals = [d[6] for d in data if d[1] in treatments and d[2] == loc]
        ax.scatter(x_vals, y_vals, c=c[loc], marker=m[loc], alpha=1.0, label=loc)

    # Plot Trendline
    x_vals = np.ones(len(treatments)).cumsum() - 1
    y_near_vals = np.ndarray(len(treatments))
    y_far_vals = np.ndarray(len(treatments))
    for treatment in treatments:
        idx = treatments[treatment][1]
        y_near_vals[idx] = np.mean(tuple(d[6] for d in data if d[1] == treatment and d[2] == 'Near'))
        y_far_vals[idx] = np.mean(tuple(d[6] for d in data if d[1] == treatment and d[2] == 'Far'))
    ax.plot(x_vals, y_near_vals, 'g--', alpha=0.6)
    ax.plot(x_vals, y_far_vals, 'r--', alpha=0.6)

    x_ticks = tuple(x[0] for x in sorted(treatments.values(), key=lambda x: x[1]))
    print(x_ticks)

    ax.set_xticks(np.arange(5))
    ax.set_xticklabels(x_ticks)
    ax.set_ylabel(r"Peak Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_title("Peak Mobility for Various Treatments")
    ax.legend()
    fig.tight_layout()
    fig.savefig("Peak Mobility Treatment.png")
    plt.show()

def plot_den_mob_scatter():
    data = np.genfromtxt("DenMobAnalysis.txt", delimiter=",", skip_header=True, converters={1: str, 2: str}, encoding="utf-8")

    fig, ax = plt.subplots()

    c = {'A': 'tab:red', 'B': 'tab:orange', 'C': 'tab:olive', 'D': 'tab:green', 'E': 'tab:blue'}
    x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/H2O", "ArH Plasma\nTMA/O3")
    treat = {a: b for a, b in zip(c.keys(), x_ticks)}
    for treatment in c:
        x_vals = [d[7] for d in data if d[1] == treatment]
        y_vals = [d[6] for d in data if d[1] == treatment]

        ax.scatter(x_vals, y_vals, c=c[treatment], alpha=0.6, label=treat[treatment])

    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Peak Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_xlim(0.18e12, 2.19e12)
    ax.set_ylim(3000, 48000)
    ax.set_title("Density v. Peak Mobility")
    ax.legend()
    fig.tight_layout()
    fig.savefig("Peak Mobility Density.png")
    plt.show()

def plot_vg_mob_scatter():
    data = np.genfromtxt("DenMobAnalysis.txt", delimiter=",", skip_header=True, converters={1: str, 2: str}, encoding="utf-8")

    fig, ax = plt.subplots()

    c = {'A': 'tab:red', 'B': 'tab:orange', 'C': 'tab:olive', 'D': 'tab:green', 'E': 'tab:blue'}
    x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/H2O", "ArH Plasma\nTMA/O3")
    treat = {a: b for a, b in zip(c.keys(), x_ticks)}
    all_xvals = []
    all_yvals = []
    for treatment in c:
        x_vals = [d[8] for d in data if d[1] == treatment]
        y_vals = [d[6] for d in data if d[1] == treatment]
        all_xvals.extend(x_vals)
        all_yvals.extend(y_vals)

        ax.scatter(x_vals, y_vals, c=c[treatment], alpha=0.6, label=treat[treatment])
    m, b, r_value, _, _ = scipy.stats.linregress(all_xvals, all_yvals)
    x_lin = np.linspace(-0.5, 0.1, 100)
    p1d = np.poly1d((m, b))
    ax.plot(x_lin, p1d(x_lin), 'k--', alpha=0.6)

    ax.set_xlabel(r"Top Gate Voltage $\left(\si{\volt}\right)$")
    ax.set_ylabel(r"Peak Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_title("Gate Voltage v. Peak Mobility")
    ax.legend()
    ax.text(0.95, 0.05, f"$r^2 = {r_value:.2}$",
            multialignment='left', bbox=dict(edgecolor='k', joinstyle='round', facecolor='w'),
            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig("Peak Mobility Voltage.png")
    plt.show()

def plot_den_vg_trace(analyzed_data):
    fig, ax = plt.subplots()

    c = {'A': 'tab:red', 'B': 'tab:orange', 'C': 'tab:olive', 'D': 'tab:blue', 'E': 'tab:green'}
    #x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/O3", "ArH Plasma\nTMA/H2O")
    x_ticks = ("Treatment A", "Treatment B", "Treatment C", "Treatment D", "Treatment E")
    treat = {a: b for a, b in zip(c.keys(), x_ticks)}
    scatters = []
    for device in c:
        density = analyzed_data[device].density
        mobility = analyzed_data[device].mobility
        current = analyzed_data[device].curr
        vg = analyzed_data[device].vg

        valid_range = np.where(np.all((current > 3.75e-9, density > 0.2e12, density < 3e12, mobility > 5000), axis=0))
        density = density[valid_range]
        mobility = mobility[valid_range]
        current = current[valid_range]
        vg = vg[valid_range]

        max_mob = np.argmax(mobility)

        s = ax.scatter(vg, density, c=c[device], alpha=0.4, label=treat[device])
        ax.scatter(vg[max_mob], density[max_mob], s=75, c='k', marker='x', alpha=1)
        scatters.append(s)

    ax.set_xlabel(r"Top Gate Voltage $\left(\si{\volt}\right)$")
    ax.set_ylabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_title("Density v. Gate Voltage")
    ax.legend(scatters, x_ticks, ncol=2)
    ax.text(0.95, 0.05, "Peak mobility for each trace\nis marked with an x",
            multialignment='left', bbox=dict(edgecolor='k', joinstyle='round', facecolor='w'),
            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig("Density Voltage.png")
    plt.show()

def plot_zero_volt_mobility_trace(analyzed_data):
    fig, ax = plt.subplots()

    c = {'A': 'tab:red', 'B': 'tab:orange', 'C': 'tab:olive', 'D': 'tab:blue', 'E': 'tab:green'}
    #x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/O3", "ArH Plasma\nTMA/H2O")
    x_ticks = ("Treatment A", "Treatment B", "Treatment C", "Treatment D", "Treatment E")
    treat = {a: b for a, b in zip(c.keys(), x_ticks)}
    scatters = []
    for device in c:
        density = analyzed_data[device].density
        mobility = analyzed_data[device].mobility
        current = analyzed_data[device].curr
        vg = analyzed_data[device].vg

        valid_range = np.where(np.all((current > 3.75e-9, density > 0.2e12, density < 3e12, mobility > 5000), axis=0))
        density = density[valid_range]
        mobility = mobility[valid_range]
        current = current[valid_range]
        vg = vg[valid_range]

        vg_zero = np.where(np.isclose(0, vg, atol=1e-3))[0][0]
        max_mob = np.argmax(mobility)

        s = ax.scatter(density[vg_zero], mobility[vg_zero], c=c[device], alpha=1, label=treat[device])
        ax.scatter(density[max_mob], mobility[max_mob], c=c[device], alpha=1, label=treat[device])
        ax.annotate("", (density[max_mob], mobility[max_mob]), (density[vg_zero], mobility[vg_zero]),
                    arrowprops=dict(headwidth=10, headlength=10, width=0.5, color=c[device], alpha=0.5))
        scatters.append(s)

    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_ylim(5000, 60000)
    ax.set_title("Density/Mobility at Zero Vg to Peak Mobility")
    ax.legend(scatters, x_ticks, ncol=2)
    fig.tight_layout()
    fig.savefig("Density Mobility Zero.png")
    plt.show()

def plot_SOI(fname):
    data = np.genfromtxt(fname, delimiter=",", skip_header=True)

    fig, ax = plt.subplots()
    m_e = 0.023 * const.electron_mass
    density = data[:, 7]
    l_so = data[:, 6] / const.nano
    alpha = data[:, 0] / const.angstrom
    l_mfp = data[:, 4] / const.nano
    l_e = data[:, 3] / const.nano

    scatter = ax.scatter(density, alpha, c=l_mfp, alpha=0.6)
    fig.colorbar(scatter, label=r"Electron Mean Free Path $l_{mfp} \left(\si{\nano\meter}\right)$")
    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Spin Orbit Parameter $\alpha \left(\si{\electronvolt\angstrom}\right)$")
    ax.set_title("Spin Orbit Coupling Parameter")
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    scatter = ax.scatter(l_e, l_so, c=density, cmap="magma", alpha=0.6)
    fig.colorbar(scatter, label=r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_xlabel(r"Electron Mean Free Path $l_{mfp} \left(\si{\nano\meter}\right)$")
    ax.set_ylabel(r"Spin Orbit Length $l_{so} \left(\si{\nano\meter}\right)$")
    ax.set_title("Spin Precession Mechanisms")
    fig.tight_layout()
    plt.show()

def plot_SOI_full(fname):
    data = np.genfromtxt(fname, delimiter=",", skip_header=True)

    m_e = 0.023 * const.electron_mass
    density = data[:, 12]
    l_so_r = data[:, 9] / const.nano
    l_so_d = data[:, 10] / const.nano
    l_so = data[:, 11] / const.nano
    alpha = data[:, 1] / const.angstrom
    gamma = data[:, 5] / const.angstrom
    l_mfp = data[:, 7] / const.nano
    l_e = data[:, 6] / const.nano
    v_f = const.hbar * np.sqrt(2 * const.pi * density/(const.centi**2))/m_e
    tau_e = l_e * const.nano / v_f / const.pico
    tau_soi = data[:, 16] / const.pico
    print(tau_e)

    # Spin orbit length
    fig, ax = plt.subplots()
    ax.scatter(density, l_so, alpha=0.6)
    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Spin Orbit Length $\l_{so} \left(\si{\nano\meter}\right)$")
    ax.set_title("Spin Orbit Length")
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    ax.scatter(density, l_so, alpha=0.6, label="Total SO Length")
    ax.scatter(density, l_so_d, alpha=0.6, label="Dresselhaus SO Length")
    ax.scatter(density, l_so_r, alpha=0.6, label="Rashba SO Length")
    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Spin Orbit Length $\l_{so} \left(\si{\nano\meter}\right)$")
    ax.set_title("Spin Orbit Length")
    ax.legend()
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    scatter = ax.scatter(density, alpha, c=l_mfp, alpha=0.6)
    fig.colorbar(scatter, label=r"Electron Mean Free Path $l_{mfp} \left(\si{\nano\meter}\right)$")
    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Rashba Spin Orbit Parameter $\alpha \left(\si{\electronvolt\angstrom}\right)$")
    ax.set_title("Rashba Spin Orbit Coupling Parameter")
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    scatter = ax.scatter(l_e, l_so, c=density, cmap="magma", alpha=0.6)
    fig.colorbar(scatter, label=r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_xlabel(r"Elastic Scattering Length $l_{mfp} \left(\si{\nano\meter}\right)$")
    ax.set_ylabel(r"Spin Relaxation Length $l_{so} \left(\si{\nano\meter}\right)$")
    ax.set_title("Spin Precession Mechanisms")
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    E_so = (np.sqrt(2)*const.hbar)/(l_so * const.nano * l_e * const.nano)
    E_so = E_so/const.electron_volt
    print(E_so)
    scatter = ax.scatter(density, E_so, c=l_mfp, alpha=0.6)
    fig.colorbar(scatter, label=r"Electron Mean Free Path $l_{mfp} \left(\si{\nano\meter}\right)$")
    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Spin Orbit Energy $E_{so} \left(\si{\electronvolt}\right)$")
    ax.set_title("Spin Precession Mechanisms")
    fig.tight_layout()
    plt.show()

def plot_parabola():
    x = np.ones(21).cumsum() - 11
    y = x**2

    fig, ax = plt.subplots()
    ax.plot(x, y, 'r-')
    fig.tight_layout()
    plt.show()

def plot_rxx_rxy(dataset, loc="top"):
    array_names = tuple(dataset.arrays.keys())

    # Check for a preamp
    preamps = tuple(array.endswith("preamp") for array in array_names if array.startswith("SR860"))
    if all(preamps):
        has_preamp = True
    elif any(preamps):
        raise ValueError("Some datasets have preamps but not all! Can't deal with that yet")
    else:
        has_preamp = False

    x_data = dataset.yoko_t_voltage_set
    if loc == "top":
        y_data_1 = dataset.SR860_1_X_preamp
        y_data_2 = dataset.SR860_2_X_preamp
    elif loc == "bot":
        y_data_1 = dataset.SR860_3_X_preamp
        y_data_2 = dataset.SR860_4_X_preamp
    else:
        raise ValueError(f"Invalid location {loc}")

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(x_data, y_data_1, 'r-')
    ax2.plot(x_data, y_data_2, 'g-')
    fig.tight_layout()
    plt.show()

def plot_dataset(dataset):
    # Check if we have a three-point (Current, Rxx, Rxy),
    # or a four-point (Rxx_1, Rxy_1, Rxx_2, Rxy_2) trace
    array_names = tuple(dataset.arrays.keys())
    if any(array.startswith("SR860_4_X") for array in array_names):
        plot_rxx_rxy(dataset, "top")
        plot_rxx_rxy(dataset, "bot")
    else:
        plot_rxx_rxy(dataset, "top")
        #plot_current(dataset)

def plot_peak_density_mobility(analyzed_data, position="any", scale="small"):
    """
    This function plots the gate voltage at which peak mobility is achieved.
    """
    exclude_ids = ()
    c = {(None, 'H2O'): 'tab:red', ('tma', 'H2O'): 'tab:orange', ('tma', 'O3'): 'tab:olive', ('plasma', 'H2O'): 'tab:blue', ('plasma', 'O3'): 'tab:green'}
    m = {(None, 'H2O'): 'x', ('tma', 'H2O'): '^', ('tma', 'O3'): 'v', ('plasma', 'H2O'): 'd', ('plasma', 'O3'): 'o'}
    # Extract a list of peak mobilities and voltages
    colors = []
    markers = []
    densities = []
    mobilities = []
    for sid, sample in analyzed_data.items():
        if sid in exclude_ids:
            continue
        if position != "any":
            if sample.position != position:
                continue
        # Get the measured values from the measurements for this sample
        for measurement in sample.measurements:
            colors.append(c[(sample.treatment, sample.ald["oxidiser"])])
            markers.append(m[(sample.treatment, sample.ald["oxidiser"])])
            densities.append(measurement.den_at_peak)
            mobilities.append(measurement.mob_at_peak)

        # Pick a mob/den trace to plot for the "large" scale
        if sample.treatment is None and sample.ald is not None and sample.ald["oxidiser"] == "H2O":
            density = sample.measurements[-1].analyzed.density
            mobility = sample.measurements[-1].analyzed.mobility

    fig, ax = plt.subplots()
    for den, mob, col, mark in zip(densities, mobilities, colors, markers):
        ax.scatter(x=den, y=mob, c=col, marker=mark, alpha=1)
    if scale == "large":
        ax.scatter(density, mobility, c="tab:grey", marker='.')
    ax.set_xlabel(r"Density at Peak $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Peak Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    if scale == "large":
        ax.axvline(0.8e12)
        ax.axvline(1e12)
        ax.set_xlim(0.3e12, 2.5e+12)
        ax.set_ylim(5000, 50000)
    elif scale == "small":
        ax.set_xlim(0.8e12, 1e+12)
        ax.set_ylim(15000, 50000)
    fig.tight_layout()
    plt.show()

def plot_peak_mobility_gate_voltage(analyzed_data, position="any"):
    """
    This function plots the gate voltage at which peak mobility is achieved.
    """
    exclude_ids = ()
    c = {(None, 'H2O'): 'tab:red', ('tma', 'H2O'): 'tab:orange', ('tma', 'O3'): 'tab:olive', ('plasma', 'H2O'): 'tab:blue', ('plasma', 'O3'): 'tab:green'}
    m = {(None, 'H2O'): 'x', ('tma', 'H2O'): '^', ('tma', 'O3'): 'v', ('plasma', 'H2O'): 'd', ('plasma', 'O3'): 'o'}
    # Extract a list of peak mobilities and voltages
    colors = []
    markers = []
    voltages = []
    mobilities = []
    for sid, sample in analyzed_data.items():
        if sid in exclude_ids:
            continue
        if position != "any":
            if sample.position != position:
                continue
        # Get the measured values from the measurements for this sample
        for measurement in sample.measurements:
            colors.append(c[(sample.treatment, sample.ald["oxidiser"])])
            markers.append(m[(sample.treatment, sample.ald["oxidiser"])])
            voltages.append(measurement.vg_at_peak)
            mobilities.append(measurement.mob_at_peak)

    fig, ax = plt.subplots()
    for volt, mob, col, mark in zip(voltages, mobilities, colors, markers):
        ax.scatter(x=volt, y=mob, c=col, marker=mark, alpha=1)
    ax.set_xlabel(r"Gate Voltage at Peak Mobility $\left(\si{\volt}\right)$")
    ax.set_ylabel(r"Peak Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_ylim(15000, 50000)
    fig.tight_layout()
    plt.show()

def plot_zero_density_peak_mobility(analyzed_data, position="any"):
    """
    This function plots the density at zero voltage vs peak density achieved
    """
    exclude_ids = ()
    c = {(None, 'H2O'): 'tab:red', ('tma', 'H2O'): 'tab:orange', ('tma', 'O3'): 'tab:olive', ('plasma', 'H2O'): 'tab:blue', ('plasma', 'O3'): 'tab:green'}
    m = {(None, 'H2O'): 'x', ('tma', 'H2O'): '^', ('tma', 'O3'): 'v', ('plasma', 'H2O'): 'd', ('plasma', 'O3'): 'o'}
    # Extract a list of peak mobilities and voltages
    colors = []
    markers = []
    densities = []
    mobilities = []
    for sid, sample in analyzed_data.items():
        if sid in exclude_ids:
            continue
        if position != "any":
            if sample.position != position:
                continue
        # Get the measured values from the measurements for this sample
        for measurement in sample.measurements:
            colors.append(c[(sample.treatment, sample.ald["oxidiser"])])
            markers.append(m[(sample.treatment, sample.ald["oxidiser"])])
            densities.append(measurement.den_at_zero)
            mobilities.append(measurement.mob_at_peak)

    linfit = np.polyfit(densities, mobilities, 1)

    fig, ax = plt.subplots()
    for dens, mob, col, mark in zip(densities, mobilities, colors, markers):
        ax.scatter(x=dens, y=mob, c=col, marker=mark, alpha=1)
    xvals = np.linspace(np.min(densities), np.max(densities), num=10)
    print(xvals)
    ax.plot(xvals, np.poly1d(linfit)(xvals), 'k--')
    ax.set_xlabel(r"Density at $V_g = 0$ $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Peak Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_ylim(15000, 50000)
    fig.tight_layout()
    plt.show()

def plot_den_mob_trace(data, title="Mobility Density"):
    fig, ax = plt.subplots()

    exclude_ids = ('ID_1', 'ID_9')

    c = {(None, 'H2O'): 'tab:red', ('tma', 'H2O'): 'tab:orange', ('tma', 'O3'): 'tab:olive', ('plasma', 'H2O'): 'tab:blue', ('plasma', 'O3'): 'tab:green'}
    m = {(None, 'H2O'): 'x', ('tma', 'H2O'): '^', ('tma', 'O3'): 'v', ('plasma', 'H2O'): 'd', ('plasma', 'O3'): 'o'}
    x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/H2O", "ArH Plasma\nTMA/O3")
    treat = {a: b for a, b in zip(c.keys(), x_ticks)}
    for sid, sample in data.items():
        if sid in exclude_ids:
            continue
        if (sample.treatment, sample.ald["oxidiser"]) not in c or sample.position == "far":
            continue
        density = sample.measurements[-1].analyzed.density
        mobility = sample.measurements[-1].analyzed.mobility
        current = sample.measurements[-1].analyzed.curr
        color = c[(sample.treatment, sample.ald["oxidiser"])]
        marker = m[(sample.treatment, sample.ald["oxidiser"])]
        treatment = treat[(sample.treatment, sample.ald["oxidiser"])]

        valid_range = np.where(np.all((current > 3.25e-9, density > 0.2e12, density < 3e12, mobility > 5000), axis=0))
        density = density[valid_range]
        mobility = mobility[valid_range]
        current = current[valid_range]

        ax.scatter(density, mobility, c=color, marker=marker, alpha=1, label=treatment)

    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig("Mobility Density.png")
    plt.show()

def plot_hall_sweep_hmob():
    datasets = (1054, 1055, 1056, 1057, 1058)

    field_data = np.zeros(0)
    R_xx_data = np.zeros(0)
    R_xy_data = np.zeros(0)

    for dataset in datasets:
        dataset = load_by_id(dataset)
        param_data = dataset.get_parameter_data("R_xx", "R_xy")
        field_data = np.append(field_data, param_data["R_xx"]["mag_GRPZ_field"])
        R_xx_data = np.append(R_xx_data, param_data["R_xx"]["R_xx"])
        R_xy_data = np.append(R_xy_data, param_data["R_xy"]["R_xy"])

    max_x, min_x = np.where(field_data > 7.2)[0][0], np.where(field_data > 6.6)[0][0]
    print(min_x, max_x)
    print(f"Rho_XY (nu = 6) = {np.average(R_xy_data[min_x:max_x])}")
    smallfield = np.where(field_data > 0.05)[0][0]
    print(smallfield)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(field_data, R_xy_data, 'b-', label="Rxy")
    ax2.plot(field_data, R_xx_data, 'r-', label="Rxx")

    # Plot plateaus
    vk_const = const.value("von Klitzing constant")
    for i in range(6, 12+1, 2):
        ax1.axhline(vk_const/i)
        print(f"At (nu = {i}): {vk_const/i:.2f}")

    # Plot theoretical locations
    density = np.polyfit(field_data[:500], R_xy_data[:500], 1)
    print(f"Fitting to field: {field_data[500]}.")
    density = 1/(const.e * density[0])
    print(f"Extracted Density: {density*1e-4:e} cm^2/Vs")
    filling_factors = (2*const.pi*density*const.hbar)/(const.e)
    for i in range(6, 12+1, 2):
        ax1.axvline(filling_factors/i)
    print(R_xx_data[11])
    sq = 1 / 5 # Width / length
    mu = 1/(np.average(R_xx_data[11-2:11+3]) * sq * const.e * density * 1e-4)
    print(f"Extracted Mobility: {mu}")

    fig.legend()
    ax1.set_xlabel(r"Field $\left(\si{\tesla}\right)$")
    ax1.set_ylabel(r"$\textrm{R}_{xy}$ $\left(\si{\ohm}\right)$")
    ax2.set_ylabel(r"$\textrm{R}_{xx}$ $\left(\si{\ohm}\right)$")
    fig.tight_layout()
    fig.show()

def dec_sin(inv_field, offs, A, f, phi, decay):
    return A * np.exp(-inv_field*decay) * np.sin(f*inv_field + phi) + offs

def plot_hall_sweep_lmob():
    datasets = (1045,)

    field_data = np.zeros(0)
    R_xx_data = np.zeros(0)
    R_xy_data = np.zeros(0)

    for dataset in datasets:
        dataset = load_by_id(dataset)
        param_data = dataset.get_parameter_data("R_xx", "R_xy")
        field_data = np.append(field_data, param_data["R_xx"]["mag_GRPZ_field"])
        R_xx_data = np.append(R_xx_data, param_data["R_xx"]["R_xx"])
        R_xy_data = np.append(R_xy_data, param_data["R_xy"]["R_xy"])
    field_data = field_data[:-100]
    R_xx_data = R_xx_data[:-100]
    R_xy_data = R_xy_data[:-100]

    field_data = 1/field_data

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(field_data, R_xy_data, 'b-', label="Rxy")
    ax2.plot(field_data, R_xx_data, 'r-', label="Rxx")

    res, _ = curve_fit(dec_sin, field_data, R_xx_data, (2000, 5000, 150, -1, 1))
    print(res)
    ax2.plot(field_data, dec_sin(field_data, *res), 'k--')

#    # Plot plateaus
#    vk_const = const.value("von Klitzing constant")
#    for i in range(12, 24+1, 2):
#        ax1.axhline(vk_const/i)
#
#    # Plot theoretical locations
#    density = np.polyfit(field_data[-100:], R_xy_data[-100:], 1)
#    print(f"Fitting to field: {field_data[-100]}.")
#    density = 1/(const.e * density[0])
#    print(f"Extracted Density: {density*1e-4:e} cm^2/Vs")
#    filling_factors = (2*const.pi*density*const.hbar)/(const.e)
#    for i in range(12, 24+1, 2):
#        ax1.axvline(filling_factors/i)

    fig.legend()
    ax1.set_xlabel(r"1/Field $\left(\si{\per\tesla}\right)$")
    ax1.set_ylabel(r"$\textrm{R}_{xy}$ $\left(\si{\ohm}\right)$")
    ax2.set_ylabel(r"$\textrm{R}_{xx}$ $\left(\si{\ohm}\right)$")
    fig.tight_layout()
    fig.show()

def plot_landau_fan():
    dataset = 1035

    dataset = load_by_id(dataset)
    param_data = dataset.get_parameter_data("R_xx", "R_xy")
    voltage_data = param_data["R_xx"]["yoko_voltage"]
    field_data = param_data["R_xx"]["mag_GRPZ_field"]
    R_xx_data = param_data["R_xx"]["R_xx"]
    R_xy_data = param_data["R_xy"]["R_xy"]

    # Reshape data
    lvoltage_data, voltage_data = np.array(detect_cycle(voltage_data))
    field_data = np.array(field_data[::lvoltage_data])
    y, x = np.mgrid[slice(voltage_data[0], voltage_data[-1], complex(0, lvoltage_data+1)),
                    slice(field_data[0], field_data[-1], complex(0, field_data.size))]
    R_xx_data = R_xx_data.reshape((field_data.size, lvoltage_data)).T
    R_xy_data = R_xy_data.reshape((field_data.size, lvoltage_data)).T
    vk_const = const.value("von Klitzing constant")

    fig, ax1 = plt.subplots()
    im = ax1.pcolormesh(y, x, np.log10(R_xx_data), vmin=1, vmax=5, cmap="inferno", rasterized=True)
    fig.colorbar(im)
    ax1.set_xlabel(r"Gate Voltage $\left(\si{\volt}\right)$")
    ax1.set_ylabel(r"$\textrm{B}_{\perp}$ $\left(\si{\tesla}\right)$")
    fig.tight_layout()
    fig.show()

    fig, ax1 = plt.subplots()
    ax1.pcolormesh(y, x, vk_const/R_xy_data, vmin=2, vmax=32, cmap="cividis_r")

    ax1.set_xlabel(r"Gate Voltage $\left(\si{\volt}\right)$")
    ax1.set_ylabel(r"$\textrm{B}_{\perp}$ $\left(\si{\tesla}\right)$")
    fig.tight_layout()
    fig.show()

def density_mobility_scattering_plot(electron_density, mobility, rho_xx=None, width=50e-6, title=""):
    fig, ax1 = plt.subplots()

    # Plot mobility
    ax1.plot(np.abs(electron_density), np.abs(mobility), 'ro', alpha=1, label='{} um device'.format(1e6*width))
    ax1.set_title(title)
    ax1.set_ylabel(r"Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")#, color='r')
    ax1.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    #ax1.tick_params('y', colors='r')
    ax1.set_xlim(0.3e12, 2.5e+12)
    ax1.set_ylim([0, 3.7e+4])

    # Plot upper and lower scattering limit, in the low density limit
    valid_range = np.where(np.abs(electron_density) < 0.65e12)
    new_ed = np.abs(electron_density)[valid_range]
    new_mob = np.abs(mobility)[valid_range]
    def power_law(density, A, alpha):
        return A * (density**alpha)
    A1 = curve_fit(partial(power_law, alpha=1.5), new_ed, new_mob, (1e-2,))[0][0]
    A2 = curve_fit(partial(power_law, alpha=0.5), new_ed, new_mob, (1e-2,))[0][0]
    A3, alpha = curve_fit(power_law, new_ed, new_mob, (1e-2, 1.5))[0]
    ax1.plot(np.abs(electron_density), power_law(np.abs(electron_density), A1, 1.5), 'k--')
    ax1.plot(np.abs(electron_density), power_law(np.abs(electron_density), A2, 0.5), 'k--')
    ax1.plot(np.abs(electron_density), power_law(np.abs(electron_density), A3, alpha), 'b-')
    print(A3, alpha)

    fig.tight_layout()
    return fig


def plot_supp_mob_den(analyzed_data, position="any"):
    """
    This function plots the density at zero voltage vs peak density achieved
    """
    exclude_ids = ()
    n = {(None, 'H2O'): 'No Pretreatment / H2O Oxidizer',
         ('tma', 'H2O'): 'TMA Reduction / H2O Oxidizer',
         ('tma', 'O3'): 'TMA Reduction / O3 Oxidizer',
         ('plasma', 'H2O'): 'H2 Plasma / H2O Oxidizer',
         ('plasma', 'O3'): 'H2 Plasma / O3 Oxidizer'}
    # Extract a list of peak mobilities and voltages
    gate_voltages = defaultdict(list)
    densities = defaultdict(list)
    mobilities = defaultdict(list)
    for sid, sample in analyzed_data.items():
        if sid in exclude_ids:
            continue
        if position != "any":
            if sample.position != position:
                continue
        # Get the measured values from the measurements for this sample
        for measurement in sample.measurements:
            if np.isclose(measurement.width, 25e-6): # Only use 50um samples
                continue
            t = (sample.treatment, sample.ald["oxidiser"])
            print(t, measurement.vg_at_peak, measurement.den_at_peak)
            gate_voltages[t].append(measurement.analyzed.vg)
            densities[t].append(measurement.analyzed.density)
            mobilities[t].append(measurement.analyzed.mobility)

    fig = plt.figure(constrained_layout=True)
    spec = gridspec.GridSpec(ncols=4, nrows=3, figure=fig, left=0.05, right=0.95, wspace=0.05)
    axes = []
    axes.append(fig.add_subplot(spec[0, 1:3]))
    axes.append(fig.add_subplot(spec[1, 0:2]))
    axes.append(fig.add_subplot(spec[1, 2:4]))
    axes.append(fig.add_subplot(spec[2, 0:2]))
    axes.append(fig.add_subplot(spec[2, 2:4]))

    for i, (treatment, name) in enumerate(n.items()):
        ax1 = axes[i]
        ax2 = ax1.twinx()
        #for vg, density, mobility in zip(gate_voltages[treatment], densities[treatment], mobilities[treatment]):
        if gate_voltages[treatment]:
            vg, density, mobility = gate_voltages[treatment][-1], densities[treatment][-1], mobilities[treatment][-1]
            # Plot electron density
            ax1.plot(vg, np.abs(density), 'mx')
            ax1.set_ylabel(r"Density $\left(\si{\per\square\centi\meter}\right)$", color='m')
            ax1.set_xlabel('applied top gate voltage (V)')
            ax1.tick_params('y', colors='m')
            ax1.set_ylim(0, 3.0e+12)

            # Filter mobility by density
            filt = np.where(np.abs(density) > 0.2e12)

            ax2.plot(vg[filt], np.abs(mobility[filt]), 'co')
            ax2.set_ylabel(r"Mobility $\left(\si{\square\centi\meter\per\volt\per\second}\right)$", color='c')
            ax2.tick_params('y', colors='c')
            ax2.set_ylim(0, 4.7e+4)

        ax1.set_xlim(-0.7, 0.5)
        ax1.set_title(name)
    fig.set_size_inches(12, 12/np.sqrt(2))
    #plt.tight_layout()

def plot_short_resistance():
    all_data = open_data_sequence("2018-01-29", 2, 4)

    points = (0, 0.5, 1.0, 1.5, 2.0)

    field = "ami_field_set"
    lockin = "SR860_2"
    param = f"{lockin}_R"

    fig, ax = plt.subplots()
    all_dp = []
    for d in all_data:
        exc_amp = d.snapshot()["station"]["instruments"][lockin]["parameters"]["amplitude"]["value"]

        input_imp = 100
        output_imp = 50
        res = (exc_amp / getattr(d, param).ndarray) - input_imp - output_imp

        zerofield = np.where(np.isclose(getattr(d, field), 0, atol=2.5e-3))
        zerofield_res = np.mean(res[zerofield])

        norm_res = res / zerofield_res
        delta_res = res - zerofield_res

        field_data = getattr(d, field)[1:-1]
        norm_res = norm_res[1:-1]
        delta_res = delta_res[1:-1]

        dp = []
        for p in points:
            loc = np.where(np.isclose(np.abs(field_data), p, atol=5e-3))
            mean_offs = np.mean(delta_res[loc])
            dp.append((p, mean_offs))
        all_dp.extend(dp)

        #ax.plot(field_data, norm_res)
        ax.scatter(*zip(*dp))

    #fit = np.polyfit(*zip(*all_dp), 1)
    #ax.plot(points, np.poly1d(fit)(points), 'k--')

    ax.set_xlabel("Field (T)")
    ax.set_ylabel("Resistance Change (Ω)")
    ax.set_ylim(-2, 20)
    fig.tight_layout()

def plot_overlaid_sheet_res(analyzed_data):
    resistances = {}

    for dev in analyzed_data:
        if not analyzed_data[dev].measurements:
            continue

        best_res = None
        best_vgs = None
        best_mob = 0

        for measurement in analyzed_data[dev].measurements:
            if measurement.mob_at_peak > best_mob:
                best_res = measurement.analyzed.res
                best_vgs = measurement.analyzed.vg
                best_mob = measurement.mob_at_peak

        resistances[dev] = (best_mob, best_res, best_vgs)

    fig, ax = plt.subplots()

    for dev in resistances:
        ax.semilogy(resistances[dev][2], resistances[dev][1], label=dev)
    ax.legend()
    ax.set_ylim(100, 10000)
    ax.set_xlabel("V_tg (V)")
    ax.set_ylabel("Sheet Resistance (Ohm)")
    fig.tight_layout()
