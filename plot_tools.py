from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import scipy.constants as const

# Initialize matplotlib
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r"\usepackage{siunitx}")
plt.ion()

def density_mobility_plot(electron_density, mobility, rho_xx=None, width=50e-6, title=""):
    fig, ax1 = plt.subplots()

    # Plot mobility
    ax1.plot(np.abs(electron_density), np.abs(mobility), 'rx', label='{} um device'.format(1e6*width))
    ax1.set_title(title)
    ax1.set_ylabel(r"Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$", color='r')
    ax1.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax1.tick_params('y', colors='r')
    ax1.set_xlim(0, 3.6e+12)
    ax1.set_ylim([0.0, 4.6e+4])

    # Plot rho_xx, if given
    if rho_xx is not None:
        ax2 = ax1.twinx()
        ax2.plot(np.abs(electron_density), np.abs(rho_xx), 'bx')
        ax2.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
        ax2.set_ylabel(r'rho_xx $\left(\textOmega \Box^{-1}$\right)', color='b')
        ax2.tick_params('y', colors='b')
        ax2.set_ylim([0, 4000])
    fig.tight_layout()
    return fig

def mobility_voltage_plot(Vg, mobility, title=""):
    fig, ax = plt.subplots()

    # Plot mobility
    ax.plot(Vg, mobility, 'mx')
    ax.set_ylabel(r"Mobility $\left(\frac{cm^2}{V \cdot s}\right)$", color='c')
    ax.set_xlabel('applied top gate voltage (V)')
    ax.tick_params('y', colors='c')
    ax.set_ylim(0, 4.5e4)
    ax.set_title(title)
    fig.tight_layout()
    return fig

def density_mobility_voltage_plot(Vg, electron_density, mobility=None, width=50e-6, title=""):
    fig, ax1 = plt.subplots()
    # Plot electron density
    ax1.plot(Vg, np.abs(electron_density), 'mx')
    ax1.set_xlabel(r"Density $\left(cm^{-2}\right)$", color='m')
    ax1.set_xlabel('applied top gate voltage (V)')
    ax1.tick_params('y', colors='m')
    #ax1.set_xlim(-0.9, 1.0)
    ax1.set_ylim(0, 3.6e+12)

    # Plot mobility
    if mobility is not None:
        ax2 = ax1.twinx()
        ax2.plot(Vg, np.abs(mobility), 'cx', label='{} um device'.format(1e6*width))
        ax2.set_ylabel(r"Mobility $\left(\frac{cm^2}{V \cdot s}\right)$", color='c')
        ax2.tick_params('y', colors='c')
        ax2.set_ylim(0, 4.5e+4)

    ax1.set_title(title)
    fig.tight_layout()
    return fig

def sheet_resistance_plot(Vg, rho_xx, width=50e-6, title=""):
    fig, ax = plt.subplots()
    ax.plot(Vg, np.abs(rho_xx), 'bx', label='{} um device'.format(1e6*width))
    ax.set_title(title)
    ax.set_ylabel(r'Sheet Resistance $\left(\textOmega \Box^{-1}$\right)', color='b')
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
    for loc in c:
        x_vals = [ord(d[1]) - ord('A') for d in data if d[2] == loc]
        y_vals = [d[6] for d in data if d[2] == loc]

        ax.scatter(x_vals, y_vals, c=c[loc], alpha=0.6, label=loc)

    x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/H2O", "ArH Plasma\nTMA/O3")

    ax.set_xticks(np.arange(5))
    ax.set_xticklabels(x_ticks)
    ax.set_ylabel(r"Peak Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_title("Peak Mobility v. Treatment")
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
    for treatment in c:
        x_vals = [d[8] for d in data if d[1] == treatment]
        y_vals = [d[6] for d in data if d[1] == treatment]

        ax.scatter(x_vals, y_vals, c=c[treatment], alpha=0.6, label=treat[treatment])

    ax.set_xlabel(r"Top Gate Voltage $\left(\si{\volt}\right)$")
    ax.set_ylabel(r"Peak Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_title("Gate Voltage v. Peak Mobility")
    ax.legend()
    fig.tight_layout()
    fig.savefig("Peak Mobility Voltage.png")
    plt.show()

def plot_den_mob_trace(analyzed_data):
    fig, ax = plt.subplots()

    c = {'A': 'tab:red', 'B': 'tab:orange', 'C': 'tab:olive', 'D': 'tab:blue', 'E': 'tab:green'}
    x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/O3", "ArH Plasma\nTMA/H2O")
    treat = {a: b for a, b in zip(c.keys(), x_ticks)}
    for device in c:
        density = analyzed_data[device].density
        mobility = analyzed_data[device].mobility
        current = analyzed_data[device].curr

        valid_range = np.where(np.all((current > 3.25e-9, density > 0.2e12, density < 3e12, mobility > 5000), axis=0))
        density = density[valid_range]
        mobility = mobility[valid_range]
        current = current[valid_range]

        ax.scatter(density, mobility, c=c[device], alpha=0.6, label=treat[device])

    ax.set_xlabel(r"Density $\left(\si{\per\square\centi\meter}\right)$")
    ax.set_ylabel(r"Mobility $\left(\si[per-mode=fraction]{\square\centi\meter\per\volt\per\second}\right)$")
    ax.set_title("Mobility Density")
    ax.legend()
    fig.tight_layout()
    fig.savefig("Mobility Density.png")
    plt.show()

def plot_den_vg_trace(analyzed_data):
    fig, ax = plt.subplots()

    c = {'A': 'tab:red', 'B': 'tab:orange', 'C': 'tab:olive', 'D': 'tab:blue', 'E': 'tab:green'}
    x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/O3", "ArH Plasma\nTMA/H2O")
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
    x_ticks = ("No Pretreatment\nTMA/H2O", "TMA Prepulse\nTMA/H2O", "TMA Prepulse\nTMA/O3", "ArH Plasma\nTMA/O3", "ArH Plasma\nTMA/H2O")
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
    ax.set_title("Density/Mobility at Zero Vg to Peak Mobility")
    ax.legend(scatters, x_ticks, ncol=2)
    fig.tight_layout()
    fig.savefig("Density Mobility Zero.png")
    plt.show()
