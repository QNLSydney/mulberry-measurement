# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 18:56:19 2018

@author: LD2007
"""

from time import sleep, clock
from collections import Iterable
import requests
import qcodes as qc
from qcodes import Parameter, MultiParameter
import math
import os, re
import warnings

class TimeParam(Parameter):
    def __init__(self, waittime):
        self._waittime = waittime
        self._start = clock()
        super().__init__("time",
             label="Time",
             unit="s",
             snapshot_get=False)
    
    def get_raw(self):
        return clock() - self._start
    
    def set_raw(self, val):
        sleep(self._waittime)
        
class FridgeTemps(MultiParameter):
    def __init__(self, fridge, url):
        self.url = url
        
        params = requests.get(url)
        if params.status_code != 200:
            raise RuntimeError("Unable to query fridge")
        params = set(params.json().keys())
        params.remove("Time")
        params = tuple(params)
        self.params = params
        
        super().__init__(
                "{}_temps".format(fridge),
                names=params,
                shapes=tuple(() for _ in params),
                units=tuple("K" for _ in params),
                snapshot_get=False)
        
    def get_raw(self):
        temps = requests.get(self.url)
        if temps.status_code != 200:
            raise RuntimeError("Unable to query fridge")
        temps = temps.json()
        temps = [temps[therm] for therm in self.params]
        return tuple(temps)

class SourceDrainVoltages(Parameter):
    def __init__(self, params):
        self.params = params
        super().__init__("v_sd",
             label="V_sd",
             unit="V")
    
    def get_raw(self):
        return self.params[0].get()
    
    def set_raw(self, value):
        for param in self.params:
            param.set(value)

class LockinResistance(Parameter):
    def __init__(self, lockin, input_imp=None, current_scale=1, voltage_scale=1, **kwargs):
        self.lockin = lockin
        self.input_imp = input_imp
        self.current_scale = current_scale
        self.voltage_scale = voltage_scale
        
        name = "{}_resistance".format(lockin.name)
        
        super().__init__(name, 
             label="Ohms",
             unit="Ohms",
             snapshot_get=False,
             get_cmd=self._get_resistance)
        
    def _get_resistance(self):
        # Figure out instrument impedances
        output_imp = 50
        if self.input_imp is None:
            input_imp = 1000 if self.lockin.input_gain() == 100e6 else 100
        else:
            input_imp = self.input_imp
        volt = self.lockin.amplitude()/self.voltage_scale
        curr = self.lockin.R()/self.current_scale
        res = volt/curr
        res = res - output_imp - input_imp
        return res

class LockinResistances(MultiParameter):
    def __init__(self, params):
        self.params = params
        super().__init__(name="Lockin_Resistances",
            names=tuple(param.name for param in params),
            shapes=tuple(() for _ in params),
            units=tuple(param.unit for param in params),
            snapshot_get=False)
    
    def get_raw(self):
        return [param.get() for param in self.params]
    
class MBSwitchedParam(Parameter):
    def __init__(self, md, switch, param, waittime):
        self.param = param
        self.switch = switch
        self.md = md
        self.waittime = waittime
        super().__init__(name="{}_{}".format(switch, param.full_name),
            label=param.label,
            unit=param.unit,
            snapshot_get=False)
    
    def get_raw(self):
        if self.md.select.get_latest() != self.switch:
            self.md.select(self.switch)
            sleep(self.waittime)
        return self.param()
    
class CurrentAmplifier(Parameter):
    def __init__(self, param, gain):
        self.param = param
        self.gain = gain
        super().__init__(name="{}_ithaco".format(param.full_name),
              label="Current",
              unit="A",
              scale=gain,
              snapshot_get=False)
        
    def get_raw(self):
        return self.param.get()

def printer(data):
    def p():
        index = data.Still.last_saved_index
        temp = data.Still[index]
        print("@Temp: {}".format(temp))
        for switch in ("B", "C", "D", "E"):
            print("\tSW: {}".format(switch), end="")
            for lockin in range(1, 4):
                d = getattr(data, "{}_SR860_{}_resistance".format(switch, lockin))[index]
                print("\t{}: {:.1f}Ω".format(lockin, d), end="")
            print()
    return p

def cooldown_loop(t, ft, resist, loops):
    loop = qc.Loop(t[0:loops:1])
    loop = loop.each(t, ft, *resist)
    data = loop.get_data_set()
    p = printer(data)
    loop.with_bg_task(p)
    loop.run()
    return data

def field_sweep(ami, voltages, start, stop, points):
    loop = qc.Loop(ami.field.sweep(start, stop, num=points))
    if isinstance(voltages, Iterable):
        loop = loop.each(*voltages)
    else:
        loop = loop.each(voltages)
    data = loop.get_data_set()
    loop.run()
    return data

def do_field_sweeps(ami, md, pairs, params):
    for pair in pairs:
        md.select(pair)
        sleep(10)
        delay = TimeParam(1)
        field_sweep(ami, params + [delay], ami.field(), -ami.field(), 2000)
    md.clear()

def find_data(date, num):
    date = str(date)
    data_dir = os.path.join("data", date)
    files = os.listdir(data_dir)
    for file in files:
        if re.match("#{:03}".format(num), file):
            return os.path.join(data_dir, file)

# Params:
#   md: mulberry driver
#   tg: top_gate voltage param
#   V_sd: source-drain voltage param
#   params: list of parameters to measure
#   bg: (start, stop) tuple of topgate voltages
#   iv: (start, stop) tuple of source-drain voltages
#   plot_params: list of strings of parameters to live plot
def do_bg_iv(md, tg_set, V_sd, states, params, tg, sd, 
             plot_params=[], bg_step=0.02, iv_step=0.00004,
             wait_time=5):
    for state in states:
        md.select(state)
        tg_set.set(tg[0])
        V_sd.set(sd[0])
        sleep(wait_time)
        
        # Inner Loop (V_sd)
        inner_loop = qc.Loop(V_sd.sweep(sd[0], sd[1], step=iv_step))
        inner_loop = inner_loop.each(*params)
        inner_loop = inner_loop.then(qc.Task(V_sd.set, sd[0]), qc.Wait(wait_time))
        
        # Outer Loop (tg)
        loop = qc.Loop(tg_set.sweep(tg[0], tg[1], step=bg_step))
        loop = loop.each(inner_loop)
        data = loop.get_data_set()
        
        # Show each plot
        plots = []
        for param in plot_params:
            param_data = getattr(data, param, None)
            if param_data is None:
                warnings.warn("Missing parameter {} not plotted".format(param),
                              RuntimeWarning)
                continue
            plot = qc.QtPlot()
            plot.add(param_data, name=param, title=param)
            plots.append(plot)
        # Define function to live update all plots
        def upd():
            for plot in plots:
                plot.update()
                
        # Update plots at the end of each line
        loop = loop.with_bg_task(upd)
        
        # Run the loop
        loop.run()
    
    # Set parameters to zero
    tg_set.set(0)
    V_sd.set(0)
    
    # Clear mulberry at end
    md.clear()
        
