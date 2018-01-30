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
import os, re

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
        

class LockinResistance(Parameter):
    def __init__(self, lockin, **kwargs):
        self.lockin = lockin
        
        name = "{}_resistance".format(lockin.name)
        
        super().__init__(name, 
             label="Ohms",
             unit="Ohms",
             snapshot_get=False,
             get_cmd=self._get_resistance)
        
    def _get_resistance(self):
        # Figure out instrument impedances
        output_imp = 50
        input_imp = 1000 if self.lockin.input_gain() == 100e6 else 100
        res = self.lockin.amplitude()/self.lockin.R()
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
    
class LockinSwitchResistances(MultiParameter):
    def __init__(self, md, pairs, param, waittime):
        self.param = param
        self.pairs = pairs
        self.md = md
        self.waittime = waittime
        super().__init__(name="Lockin_Switch_Resistances",
            names=tuple("{}_{}".format(pair, param.name) for pair in pairs),
            shapes=tuple((len(param.shapes),) for _ in pairs),
            units=tuple(param.units[0] for _ in pairs),
            snapshot_get=False)
    
    def get_raw(self):
        items = []
        for pair in self.pairs:
            self.md.select(pair)
            sleep(self.waittime)
            items.append(self.param.get())
        self.md.clear()
        return tuple(items)

class LockinsComb(MultiParameter):
    def __init__(self, params):
        self.params = params
        self.param_name = params[0].name
        super().__init__(name="Lockins_{}".format(self.param_name),
              names=tuple("{}".format(param.full_name) for param in self.params),
              shapes=tuple(() for _ in self.params),
              units=tuple(param.unit for param in self.params),
              snapshot_get=False)
    
    def get_raw(self):
        return [param.get() for param in self.params]

class MBLockins(MultiParameter):
    def __init__(self, md, pairs, params, waittime):
        self.pairs = pairs
        self.md = md
        self.waittime = waittime
        if not isinstance(params, Iterable):
            params = (params,)
        self.params = params
        
        shapes = []
        names = []
        units = []
        for pair in pairs:
            for param in params:
                names.append("{}_{}".format(pair, param.full_name))
                if isinstance(param, MultiParameter):
                    shapes.append((len(param.shapes),))
                    units.append(param.units)
                else:
                    shapes.append(())
                    units.append(param.unit)
        shapes = tuple(shapes)
        names = tuple(names)
        units = tuple(units)
        
        super().__init__(name="MB_Switched",
            names=names,
            shapes=shapes,
            units=units,
            snapshot_get=False)
    
    def get_raw(self):
        items = []
        i = 0
        for pair in self.pairs:
            self.md.select(pair)
            sleep(self.waittime)
            for param in self.params:
                items.append(param.get())
                i += 1
        self.md.clear()
        return tuple(items)
    
def gen_resistances_param(lockins):
    params = [LockinResistance(lockin) for lockin in lockins]
    
    return LockinResistances(params)

def meas_all(md, pairs, param):
    lsr = MBLockins(md, pairs, param, 10)
    m = qc.Measure(lsr)
    data = m.get_data_set()
    m.run()
    return data

def printer(data):
    def p():
        index = data.Still.last_saved_index
        temp = data.Still[index]
        print("@Temp: {}".format(temp))
        for switch in ("B", "C", "D", "E"):
            d = data.__getattr__("{}_Lockin_Resistances".format(switch))[index]
            print("\tSW {}: Thru: {:.2e} Ohms, 25um: {:.2e} Ohms, 50um: {:.2e} Ohms".format(switch, *d))
    return p

def cooldown_loop(t, ft, resist, loops):
    loop = qc.Loop(t[0:loops:1])
    loop = loop.each(t, ft, resist)
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
        delay = TimeParam(1)
        field_sweep(ami, params + [delay], ami.field(), -ami.field(), 1000)
    md.clear()