# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 18:56:19 2018

@author: LD2007
"""

from time import sleep, clock
import requests
import qcodes as qc
from qcodes import Parameter, MultiParameter

class TimeParam(Parameter):
    def __init__(self, waittime):
        self._waittime = waittime
        self._start = clock()
        super().__init__("time",
             label="Time",
             unit="s",
             snapshot_get=False,
             get_cmd=self._get,
             set_cmd=self._set)
    
    def _get(self):
        return clock() - self._start
    
    def _set(self, val):
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
        
    def get(self):
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
    
    def get(self):
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
    
    def get(self):
        items = []
        for pair in self.pairs:
            self.md.select(pair)
            sleep(self.waittime)
            items.append(self.param.get())
        self.md.clear()
        return tuple(items)
        

def gen_resistances_param(lockins):
    params = [LockinResistance(lockin) for lockin in lockins]
    
    return LockinResistances(params)

def meas_all(md, pairs, param):
    lsr = LockinSwitchResistances(md, pairs, param, 10)
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