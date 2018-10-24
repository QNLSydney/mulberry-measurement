# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 18:56:19 2018

@author: LD2007
"""

from time import sleep, clock
from collections import Iterable
import requests
import qcodes as qc
from qcodes import Parameter, MultiParameter, new_data_set
import math
import os, re
import warnings
import numpy as np

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
        curr = self.lockin.X()/self.current_scale
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
    
class VoltageAmplifier(Parameter):
    def __init__(self, param, gain):
        self.param = param
        self.gain = gain
        super().__init__(name="{}_preamp".format(param.full_name),
              label="Voltage",
              unit="V",
              scale=gain,
              snapshot_get=False)
        
    def get_raw(self):
        return self.param.get()

def printer(data):
    def p():
        index = data.Still.last_saved_index
        temp = data.Still[index]
        print("@Temp: {}".format(temp))
        for switch in ("A", "B", "C", "D", "E"):
            print("\tSW: {}".format(switch), end="")
            for lockin in range(1, 3):
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

def field_sweep(ami, parameters, start, stop, points):
    loop = qc.Loop(ami.field.sweep(start, stop, num=points))
    if isinstance(parameters, Iterable):
        loop = loop.each(*parameters)
    else:
        loop = loop.each(parameters)
    data = loop.get_data_set()
    loop.run()
    return data

def do_field_sweeps(ami, md, pairs, params, start=-2, stop=2, points=1000):
    for i, pair in enumerate(pairs):
        md.select(pair)
        sleep(10)
        delay = TimeParam(1)
        field_sweep(ami, params + [delay], start, stop, points)
        start, stop = stop, start
    md.clear()
    ami.field.set(0)
    
def do_field_sweeps_at_V(ami, md, pairs, params, V, start=-2, stop=2, points=1000):
    
    yoko_t.voltage.post_delay =  0.1
    yoko_t.voltage.step = 0.01
    yoko_t.voltage.inter_delay = 0.1
    yoko_t.voltage(0)
    yoko_t.output('on')
    yoko_t.voltage(V)
    
    for i, pair in enumerate(pairs):
        md.select(pair)
        sleep(10)
        delay = TimeParam(1)
        field_sweep(ami, params + [delay], start, stop, points)
        start, stop = stop, start
    md.clear()
    ami.field.set(0)
    yoko_t.voltage(0)
    yoko_t.output('off')
    


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

"""
    Top gate measurement functions
"""

def continuity():
    global md, resistances
    for switch in ("A", "B", "C", "D", "E"):
        md.select(switch)
        sleep(5)
        print(list(resistance() for resistance in resistances))
    md.clear()
    
    
    
        
def top_gate_leak_test(pairs, parameters, start=-0.8, stop=1.0):  # value in volts
    global md, ami, yoko_t, lockins
    E_field_start = start      
    E_field_stop = stop
    
    yoko_t.voltage.post_delay =  0.1
    yoko_t.voltage.step = 0.01
    yoko_t.voltage.inter_delay = 0.1
    yoko_t.voltage(0)
    yoko_t.output('on')
    
    for pair in pairs:
        md.select(pair) 
        sleep(10)
        yoko_t.voltage(E_field_start)
        yoko_t.voltage.post_delay =  7*lockins[0].time_constant()
        loop = qc.Loop(yoko_t.voltage.sweep(E_field_start,E_field_stop,num=50))
        if isinstance(parameters, Iterable):
            loop = loop.each(*parameters)
        else:
            loop = loop.each(parameters)
        loop.run(name="Leak_test_MD_{}_".format(pair))
        yoko_t.voltage.post_delay =  0.1
        yoko_t.voltage(0)
    md.clear()
    yoko_t.output('off') 
    ami.field.set(0)
    
    
        
            
def top_gate_step_B_field_sweep(pairs, parameters, start=-0.015, stop=0.015, points=301):
    global md, ami, yoko_t, lockins
    test_values = np.linspace(-0.6,-0.2,30)
    yoko_t.voltage(0)
    yoko_t.output('on')
    yoko_t.voltage.step = 0.01
    yoko_t.voltage.inter_delay = 0.1
    tau = 7*lockins[0].time_constant()
    for pair in pairs:
        md.select(pair)
        for item in test_values:        
            yoko_t.voltage(item)
            ami.field.set(start)
            sleep(600)
            loop = qc.Loop(ami.field.sweep(start, stop, num=points), delay=tau)
            if isinstance(parameters, Iterable):
                loop = loop.each(*parameters)
            else:
                loop = loop.each(parameters)
            loop.run(name="field_sweep_MD_{}_".format(pair))
        yoko_t.voltage(0)
    yoko_t.output('off')
    yoko_t.voltage(0) 
    ami.field.set(0)
    
    
    
    
def top_gate_sweep_B_field_step(md, pairs, parameters):
    global ami, yoko_t, lockins
    B_field_values =  [-1.0,-0.5,-0.25,-0.05,0.25,0.5,1.0]  # if you change these you need to change the values in the analysis
    E_field_start = -0.7     
    E_field_stop = 0.5
    yoko_t.voltage.post_delay =  0.1
    yoko_t.voltage.step = 0.01
    yoko_t.voltage.inter_delay = 0.1
    yoko_t.voltage(0)
    yoko_t.output('on')
    for field in B_field_values:
        ami.field.set(field)
        for pair in pairs:
            md.select(pair)
            sleep(10)
            yoko_t.voltage(E_field_start)
            yoko_t.voltage.post_delay =  5*lockins[0].time_constant()
            loop = qc.Loop(yoko_t.voltage.sweep(E_field_start,E_field_stop,num=200))   # same as above, here
            if isinstance(parameters, Iterable):
                loop = loop.each(*parameters)
            else:
                loop = loop.each(parameters)
            loop.run(name="MD_{}_FIELD_{}".format(pair, field))
            yoko_t.voltage.post_delay =  0.1
            yoko_t.voltage(0)
    md.clear()
    yoko_t.output('off') 
    ami.field.set(0)
    
    
  
    
def inverse_B_field_sweep(md, pairs, parameters):
    global ami, yoko_t, lockins
    
    start_field = 1.2
    stop_field = 2.0   
    step_num = 200
    steps = np.linspace(1/start_field,1/stop_field,step_num)
    field_vals = np.array([1/steps])
    
    tg_values = [-0.4,-0.2,-0.15,-0.1,-0.05,0] # get n_e at these exact biases - fit?
    yoko_t.voltage.step = 0.01
    yoko_t.voltage.inter_delay = 0.1
    yoko_t.voltage(0)
    yoko_t.output('on')
            
    for pair in pairs:
        md.select(pair)
        tau = 7*lockins[0].time_constant()
        
        for item in tg_values:            
            
            yoko_t.voltage(item)
            ami.field.set(start_field)
            sleep(300)
            
            loop = qc.Loop(ami.field[field_vals], delay=tau)
            if isinstance(parameters, Iterable):
                loop = loop.each(*parameters)
            else:
                loop = loop.each(parameters)
            loop.run()
            
        yoko_t.voltage(0)
        
    md.clear()
    yoko_t.output('off') 
    ami.field.set(0)


    
    
def lockin_IV(pairs, parameters, start=-100, stop=100):
    global md, ami, yoko_t, lockins
    
    resistance_value = 1.0e+7                           # eries resistor Ohms
    excitation = 1.0                                    # nA
    points = int((np.abs(start)+np.abs(stop)+1)/excitation)
    I_range = np.linspace(start,stop,points)               # nA
    print(points)
    
    lockins[0].amplitude.set(resistance_value*excitation*1.e-9)
    tau = 7*lockins[0].time_constant()
    DC_out = (resistance_value*I_range*1.e-9)                       
    
    for pair in pairs:
        md.select(pair)
        sleep(10)
                        
        loop = qc.Loop(lockins[0].sine_outdc[DC_out], delay=tau)
        if isinstance(parameters, Iterable):
            loop = loop.each(*parameters)
        else:
            loop = loop.each(parameters)
        loop.run(name="IV_MD_{}_".format(pair))          
    md.clear()
    
    
    
def top_gate_step_B_field_sweep_var(pairs, parameters, le_array, E_start, E_stop, E_steps):
    '''with variable field limits dependent on le'''
    global md, ami, yoko_t, lockins
    test_values = np.linspace(E_start,E_stop,E_steps)
    yoko_t.voltage.step = 0.01
    yoko_t.voltage.inter_delay = 0.1
    yoko_t.voltage(0)
    yoko_t.output('on')

    tau = 7*lockins[0].time_constant()
    for pair in pairs:
        md.select(pair)
        for item in test_values:        
            yoko_t.voltage(item)
            
            h_bar = 1.0545718e-34
            e = 1.60217662e-19
            
            poly = np.poly1d(le_array)
            l_e = poly(item)
            B_e = np.abs(h_bar/(4*e*l_e**2))     
            print(B_e)
            if B_e > 0.010:
                B_e = 0.010
            start = B_e + 0.005
            stop = -1*start
            
            points = int(round((start - stop)/0.00005))
            print(points)
            
            ami.field.set(start)
            sleep(60)
            loop = qc.Loop(ami.field.sweep(start, stop, num=points), delay=tau)
            if isinstance(parameters, Iterable):
                loop = loop.each(*parameters)
            else:
                loop = loop.each(parameters)
            loop.run(name="field_sweep_MD_{}_".format(pair))
        yoko_t.voltage(0)
    yoko_t.output('off')
    yoko_t.voltage(0) 
    ami.field.set(0)    
    
    
def top_gate_step_B_field_sweep_var(pairs, parameters, le_array, E_start, E_stop, E_steps):
    '''with variable field limits dependent on le
    WITH YOKO MAGNET POWER SUPPLY'''    
    global md, yoko_mag, yoko_t, lockins
    test_values = np.linspace(E_start,E_stop,E_steps)
    yoko_t.voltage.step = 0.01
    yoko_t.voltage.inter_delay = 0.1
    yoko_t.voltage(0)
    yoko_t.output('on')

    tau = 7*lockins[0].time_constant()
    for pair in pairs:
        md.select(pair)
        
        yoko_t.voltage(test_values[0])
        sleep(60*5)
        
        for item in test_values:        
            
            h_bar = 1.0545718e-34
            e = 1.60217662e-19
            
            poly = np.poly1d(le_array)
            l_e = poly(item)
            print(l_e)
            B_e = np.abs(h_bar/(4*e*l_e**2))     
            print(B_e)
            if B_e > 0.010:
                B_e = 0.010
            if B_e < 0.001:
                B_e = 0.001
            
            step_f = 0.000005
            step_c = 0.0001
            coarse_range = B_e + 0.001
            fine_range = 0.0015
            
            # Round coarse_range to the smaller multiple of step_c
            coarse_range = int(coarse_range/step_c)*step_c

            # Convert from tesla to A
            coarse_range = coarse_range*(1/0.0641)
            fine_range = fine_range*(1/0.0641)
            step_f = step_f*(1/0.0641)
            step_c = step_c*(1/0.0641)
            
            # Create Segments
            seg1 = yoko_mag.current.sweep(coarse_range, fine_range, step=step_c)
            seg2 = yoko_mag.current.sweep(fine_range, -fine_range, step=step_f)
            seg3 = yoko_mag.current.sweep(-fine_range, -coarse_range, step=step_c)
            
            yoko_t.voltage(item)
            yoko_mag.current.set(coarse_range)
            sleep(20)
            
            loop = qc.Loop(seg1+seg2+seg3, delay=tau)
            if isinstance(parameters, Iterable):
                loop = loop.each(*parameters)
            else:
                loop = loop.each(parameters)
            loop.run(name="field_sweep_MD_{}_".format(pair))
        yoko_t.voltage(0)
    yoko_t.output('off')
    yoko_t.voltage(0) 
    yoko_mag.current.set(0)


def combo(): 
    le_array = [3.8550800197375135, 3.44431569778274, -3.2755097452342588, -3.3781185463555525, 0.97781231908574107, 1.2934599991451401, -0.089643803990583998, -0.23432950291634186, -0.011830652354101209, 0.017386704364718043, 0.0029019847014007875, 0.00021904412474488144, -9.3042687289569489e-05, -7.0584820377935181e-05, -2.1039052974546231e-05, -2.7884246144565719e-06, 2.0623194846359774e-06, 5.5010591855292152e-07]
    top_gate_step_B_field_sweep_var(('E'), voltages, le_array, -0.25, 0.5, 50)
    
    le_array = [2.8872147858341384, 1.4844142937134688, -3.5681356758556428, -1.6089839132938182, 1.9261763124050821, 0.70851885882605226, -0.59511327902472388, -0.15933279716703977, 0.11499147895017386, 0.018023640714713227, -0.014041338877701726, -0.00058320985694616453, 0.00099677077156168988, -7.5618243768608198e-05, -2.6217166246436922e-05, 5.6634426623554682e-06, -8.5841828224469601e-07, 2.7414494188420501e-07]
    top_gate_step_B_field_sweep_var(('C'), voltages, le_array, -0.5, 0.25, 50)
    
    le_array = [1.8011214225584711, 0.34157390962799689, -2.8829700363255677, -0.41284554900496695, 2.0071658424581837, 0.23348581280264394, -0.75786221699657996, -0.067503576915902919, 0.16638449147669321, 0.0085250857417301928, -0.021492254894203089, 2.7030333104555563e-05, 0.0015345987456724555, -0.00010842897836591451, -4.2847045012635924e-05, 6.8398895644934003e-06, -8.0950405369759159e-07, 2.7428020240829887e-07]
    top_gate_step_B_field_sweep_var(('A'), voltages, le_array, -0.5, 0.25, 50)
        
    le_array = [-0.24837429751492376, 0.74662960133443812, 0.75388756574432936, -0.86906350499506579, -0.59934361330575547, 0.44468130852990767, 0.21384300348526536, -0.12952482337124993, -0.037544000848299335, 0.022710586854185307, 0.0026215841557996229, -0.0022416455280034318, 7.8426620036530051e-05, 8.8561524754453101e-05, -1.4123321118593649e-05, 9.7357812063617862e-07, -4.2477725456855097e-07, 2.0958445297615039e-07]
    top_gate_step_B_field_sweep_var(('B'), voltages, le_array, -0.5, 0.25, 50)
    
    le_array = [-0.17855665406658161, -0.28068666598167746, 0.079690265024948456, 0.2910348391599239, 0.020462818897240304, -0.12683859091491673, -0.018028380664479506, 0.029932096888560911, 0.0035555259484676212, -0.0039263623848528493, -0.00019247919509747386, 0.00023010328643151019, -1.1719566440806479e-07, -1.2488933199170218e-06, -1.275617054419993e-06, 4.0488132826110289e-07, -3.0834032428491116e-07, 1.4146975152802756e-07]
    top_gate_step_B_field_sweep_var(('D'), voltages, le_array, -0.75, 0.0, 50)
    
def combo_p1(): 
    le_array = [3.8550800197375135, 3.44431569778274, -3.2755097452342588, -3.3781185463555525, 0.97781231908574107, 1.2934599991451401, -0.089643803990583998, -0.23432950291634186, -0.011830652354101209, 0.017386704364718043, 0.0029019847014007875, 0.00021904412474488144, -9.3042687289569489e-05, -7.0584820377935181e-05, -2.1039052974546231e-05, -2.7884246144565719e-06, 2.0623194846359774e-06, 5.5010591855292152e-07]
    top_gate_step_B_field_sweep_var(('E'), voltages, le_array, -0.20, 0.55, 6)
    
    le_array = [2.8872147858341384, 1.4844142937134688, -3.5681356758556428, -1.6089839132938182, 1.9261763124050821, 0.70851885882605226, -0.59511327902472388, -0.15933279716703977, 0.11499147895017386, 0.018023640714713227, -0.014041338877701726, -0.00058320985694616453, 0.00099677077156168988, -7.5618243768608198e-05, -2.6217166246436922e-05, 5.6634426623554682e-06, -8.5841828224469601e-07, 2.7414494188420501e-07]
    top_gate_step_B_field_sweep_var(('C'), voltages, le_array, -0.45, 0.3, 6)
    
    le_array = [1.8011214225584711, 0.34157390962799689, -2.8829700363255677, -0.41284554900496695, 2.0071658424581837, 0.23348581280264394, -0.75786221699657996, -0.067503576915902919, 0.16638449147669321, 0.0085250857417301928, -0.021492254894203089, 2.7030333104555563e-05, 0.0015345987456724555, -0.00010842897836591451, -4.2847045012635924e-05, 6.8398895644934003e-06, -8.0950405369759159e-07, 2.7428020240829887e-07]
    top_gate_step_B_field_sweep_var(('A'), voltages, le_array, -0.45, 0.3, 6)
        
    le_array = [-0.24837429751492376, 0.74662960133443812, 0.75388756574432936, -0.86906350499506579, -0.59934361330575547, 0.44468130852990767, 0.21384300348526536, -0.12952482337124993, -0.037544000848299335, 0.022710586854185307, 0.0026215841557996229, -0.0022416455280034318, 7.8426620036530051e-05, 8.8561524754453101e-05, -1.4123321118593649e-05, 9.7357812063617862e-07, -4.2477725456855097e-07, 2.0958445297615039e-07]
    top_gate_step_B_field_sweep_var(('B'), voltages, le_array, -0.45, 0.3, 6)
    
    le_array = [-0.17855665406658161, -0.28068666598167746, 0.079690265024948456, 0.2910348391599239, 0.020462818897240304, -0.12683859091491673, -0.018028380664479506, 0.029932096888560911, 0.0035555259484676212, -0.0039263623848528493, -0.00019247919509747386, 0.00023010328643151019, -1.1719566440806479e-07, -1.2488933199170218e-06, -1.275617054419993e-06, 4.0488132826110289e-07, -3.0834032428491116e-07, 1.4146975152802756e-07]
    top_gate_step_B_field_sweep_var(('D'), voltages, le_array, -0.65, 0.1, 6)
    
def field_sweep_combo():
    do_field_sweeps_at_V(ami, md, ('A'), voltages, V=-0.2236 , start=-2, stop=2, points=1000)
    do_field_sweeps_at_V(ami, md, ('B'), voltages, V=-0.2719 , start=-2, stop=2, points=1000)
    do_field_sweeps_at_V(ami, md, ('C'), voltages, V=-0.2055 , start=-2, stop=2, points=1000)
    do_field_sweeps_at_V(ami, md, ('D'), voltages, V=-0.3804 , start=-2, stop=2, points=1000)
    do_field_sweeps_at_V(ami, md, ('E'), voltages, V=0.1201 , start=-2, stop=2, points=1000)
    
    