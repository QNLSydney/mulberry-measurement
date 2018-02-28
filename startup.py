# -*- coding: utf-8 -*-

import qcodes as qc
from qcodes.instrument_drivers.stanford_research import SR860
from qcodes.instrument_drivers.american_magnetics import AMI430
from qcodes.instrument_drivers.agilent import Agilent_34400A
from qcodes.instrument_drivers.yokogawa import GS200
from MulberryDriver import MulberryDriver

# Close any instruments that may already be open
instruments = list(qc.Instrument._all_instruments.keys())
for instrument in instruments:
    instr = qc.Instrument._all_instruments.pop(instrument)
    instr = instr()
    instr.close()
    
# Create a new QCodes Station
station = qc.Station()

# Connect to all lockins
lockins = []
for i in range(4):
    addr = "SR860_{}".format(i+1)
    lockins.append(SR860.SR860(addr, addr))
    station.add_component(lockins[-1])
    
# Connect to all dmms
dmms = []
for i in range(2):
    addr = "dmm_{}".format(i+1)
    dmms.append(Agilent_34400A.Agilent_34400A(addr, addr))
    station.add_component(dmms[-1])
    
# Connect to yokogawa
#yoko_tg = GS200.GS200("yoko_tg", "yoko_tg")
#station.add_component(yoko_tg)

# Connect to mulberry
md = MulberryDriver("md", "ASRL6::INSTR")
station.add_component(md)

# Connect to magnet backend
ami = AMI430.AMI430("ami", "192.168.0.20", 7180)
station.add_component(ami)

from measurements import *

ft = FridgeTemps("BlueFors_LD", 
     "https://qphys1114.research.ext.sydney.edu.au/therm_flask/BlueFors_LD/data/?current")
t = TimeParam(20)

# Source_Drain Currents
for lockin in lockins[::2]:
    lockin.sine_outdc.scale = 100
V_sd_params = [lockin.sine_outdc for lockin in lockins[::2]]
V_sd = SourceDrainVoltages(V_sd_params)

# Lockin voltage/current measurements
currents_X = [CurrentAmplifier(lockin.X, 1e6) for lockin in lockins[::2]]
currents_Y = [CurrentAmplifier(lockin.Y, 1e6) for lockin in lockins[::2]]
currents_R = [CurrentAmplifier(lockin.R, 1e6) for lockin in lockins[::2]]

currents_X += [lockins[1].X]
currents_Y += [lockins[1].Y]
currents_R += [lockins[1].R]

currents = [item for sublist in zip(currents_X, currents_Y, currents_R) for item in sublist]

resistances = [LockinResistance(lockin, 
                                input_imp=20, 
                                current_scale=1e6, 
                                voltage_scale=1) for lockin in lockins[::2]]
resistances += [LockinResistance(lockins[1])]
    
switched_resistances = []
for switch in ("A", "B", "C", "D", "E"):
    switched_resistances.extend(MBSwitchedParam(md, switch, resistance, 20) 
        for resistance in resistances)
switched_resistances = tuple(switched_resistances)

voltages_X = [lockin.X for lockin in lockins]
voltages_Y = [lockin.Y for lockin in lockins]
voltages_R = [lockin.R for lockin in lockins]
voltages = [item for sublist in zip(voltages_X, voltages_Y, voltages_R) for item in sublist]

#h5fmt = hdf5_format.HDF5Format()
#qc.DataSet.default_formatter = h5fmt
