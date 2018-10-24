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
for i in range(2):
    addr = "SR860_{}".format(i+1)
    lockins.append(SR860.SR860(addr, addr))
    station.add_component(lockins[-1])
    
## Connect to all dmms
#dmms = []
#for i in range(2):
#    addr = "dmm_{}".format(i+1)
#    dmms.append(Agilent_34400A.Agilent_34400A(addr, addr))
#    station.add_component(dmms[-1])
    
# Connect to yokogawa
yoko_t = GS200.GS200("yoko_t", "TCPIP::192.168.0.28::INSTR")
station.add_component(yoko_t)

# Connect to mulberry
md = MulberryDriver("md", "ASRL6::INSTR")
station.add_component(md)

# Connect to magnet backend
#ami = AMI430.AMI430("ami", "192.168.0.20", 7180)
#station.add_component(ami)

# Connect to yokogawa
yoko_mag = GS200.GS200("yoko_mag", "TCPIP::192.168.0.21::INSTR")
station.add_component(yoko_mag)
yoko_mag.current.step = 0.00001
yoko_mag.current.inter_delay = 1e-2

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
currents_X = [CurrentAmplifier(lockin.X, 1e6) for lockin in lockins[2:3]]
currents_Y = [CurrentAmplifier(lockin.Y, 1e6) for lockin in lockins[2:3]]
currents_R = [CurrentAmplifier(lockin.R, 1e6) for lockin in lockins[2:3]]

#currents_X =+ [lockins[0].X]
#currents_Y =+ [lockins[0].Y]
#currents_R =+ [lockins[0].R]

currents = [item for sublist in zip(currents_X, currents_Y, currents_R) for item in sublist]

resistances = [LockinResistance(lockin, 
                                input_imp=20, 
                                current_scale=1e6, 
                                voltage_scale=1) for lockin in lockins[2:3]]
#resistances += [LockinResistance(lockins[1])]
#resistances = [resistances[0], resistances[2], resistances[1]]
    
switched_resistances = []
for switch in ("A", "B", "C", "D", "E"):
    switched_resistances.extend(MBSwitchedParam(md, switch, resistance, 20) 
        for resistance in resistances)
switched_resistances = tuple(switched_resistances)

voltages_X = [VoltageAmplifier(lockin.X, 1) for lockin in lockins[0:2:1]]
voltages_Y = [VoltageAmplifier(lockin.Y, 1) for lockin in lockins[0:2:1]]
voltages_R = [VoltageAmplifier(lockin.R, 1) for lockin in lockins[0:2:1]]
voltages = [item for sublist in zip(voltages_X, voltages_Y, voltages_R) for item in sublist]
