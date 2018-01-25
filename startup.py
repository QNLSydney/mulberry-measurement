# -*- coding: utf-8 -*-

import qcodes as qc
from qcodes.instrument_drivers.stanford_research import SR860
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

# Connect to mulberry
md = MulberryDriver("md", "ASRL6::INSTR")
station.add_component(md)

from measurements import *

resistances = gen_resistances_param(lockins[:3])
lsr = LockinSwitchResistances(md, ("B", "C", "D", "E"), resistances, 10)
ft = FridgeTemps("BlueFors_LD", 
     "https://qphys1114.research.ext.sydney.edu.au/therm_flask/BlueFors_LD/data/?current")
t = TimeParam(60)



#h5fmt = hdf5_format.HDF5Format()
#qc.DataSet.default_formatter = h5fmt