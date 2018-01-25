from qcodes import VisaInstrument
from qcodes.utils.validators import Enum
from pyvisa.resources.serial import constants

class MulberryDriver(VisaInstrument):


    def __init__(self, name, address, **kwargs):

        super().__init__(name, address, terminator='\r', **kwargs)

        self.num_channels = 20
#
#        self.add_parameter('write',
#                           label='Write',
#                           get_parser=float,
#                           set_cmd='WRITE' + '{}',
#                           vals=Enum('somestring.....'))

        self.add_parameter('select',
                            label='Select',
                            set_cmd='SELECT {}',
                            vals=Enum('A', 'B', 'C', 'D', 'E'))


        self.add_function('clear', call_cmd='CLEAR')
        self.add_function('load', call_cmd='LOAD')
        self.add_function('clock', call_cmd='CLOCK')
        self.add_function('stop', call_cmd='STOP')
        
    def write(self, cmd):
        super().write(cmd)
        self.visa_handle.flush(64)