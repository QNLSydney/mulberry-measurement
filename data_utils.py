import re
import os, os.path
import datetime

import qcodes as qc

def find_data(date, num):
    """
    Find the location of the data file taken on a given date
    and with a given ID
    """
    date = str(date)
    data_dir = os.path.join("data", date)
    files = os.listdir(data_dir)
    for file in files:
        if re.match("#{:03}".format(num), file):
            return os.path.join(data_dir, file)
    return None

def open_data(date, num):
    """
    Open a data file taken on a given date and with a given ID,
    and read in the data/metadata.
    """
    loc = find_data(date, num)
    if loc is None:
        return None
    data = qc.DataSet(location=loc)
    data.read()
    data.read_metadata()
    return data

def open_data_sequence(date, start_num, num_sweeps):
    """
    Open a series of data that may have been taken over multiple days

    Returns a list of datasets that have had their data loaded
    """

    data_list = []
    offs = 0
    for i in range(num_sweeps):
        data = open_data(date, start_num+i-offs)

        if data is None:
            date = datetime.datetime.strptime(date, "%Y-%m-%d")
            date += datetime.timedelta(days=1)
            date = date.strftime("%Y-%m-%d")
            offs = start_num+i-1
            data = open_data(date, start_num+i-offs)
            if data is None:
                print(f"Error opening data at: {date}, {start_num+i-offs}")
                break
        data_list.append(data)
    return data_list

def inst_param_val(data, instr, param):
    instrs_snapshot = data.metadata['station']['instruments']
    if instr not in instrs_snapshot:
        raise KeyError(f"Instrument {instr} not saved in snapshot")

    instr_snapshot = instrs_snapshot[instr]['parameters']
    if param not in instr_snapshot:
        raise KeyError(f'Parameter {param} not found in instrument {instr}')

    param_val = instr_snapshot[param]
    if 'value' not in param_val:
        raise KeyError(f'Parameter {param} found in {instr}, but value was not saved.')

    return param_val['value']

def detect_cycle(data):
    cycle_length = 0
    first_element = data[0]
    cycle_elements = [first_element]

    # Detect the cycle, and check that the data does truly cycle
    for i, elem in enumerate(data[1:], start=1):
        if cycle_length == 0:
            if elem == first_element:
                cycle_length = i
            else:
                cycle_elements.append(elem)
        else:
            if cycle_elements[i%cycle_length] != elem:
                raise ValueError(f"Unexpected value in cycle. "
                                 f"Expecting {cycle_elements[i%cycle_length]}, got {elem}.")

    if len(cycle_elements)%cycle_length != 0:
        raise ValueError("Last cycle incomplete? Perhaps the sweep was aborted!")

    return cycle_length, cycle_elements