# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:02:06 2018

@author: spauka
"""

import os, re

def find_data(date, num):
    date = str(date)
    data_dir = os.path.join("data", date)
    files = os.listdir(data_dir)
    for file in files:
        if re.match("#{:03}".format(num), file):
            return os.path.join(data_dir, file)