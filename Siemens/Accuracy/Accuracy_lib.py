# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:47:03 2015

@author: home
"""

import numpy as np
import datetime
import time
from scipy import interpolate

class InterfileHeader:
    def __init__(self, filename):
        self.filename = filename
        self.tag = []
        self.value = []
        
        with open(self.filename) as f:
            for line in f:
                line = line.strip()
                parts = line.split("=")
                if parts[0]:
                    self.tag.append(parts[0])
                    if parts[1]:
                        self.value.append(parts[1])
                    else:
                        self.value.append('')
    
    def get(self, gettag):
        for x in range(len(self.tag)):
            if self.tag[x] == gettag:
                return self.value[x]

# http://stackoverflow.com/questions/2788871/python-date-difference-in-minutes
def diff_times_in_mins(CalTime, AcqTime):
    # caveat emptor - assumes t1 & t2 are python times, on the same day and t2 is after t1
    t1 = datetime.datetime.strptime( str(CalTime), '%Y:%m:%d:%H:%M:%S')
    t2 = datetime.datetime.strptime( str(AcqTime), '%Y:%m:%d:%H:%M:%S')
    d1_ts = time.mktime(t1.timetuple())
    d2_ts = time.mktime(t2.timetuple())
    return int(d2_ts-d1_ts) / 60

def CalDecayAndAverageActivity(InitialActivity, TimeDif, halflife, ImageDuration):
    Ao = InitialActivity * np.exp(-1*np.log(2) * TimeDif / halflife)
    Ave = ( Ao / np.log(2) ) * ( halflife / ImageDuration ) * ( 1 - np.exp(-1*np.log(2) * ImageDuration / halflife) )
    return Ave

def GetDeltaRAtNECRpeak(ActConc, DeltaRij, NECRpeakconc):
    f = interpolate.interp1d(ActConc[::-1], DeltaRij[::-1])
    return f(NECRpeakconc)   