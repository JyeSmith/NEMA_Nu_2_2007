# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:47:03 2015

@author: home
"""

import numpy as np
from scipy import interpolate
import datetime
import time

def SSRB(SinogramFile):
    ## dtype='<i2' < = LittleEndian, signed 16 bit. 14 = 13 ToF bins + randoms
    sinogramAll = np.fromfile(SinogramFile, dtype=np.dtype('<i2'), count=(14*621*168*400))
    
    ### [Plane, Angle, Bin]
    #sinogram.resize(621*n, 168, 400)
    #number of rings:=55
    #%number of segments:=9
    #%segment table:={109,97,97,75,75,53,53,31,31}
    
    SSRB_Prompts = np.zeros([109, 168, 400], dtype=np.dtype('<i2'))
    SSRB_Randoms = np.zeros([109, 168, 400], dtype=np.dtype('<i2'))
    
    for ToF in range(14):
        #print SinogramFile, 'ToF bin', ToF
        sinogram = np.zeros(621*168*400, dtype=np.dtype('<i2'))
        sinogram[:] = sinogramAll[(ToF*621*168*400):((ToF+1)*621*168*400)]  
        sinogram.resize(621, 168, 400)
        
        if ToF != 13:
            SSRB_Prompts[0:108, :, :] += sinogram[0:108, :, :]
            SSRB_Prompts[6:102, :, :] += sinogram[109:205:, :, :]
            SSRB_Prompts[6:102, :, :] += sinogram[206:302:, :, :]
            SSRB_Prompts[17:91, :, :] += sinogram[303:377:, :, :]
            SSRB_Prompts[17:91, :, :] += sinogram[378:452:, :, :]
            SSRB_Prompts[28:80, :, :] += sinogram[453:505:, :, :]
            SSRB_Prompts[28:80, :, :] += sinogram[506:558:, :, :]
            SSRB_Prompts[39:69, :, :] += sinogram[559:589:, :, :]
            SSRB_Prompts[39:69, :, :] += sinogram[590:620:, :, :]
        
        # subtract randoms bin
        if ToF == 13:
            SSRB_Randoms[0:108, :, :] += sinogram[0:108, :, :]
            SSRB_Randoms[6:102, :, :] += sinogram[109:205:, :, :]
            SSRB_Randoms[6:102, :, :] += sinogram[206:302:, :, :]
            SSRB_Randoms[17:91, :, :] += sinogram[303:377:, :, :]
            SSRB_Randoms[17:91, :, :] += sinogram[378:452:, :, :]
            SSRB_Randoms[28:80, :, :] += sinogram[453:505:, :, :]
            SSRB_Randoms[28:80, :, :] += sinogram[506:558:, :, :]
            SSRB_Randoms[39:69, :, :] += sinogram[559:589:, :, :]
            SSRB_Randoms[39:69, :, :] += sinogram[590:620:, :, :]
            
            return SSRB_Prompts, SSRB_Randoms
            
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
         
def GetScatterPlusRandom(PromptsLineRsponse):
    # scale factor (mm/pixel) [1]:=2.005
    
    CL, CL1, CL2 = [200-20/2.005, np.floor(200-20/2.005), np.ceil(200-20/2.005)] 
    CL1y = PromptsLineRsponse[CL1]
    CL2y = PromptsLineRsponse[CL2]
    PolyCoeff = np.polyfit([CL1, CL2], [CL1y, CL2y], 1)        
    CLy = CL*PolyCoeff[0]+PolyCoeff[1]
    
    CR, CR1, CR2 = [200+20/2.005, np.floor(200+20/2.005), np.ceil(200+20/2.005)]    
    CR1y = PromptsLineRsponse[CR1]
    CR2y = PromptsLineRsponse[CR2]
    PolyCoeff = np.polyfit([CR1, CR2], [CR1y, CR2y], 1)        
    CRy = CR*PolyCoeff[0]+PolyCoeff[1]
    
    SR_Under_Peak = np.average([CLy, CRy]) * (CR - CL)
    ScatterLeft  = np.sum(PromptsLineRsponse[0:CL1]) + np.average([CL1y, CLy]) * (CL - CL1)
    Scatterright = np.sum(PromptsLineRsponse[CR2:399]) + np.average([CR2y, CRy]) * (CR2 - CR)
    
    return ScatterLeft + SR_Under_Peak + Scatterright

def GetPeak(x, y):
    MaximumIndices = y.argmax()
    PolyCoeff = np.polyfit([x[MaximumIndices-1], x[MaximumIndices], x[MaximumIndices+1]], [y[MaximumIndices-1], y[MaximumIndices], y[MaximumIndices+1]], 2)
    p = np.poly1d(PolyCoeff)
    
    ## Differentiate the polynomial and solve for zero. This is the maximum index.
    MaxIndex = -PolyCoeff[1] / (2 * PolyCoeff[0])
    
    return p(MaxIndex), MaxIndex

def GetScatterFractionAtNECRpeak(ActConc, SFj, NECRpeakconc):
    f = interpolate.interp1d(ActConc[::-1], SFj[::-1])
    return f(NECRpeakconc)    
    
def CalDecayAndAverageActivity(InitialActivity, TimeDif, halflife, ImageDuration):
    Ao = InitialActivity * np.exp(-1*np.log(2) * TimeDif / halflife)
    Ave = ( Ao / np.log(2) ) * ( halflife / ImageDuration ) * ( 1 - np.exp(-1*np.log(2) * ImageDuration / halflife) )
    return Ave
    
# http://stackoverflow.com/questions/2788871/python-date-difference-in-minutes
def diff_times_in_mins(CalTime, AcqTime):
    # caveat emptor - assumes t1 & t2 are python times, on the same day and t2 is after t1
    t1 = datetime.datetime.strptime( str(CalTime), '%Y:%m:%d:%H:%M:%S')
    t2 = datetime.datetime.strptime( str(AcqTime), '%Y:%m:%d:%H:%M:%S')
    d1_ts = time.mktime(t1.timetuple())
    d2_ts = time.mktime(t2.timetuple())
    return int(d2_ts-d1_ts) / 60