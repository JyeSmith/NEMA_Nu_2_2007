# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 13:44:49 2015

@author: home
"""

import datetime
import numpy as np

def diff_times_in_mins(CalTime, AcqTime):
    # caveat emptor - assumes t1 & t2 are python times, on the same day and t2 is after t1
    t1 = datetime.datetime.strptime( str(CalTime), '%H:%M:%S').time()
    t2 = datetime.datetime.strptime( str(AcqTime), '%H:%M:%S').time()
    h1, m1, s1 = t1.hour, t1.minute, t1.second
    h2, m2, s2 = t2.hour, t2.minute, t2.second
    t1_secs = s1 + 60 * (m1 + 60*h1)
    t2_secs = s2 + 60 * (m2 + 60*h2)
    return( t2_secs - t1_secs) / 60. # return time difference in minutes
    
def SSRB(SinogramFile):
    
    sinogramAll = np.fromfile(SinogramFile, dtype=np.dtype('<i2'), count=(14*621*168*400))
    
    SSRB = np.zeros([109, 168, 400], dtype=np.dtype('<i2'))

    for ToF in range(14):
        sinogram = np.zeros(621*168*400, dtype=np.dtype('<i2'))
        sinogram[:] = sinogramAll[(ToF*621*168*400):((ToF+1)*621*168*400)]  
        sinogram.resize(621, 168, 400)
        
        if ToF != 13:
            SSRB[0:108, :, :] += sinogram[0:108, :, :]
            SSRB[6:102, :, :] += sinogram[109:205:, :, :]
            SSRB[6:102, :, :] += sinogram[206:302:, :, :]
            SSRB[17:91, :, :] += sinogram[303:377:, :, :]
            SSRB[17:91, :, :] += sinogram[378:452:, :, :]
            SSRB[28:80, :, :] += sinogram[453:505:, :, :]
            SSRB[28:80, :, :] += sinogram[506:558:, :, :]
            SSRB[39:69, :, :] += sinogram[559:589:, :, :]
            SSRB[39:69, :, :] += sinogram[590:620:, :, :]
        
        # subtract randoms bin
        if ToF == 13:
            SSRB[0:108, :, :] -= sinogram[0:108, :, :]
            SSRB[6:102, :, :] -= sinogram[109:205:, :, :]
            SSRB[6:102, :, :] -= sinogram[206:302:, :, :]
            SSRB[17:91, :, :] -= sinogram[303:377:, :, :]
            SSRB[17:91, :, :] -= sinogram[378:452:, :, :]
            SSRB[28:80, :, :] -= sinogram[453:505:, :, :]
            SSRB[28:80, :, :] -= sinogram[506:558:, :, :]
            SSRB[39:69, :, :] -= sinogram[559:589:, :, :]
            SSRB[39:69, :, :] -= sinogram[590:620:, :, :]
            
        return SSRB