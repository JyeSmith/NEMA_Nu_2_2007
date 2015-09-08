# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 20:40:46 2015

@author: Jye Smith
"""

import numpy as np
import matplotlib.pyplot as plot
import datetime

def diff_times_in_mins(CalTime, AcqTime):
    # caveat emptor - assumes t1 & t2 are python times, on the same day and t2 is after t1
    t1 = datetime.datetime.strptime( str(CalTime), '%H:%M:%S').time()
    t2 = datetime.datetime.strptime( str(AcqTime), '%H:%M:%S').time()
    h1, m1, s1 = t1.hour, t1.minute, t1.second
    h2, m2, s2 = t2.hour, t2.minute, t2.second
    t1_secs = s1 + 60 * (m1 + 60*h1)
    t2_secs = s2 + 60 * (m2 + 60*h2)
    return( t2_secs - t1_secs) / 60. # return time difference in minutes

    
####################################################
Dir = "0cm Sinograms/" ## images need to be seperate folders numbered 1 through 5 (# of sleves)
## ordered 1 sleeve through to 5
SinogramFiles = ['uncompnemasensitivity_5.s.hdr', 'uncompnemasensitivity_4.s.hdr', \
'uncompnemasensitivity_3.s.hdr', 'uncompnemasensitivity_2.s.hdr', 'uncompnemasensitivity_1.s.hdr']
## Uncompressed sinogram for axial profile
AxialSinogram = 'uncompnemasensitivity_5.s'
####################################################

## %tracer injection time (hh:mm:ss GMT+00:00):=06:15:00
TracerInjectionTime = ''
## tracer activity at time of injection (Bq):=6.3e+006
TracerActivityAtTimeOfInjection = 0

# %total net trues:=11144102
TotalNetTrues = np.zeros(5, dtype=int)
## %study time (hh:mm:ss GMT+00:00):=07:10:10
AcqTime = ['', '', '', '', '']
## !image duration (sec):=300
ImageDuration = np.zeros(5, dtype=int) 

## Empty Rcorrj vector to be filled
Rcorrj = np.zeros(5, dtype=float)
              
for i in range(5):
    with open(Dir+SinogramFiles[i]) as f:
        
        for line in f:
            line = line.strip()
            parts = line.split("=")
            
            if parts[0] == '%total net trues:':
                TotalNetTrues[i] = parts[1]
                
            if parts[0] == '%study time (hh:mm:ss GMT+00:00):':
                AcqTime[i] = parts[1].strip(':')
                
            if parts[0] == '!image duration (sec):':
                ImageDuration[i] = parts[1]
                
            if parts[0] == '%tracer injection time (hh:mm:ss GMT+00:00):':
                TracerInjectionTime = parts[1]
                
            if parts[0] == 'tracer activity at time of injection (Bq):':
                moreparts = parts[1].split("e+")
                TracerActivityAtTimeOfInjection = float(moreparts[0]) * np.power(10, int(moreparts[1])) / 1000000

for i in range(5):
    Rcorrj[i] = ( TotalNetTrues[i] / ImageDuration[i] ) * np.exp(np.log(2) * diff_times_in_mins(TracerInjectionTime, AcqTime[i]) / 109.)
    
## Take the log to make the fit linear   
PolyCoeff = np.polyfit([1,2,3,4,5], np.log(Rcorrj), 1)
p = np.poly1d(PolyCoeff)

Stot = np.exp( PolyCoeff[1] ) / TracerActivityAtTimeOfInjection
print 'Sensitivity =', int(Stot), 'cps/MBq'


################# Axial Sensitivity Profile #########################

## dtype='<i2' < = LittleEndian, signed 16 bit. 14 = 13 ToF bins + randoms
sinogramAll = np.fromfile(Dir+AxialSinogram, dtype=np.dtype('<i2'), count=(14*621*168*400))

### [Plane, Angle, Bin]
#sinogram.resize(621*n, 168, 400)
#number of rings:=55
#%number of segments:=9
#%segment table:={109,97,97,75,75,53,53,31,31}

SSRB = np.zeros([109, 168, 400], dtype=np.dtype('<i2'))

for ToF in range(14):
    print ToF
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

AxialSensitivityProfile = np.sum(SSRB, axis=2)
AxialSensitivityProfile = np.sum(AxialSensitivityProfile, axis=1)
AxialSensitivityProfile *= Stot / np.sum(AxialSensitivityProfile)

plot.figure(1)

# Rate vs Thickness Fit
plot.subplot(121)
plot.plot([1,2,3,4,5], np.log(Rcorrj), 'bs', [0, 5], [PolyCoeff[1], 5 * PolyCoeff[0] + PolyCoeff[1]])
plot.title('Rate vs Thickness Fit')
plot.xlabel('Number of Sleeves')
plot.ylabel('Ln(Rcorr,j)')
plot.xlim([0, 6])

plot.annotate('Stot = '+str(int(Stot))+' cps/MBq', xy=(0, PolyCoeff[1]), xytext=(0.5, np.log(Rcorrj)[4]))
plot.annotate(p, xy=(0, PolyCoeff[1]), xytext=(2, np.log(Rcorrj)[0]))

# cps/MBq vs plane 
plot.subplot(122)
plot.plot(range(109), AxialSensitivityProfile)
plot.title('Axial Sensitivity Profile')
plot.xlabel('Plane Number (zero based)')
plot.ylabel('cps / MBq')
plot.xlim([0, 108])
plot.ylim([0, 1.1 * np.max(AxialSensitivityProfile)])
plot.show()





