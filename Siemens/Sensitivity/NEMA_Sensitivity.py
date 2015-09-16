# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 20:40:46 2015

@author: Jye Smith
"""

import numpy as np
import matplotlib.pyplot as plot
import NEMA_Sensitivity_lib as lib

   
####################################################
#Dir = "0cm Sinograms/" ## images need to be seperate folders numbered 1 through 5 (# of sleves)
Dir = "10cm Sinograms/" ## images need to be seperate folders numbered 1 through 5 (# of sleves)
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
    
    Rcorrj[i] = ( TotalNetTrues[i] / ImageDuration[i] ) * np.exp(np.log(2) * lib.diff_times_in_mins(TracerInjectionTime, AcqTime[i]) / 109.)
    
## Take the log to make the fit linear   
PolyCoeff = np.polyfit([1,2,3,4,5], np.log(Rcorrj), 1)
p = np.poly1d(PolyCoeff)

Newx = np.arange(0, 6, 0.1)
Newy = np.exp( PolyCoeff[1] ) * np.exp( PolyCoeff[0] * Newx)

um = round(PolyCoeff[0] / (-2), 3)

Stot = np.exp( PolyCoeff[1] ) / TracerActivityAtTimeOfInjection
print 'Sensitivity =', int(Stot), 'cps/MBq'


################# Axial Sensitivity Profile #########################

SSRB = lib.SSRB(Dir+AxialSinogram)

AxialSensitivityProfile = np.sum(SSRB, axis=2)
AxialSensitivityProfile = np.sum(AxialSensitivityProfile, axis=1)
AxialSensitivityProfile *= Stot / np.sum(AxialSensitivityProfile)

################# Plot #############################################
plot.figure(2)

# Rate vs Thickness Fit
plot.subplot(121)
plot.plot([1,2,3,4,5], Rcorrj / TracerActivityAtTimeOfInjection, 'bs')
plot.plot(Newx, Newy / TracerActivityAtTimeOfInjection)
plot.title('Rate vs Thickness Fit')
plot.xlabel('Number of Sleeves')
plot.ylabel('cps/MBq')
plot.xlim([0, 6])

plot.annotate('Stot = '+str(int(Stot))+' cps/MBq \n'+ \
                r'$\mu_M$ = '+str(um)+r' $ mm^-$'+r'$^1$', xy=(1.5, Rcorrj[0] / TracerActivityAtTimeOfInjection))

# cps/MBq vs plane 
plot.subplot(122)
plot.plot(range(109), AxialSensitivityProfile)
plot.title('Axial Sensitivity Profile')
plot.xlabel('Plane Number (zero based)')
plot.ylabel('cps / MBq')
plot.xlim([0, 108])
plot.ylim([0, 1.1 * np.max(AxialSensitivityProfile)])
plot.show()
