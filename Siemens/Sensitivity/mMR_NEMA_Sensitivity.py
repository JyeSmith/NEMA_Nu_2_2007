# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 20:40:46 2015

@author: Jye Smith
"""

import numpy as np
import matplotlib.pyplot as plot
import NEMA_Sensitivity_lib as lib

   
####################################################

Dir = "mMR 10cm Sinograms/" 
## ordered 1 sleeve through to 5
SinogramFiles = ['2007Sensitivity_5.s.hdr', '2007Sensitivity_4.s.hdr', \
'2007Sensitivity_3.s.hdr', '2007Sensitivity_2.s.hdr', '2007Sensitivity_1.s.hdr']

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
    ImageHeader = lib.InterfileHeader(Dir+SinogramFiles[i])
    TotalNetTrues[i] = ImageHeader.get('%total net trues:')
    AcqTime[i] = ImageHeader.get('%study time (hh:mm:ss GMT+00:00):')#.strip(':')
    ImageDuration[i] = ImageHeader.get('!image duration (sec):')
    TracerInjectionTime = ImageHeader.get('%tracer injection time (hh:mm:ss GMT+00:00):')
    TracerActivityAtTimeOfInjection = float(ImageHeader.get('tracer activity at time of injection (Bq):')) / 1000000 # convert to MBq
    
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

#SSRB = lib.SSRB(Dir+AxialSinogram)
#
#AxialSensitivityProfile = np.sum(SSRB, axis=2)
#AxialSensitivityProfile = np.sum(AxialSensitivityProfile, axis=1)
#AxialSensitivityProfile *= Stot / np.sum(AxialSensitivityProfile)

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

## cps/MBq vs plane 
#plot.subplot(122)
#plot.plot(range(109), AxialSensitivityProfile)
#plot.title('Axial Sensitivity Profile')
#plot.xlabel('Plane Number (zero based)')
#plot.ylabel('cps / MBq')
#plot.xlim([0, 108])
#plot.ylim([0, 1.1 * np.max(AxialSensitivityProfile)])
plot.show()
