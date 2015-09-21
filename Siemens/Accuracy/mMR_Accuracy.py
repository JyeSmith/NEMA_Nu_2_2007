# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:44:39 2015

@author: Jye Smith
"""

import numpy as np
import matplotlib.pyplot as plot
import Accuracy_lib as lib

####################################################
Dir = 'mMR images/' 
j = 45 # Number of acquisitions
NECRpeakconc = 0.0212
####################################################

RefImageHeader = lib.InterfileHeader(Dir+'2007scatter_1_000_000.v.hdr')

InitialActivity = float(RefImageHeader.get('tracer activity at time of injection (Bq):')) / 1000000
ActivitCalDateTime = RefImageHeader.get('%tracer injection date (yyyy:mm:dd):')+':'+\
                     RefImageHeader.get('%tracer injection time (hh:mm:ss GMT+00:00):')
halflife = int(RefImageHeader.get('isotope gamma halflife (sec):')) / 60.0

MatrixSizex = int(RefImageHeader.get('matrix size[1]:'))
MatrixSizey = int(RefImageHeader.get('matrix size[2]:'))
i = int(RefImageHeader.get('matrix size[3]:'))
x_scale = float(RefImageHeader.get('scale factor (mm/pixel) [1]:'))
y_scale = float(RefImageHeader.get('scale factor (mm/pixel) [2]:'))
slice_thickness = float(RefImageHeader.get('scale factor (mm/pixel) [3]:'))

r = int( 180 / x_scale )

y,x = np.ogrid[-MatrixSizex/2: MatrixSizex/2, -MatrixSizey/2: MatrixSizey/2]
mask = x**2+y**2 <= r**2
mask *= 1

ActConc = np.zeros([j], dtype=float)
RROIij  = np.zeros([i, j], dtype=int)
RExtrij = np.zeros([i, j], dtype=int)
DeltaRij = np.zeros([i, j], dtype=float)
HighestDeltaRij = np.zeros([j], dtype=float)
LowestDeltaRij  = np.zeros([j], dtype=float)
MaxCentralEight = np.zeros([j], dtype=float)

for acquisition in range(j):
    print 'acquisition', acquisition
    ImageHeader = lib.InterfileHeader(Dir+'2007scatter_'+str(acquisition+1)+'_000_000.v.hdr')
    ImageDateTime = ImageHeader.get('%study date (yyyy:mm:dd):')+':'+\
                     ImageHeader.get('%study time (hh:mm:ss GMT+00:00):')
    print ImageDateTime
    TimeDif = lib.diff_times_in_mins(ActivitCalDateTime, ImageDateTime)
    ImageDuration = int(ImageHeader.get('!image duration (sec):'))
    ActConc[acquisition] = lib.CalDecayAndAverageActivity(InitialActivity, TimeDif, halflife, ImageDuration / 60) / 22000
    DataFile = ImageHeader.get('!name of data file:')
    ImageVolume = np.fromfile(Dir+DataFile, dtype=np.dtype('<f4'), count=(MatrixSizex*MatrixSizey*i))
    ImageVolume.resize(i, MatrixSizex, MatrixSizey)
    for image in range(i):
        RROIij[image, acquisition] = np.sum(mask * ImageVolume[image,:,:]) / ImageDuration

for acquisition in range(j):
    for image in range(i):
        
        RExtrij[image, acquisition] = ( ActConc[acquisition] / 3) * \
            ((RROIij[image, j-1]/ActConc[j-1]) + (RROIij[image, j-2]/ActConc[j-2]) + (RROIij[image, j-3]/ActConc[j-3]))
            
        DeltaRij[image, acquisition] = 100 * ((RROIij[image, acquisition] / float(RExtrij[image, acquisition])) - 1)
    
    HighestDeltaRij[acquisition] = np.max(DeltaRij[:, acquisition])
    LowestDeltaRij[acquisition]  = np.min(DeltaRij[:, acquisition])
    MaxCentralEight[acquisition]  = np.max(np.abs(DeltaRij[int(0.1*i):int(0.9*i), acquisition]))

MaxDeltaRAtNNECR = np.max(np.abs([lib.GetDeltaRAtNECRpeak(ActConc, HighestDeltaRij, NECRpeakconc), \
                                  lib.GetDeltaRAtNECRpeak(ActConc, LowestDeltaRij, NECRpeakconc)]))

plot.figure(1)

plot.subplot(121)
plot.plot(ActConc, HighestDeltaRij)
plot.plot(ActConc, LowestDeltaRij)
plot.plot(ActConc, MaxCentralEight, label='Max over central 80% of axial FOV')
plot.title('Min/Max Error vs Average Concentration')
plot.xlabel('Average Concentration (MBq/cc)')
plot.ylabel('Percentage (%)')
plot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plot.annotate('Maximum bias at peak NECR = '+str(round(MaxDeltaRAtNNECR, 2))+' %', xy=(0.2*np.max(ActConc), 0))

plot.subplot(224)
plot.plot(ActConc, np.sum(RROIij, axis=0))
plot.plot(ActConc, np.sum(RExtrij, axis=0))
plot.title('Respose vs Average Concentration')
plot.xlabel('Average Concentration (MBq/cc)')
plot.ylabel('Response (cps)')
