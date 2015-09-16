# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:44:39 2015

@author: Jye Smith
"""

import numpy as np
import matplotlib.pyplot as plot
import Scatter_and_Randoms_lib as SR_lib
    
####################################################
Dir = '/media/home/TOSHIBA EXT/NEMAScatter/' ## images need to be seperate folders numbered 1 through 5 (# of sleves)
#SSRB_Prompts, SSRB_Randoms = SR_lib.SSRB(Dir+'uncompnemascatter_'+str(1)+'.s')
InitialActivity = 1052
####################################################

# [Plane, Angle, Bin]
# [621, 168, 400]
i = 109 # images (planes)
j = 45  # j = acquisition

RefImageHeader = SR_lib.InterfileHeader(Dir+'uncompnemascatter_1.s.hdr')

InitialActivity = float(RefImageHeader.get('tracer activity at time of injection (Bq):')) / 1000000
ActivitCalDateTime = RefImageHeader.get('%tracer injection date (yyyy:mm:dd):')+':'+\
                     RefImageHeader.get('%tracer injection time (hh:mm:ss GMT+00:00):')
halflife = int(RefImageHeader.get('isotope gamma halflife (sec):')) / 60.0

CTotij = np.zeros([i, j], dtype=int)
Crsij = np.zeros([i, j], dtype=int)
Crij = np.zeros([i, j], dtype=int)
SFij = np.zeros([i, j], dtype=float)

#ActConc = InitialActivity * np.exp(-1*np.log(2)*(np.arange(45)*20+5)/109) / (22000) # volume is defined as 22,000 cm3 as in the NEMA doc
ActConc = np.zeros([j], dtype=float)

for acquisition in range(j): #j
    print 'acquisition', acquisition
    
    ImageHeader = SR_lib.InterfileHeader(Dir+'uncompnemascatter_'+str(acquisition+1)+'.s.hdr')
    ImageDateTime = ImageHeader.get('%study date (yyyy:mm:dd):')+':'+\
                     ImageHeader.get('%study time (hh:mm:ss GMT+00:00):')
    TimeDif = SR_lib.diff_times_in_mins(ActivitCalDateTime, ImageDateTime)
    ImageDuration = int(ImageHeader.get('!image duration (sec):'))
    ActConc[acquisition] = SR_lib.CalDecayAndAverageActivity(InitialActivity, TimeDif, halflife, ImageDuration / 60) / 22000
    
    # Read sinogram and SSRB
    SSRB_Prompts, SSRB_Randoms = SR_lib.SSRB(Dir+'uncompnemascatter_'+str(acquisition+1)+'.s')
    
    for image in range(i): #i
        #print 'image', image
        # extract image plane to sum
        SSRB_Prompts_slice = SSRB_Prompts[image,:,:]
        SSRB_Randoms_slice = SSRB_Randoms[image,:,:]
        
        for projection in range(168):
            # Set all pixels outside +-12cm to zero
            # scale factor (mm/pixel) [1]:=2.005
            SSRB_Prompts_slice[projection, 0:200-int(120/2.005)]   *= 0
            SSRB_Prompts_slice[projection, 200+int(120/2.005):399] *= 0
            SSRB_Randoms_slice[projection, 0:200-int(120/2.005)]   *= 0
            SSRB_Randoms_slice[projection, 200+int(120/2.005):399] *= 0          
            # Center the max pixel
            line_Prompts = SSRB_Prompts_slice[projection, :]
            MaxIndex = np.argmax(line_Prompts)
            SSRB_Prompts_slice[projection,:] = np.roll(line_Prompts, 200 - MaxIndex)

        PromptsLineRsponse = np.squeeze(np.sum(SSRB_Prompts_slice, axis=0))
        RandomsLineRsponse = np.squeeze(np.sum(SSRB_Randoms_slice, axis=0))
        
        CTotij[image, acquisition] = np.sum(PromptsLineRsponse)
        
        Crsij[image, acquisition] = SR_lib.GetScatterPlusRandom(PromptsLineRsponse)
        
        Crij[image, acquisition] = np.sum(RandomsLineRsponse)
        
        SFij[image, acquisition] = ( Crsij[image, acquisition] - Crij[image, acquisition] ) / ( CTotij[image, acquisition] - Crij[image, acquisition] )


SFj = ( np.sum(Crsij, axis=0) - np.sum(Crij, axis=0) ) / ( np.sum(CTotij, axis=0) - np.sum(Crij, axis=0) )
RTotj = np.sum(CTotij, axis=0) / 600   
Rtj = np.sum(CTotij - Crsij, axis=0) / 600     
Rrj = np.sum(Crij, axis=0) / 600     
Rsj = np.sum(Crsij - Crij, axis=0) / 600     
RNECj = pow(np.int64(Rtj), 2) / (1.0 * RTotj)
    
NECRpeak, NECRpeakconc = SR_lib.GetPeak(ActConc, RNECj)
Rtpeak, Rtpeakconc = SR_lib.GetPeak(ActConc, Rtj)
SF_at_NECRpeak = SR_lib.GetScatterFractionAtNECRpeak(ActConc, SFj, NECRpeakconc)


plot.figure(1)

plot.subplot(121)
plot.plot(ActConc, RTotj, label='Total')
plot.plot(ActConc, Rtj, label='True')
plot.plot(ActConc, Rrj, label='Random')
plot.plot(ActConc, Rsj, label='Scatter')
plot.plot(ActConc, RNECj, label='NECR')
plot.title('Total, True, Random, Scatter and NEC Rates')
plot.xlabel('Average Concentration (MBq/cc)')
plot.ylabel('cps')
plot.ylim([0, 1.1 * np.max(Rtj)])
plot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plot.annotate('Peak Trues ['+str(int(Rtpeak/1000))+' kcps, '+str(round(Rtpeakconc, 4))+' MBq/cc] \n'+ \
            'Peak NECR ['+str(int(NECRpeak/1000))+' kcps, '+str(round(NECRpeakconc, 4))+' MBq/cc]', \
            xy=(0.2*np.max(ActConc), 0.25*np.max(RNECj)))

# cps/MBq vs plane 
plot.subplot(224)
plot.plot(ActConc, 100*SFj)
plot.title('Scatter Fraction')
plot.xlabel('Average Concentration (MBq/cc)')
plot.ylabel('Percentage (%)')
plot.ylim([0, 110 * np.max(SFj)])

plot.annotate('SF at peak NECR = '+str(round(100*SF_at_NECRpeak, 1))+' %', xy=(0.1 * np.max(ActConc), 50*SF_at_NECRpeak), xytext=(0.05 * np.max(ActConc), 50*SF_at_NECRpeak))

plot.show()

