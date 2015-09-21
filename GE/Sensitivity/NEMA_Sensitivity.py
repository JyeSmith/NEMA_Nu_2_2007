# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 20:40:46 2015

@author: Jye Smith
"""

import numpy as np
import os
import matplotlib.pyplot as plot
import datetime
import dicom

def diff_times_in_mins(CalTime, AcqTime):
    # caveat emptor - assumes t1 & t2 are python times, on the same day and t2 is after t1
    t1 = datetime.datetime.strptime( str(CalTime), '%Y%m%d%H%M%S.%f').time()
    t2 = datetime.datetime.strptime( str(AcqTime), '%H%M%S.%W').time()
    h1, m1, s1 = t1.hour, t1.minute, t1.second
    h2, m2, s2 = t2.hour, t2.minute, t2.second
    t1_secs = s1 + 60 * (m1 + 60*h1)
    t2_secs = s2 + 60 * (m2 + 60*h2)
    return( t2_secs - t1_secs) / 60. # return time difference in minutes

###########################################################################

PathDicom = "10cm Sinograms/"
#PathDicom = "10cm Sinograms/"

lstFilesDCM = []  # create an empty list
for dirName, subdirList, fileList in os.walk(PathDicom):
    for filename in fileList:
            lstFilesDCM.append(os.path.join(dirName,filename))
			
Rcorrj = np.zeros(5, dtype=int)

for filenameDCM in lstFilesDCM:
    ds = dicom.read_file(filenameDCM)
    print ds[0x9,0x1062].value
    SeriesDescription = ds[0x8,0x103e].value
    frame = int(SeriesDescription.split(' ')[1])
    post_inj_activity = float( ds[0x9,0x103c].value )
    tracer_activity = float( ds[0x9,0x1038].value ) - post_inj_activity
    meas_datetime = ds[0x9,0x1039].value
    AcquisitionTime = ds[0x8,0x32].value
    acq_duration = int( ds[0x9,0x106d].value )
    total_delays = int( ds[0x9,0x1072].value )
    total_prompts = int( ds[0x9,0x1071].value ) - total_delays
    Rcorrj[frame-1] = (total_prompts / acq_duration) \
        * np.exp(np.log(2) * diff_times_in_mins(meas_datetime, AcquisitionTime) / 109.)

print 'Rcorrj', Rcorrj
    
## Take the log to make the fit linear   
PolyCoeff = np.polyfit([1,2,3,4,5], np.log(Rcorrj), 1)
p = np.poly1d(PolyCoeff)

Newx = np.arange(0, 6, 0.1)
Newy = np.exp( PolyCoeff[1] ) * np.exp( PolyCoeff[0] * Newx)

um = round(PolyCoeff[0] / (-2), 3)

Stot = np.exp( PolyCoeff[1] ) / tracer_activity
print 'Sensitivity =', int(Stot), 'cps/MBq'

plot.figure(1)

# Rate vs Thickness Fit
plot.subplot(121)
plot.plot([1,2,3,4,5], Rcorrj / tracer_activity, 'bs')
plot.plot(Newx, Newy / tracer_activity)
plot.title('Rate vs Thickness Fit')
plot.xlabel('Number of Sleeves')
plot.ylabel('cps/MBq')
plot.xlim([0, 6])

plot.annotate('Stot = '+str(int(Stot))+' cps/MBq \n'+ \
                r'$\mu_M$ = '+str(um)+r' $ mm^-$'+r'$^1$', xy=(1.5, Rcorrj[0] / tracer_activity))

plot.show()
        
#(0009, 103c) [PET post_inj_activity]             FL: 1.0099999904632568
        
#(0008, 0032) Acquisition Time                    TM: '111617.00'
#(0008, 103e) Series Description                  LO: 'frame 5'
#(0009, 102b) [PET scan_fov]                      SL: 70
#(0009, 102c) [PET axial_fov]                     SL: 153
#(0009, 1038) [PET tracer_activity]               FL: 6.739999771118164
#(0009, 1039) [PET meas_datetime]                 DT: '20150910104500.00'
#(0009, 1069) [PET slice_count]                   SL: 553
#(0009, 106d) [PET acq_duration]                  SL: 60
#(0009, 1071) [PET total_prompts]                 FD: 1821459.0
#(0009, 1072) [PET total_delays]                  FD: 221912.0
#(0009, 10df) [PET num_of_slices]                 US: 47
#(0021, 1001) [PET raw_data_type]                 US: 3
#(0021, 1002) [PET raw_data_size]                 UL: 0
#(0023, 1002) [PET raw_data_blob]                 OB: Array of 1586944 bytes
        


################## Axial Sensitivity Profile #########################
#
### dtype='<i2' < = LittleEndian, signed 16 bit. 14 = 13 ToF bins + randoms
#sinogramAll = np.fromfile(Dir+AxialSinogram, dtype=np.dtype('<i2'), count=(14*621*168*400))
#
#### [Plane, Angle, Bin]
##sinogram.resize(621*n, 168, 400)
##number of rings:=55
##%number of segments:=9
##%segment table:={109,97,97,75,75,53,53,31,31}
#
#SSRB = np.zeros([109, 168, 400], dtype=np.dtype('<i2'))
#
#for ToF in range(14):
#    print ToF
#    sinogram = np.zeros(621*168*400, dtype=np.dtype('<i2'))
#    sinogram[:] = sinogramAll[(ToF*621*168*400):((ToF+1)*621*168*400)]  
#    sinogram.resize(621, 168, 400)
#    
#    if ToF != 13:
#        SSRB[0:108, :, :] += sinogram[0:108, :, :]
#        SSRB[6:102, :, :] += sinogram[109:205:, :, :]
#        SSRB[6:102, :, :] += sinogram[206:302:, :, :]
#        SSRB[17:91, :, :] += sinogram[303:377:, :, :]
#        SSRB[17:91, :, :] += sinogram[378:452:, :, :]
#        SSRB[28:80, :, :] += sinogram[453:505:, :, :]
#        SSRB[28:80, :, :] += sinogram[506:558:, :, :]
#        SSRB[39:69, :, :] += sinogram[559:589:, :, :]
#        SSRB[39:69, :, :] += sinogram[590:620:, :, :]
#    
#    # subtract randoms bin
#    if ToF == 13:
#        SSRB[0:108, :, :] -= sinogram[0:108, :, :]
#        SSRB[6:102, :, :] -= sinogram[109:205:, :, :]
#        SSRB[6:102, :, :] -= sinogram[206:302:, :, :]
#        SSRB[17:91, :, :] -= sinogram[303:377:, :, :]
#        SSRB[17:91, :, :] -= sinogram[378:452:, :, :]
#        SSRB[28:80, :, :] -= sinogram[453:505:, :, :]
#        SSRB[28:80, :, :] -= sinogram[506:558:, :, :]
#        SSRB[39:69, :, :] -= sinogram[559:589:, :, :]
#        SSRB[39:69, :, :] -= sinogram[590:620:, :, :]
#
#AxialSensitivityProfile = np.sum(SSRB, axis=2)
#AxialSensitivityProfile = np.sum(AxialSensitivityProfile, axis=1)
#AxialSensitivityProfile *= Stot / np.sum(AxialSensitivityProfile)
#
## cps/MBq vs plane 
#plot.subplot(122)
#plot.plot(range(109), AxialSensitivityProfile)
#plot.title('Axial Sensitivity Profile')
#plot.xlabel('Plane Number (zero based)')
#plot.ylabel('cps / MBq')
#plot.xlim([0, 108])
#plot.ylim([0, 1.1 * np.max(AxialSensitivityProfile)])
#plot.show()
#




