# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 18:28:58 2015

@author: Jye Smith

NEMA NU 2-2007

Set 'PathDicom' to dir with dicom files. Can calculate FWHM of up to 3 points sources in a image.
"""

## https://pyscience.wordpress.com/2014/09/08/dicom-in-python-importing-medical-image-data-into-numpy-with-pydicom-and-vtk/

import dicom
import os
import numpy
from numpy import unravel_index
import matplotlib.pyplot as pyplot
import NEMA_Resolution_lib

#PathDicom = "half_FOV/"
#PathDicom = "half_FOV_HD/"
#PathDicom = "quater_FOV/"
PathDicom = "quater_FOV_HD/"

lstFilesDCM = []  # create an empty list
for dirName, subdirList, fileList in os.walk(PathDicom):
    for filename in fileList:
            lstFilesDCM.append(os.path.join(dirName,filename))
			
# Read the first file to get header information
RefDs = dicom.read_file(lstFilesDCM[0])

# Load dimensions based on the number of rows, columns, and slices (along the Z axis)
ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))
print 'ConstPixelDims = ', ConstPixelDims[0], ConstPixelDims[1], ConstPixelDims[2]

# Load spacing values (in mm)
ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))
print 'ConstPixelSpacing = ', ConstPixelSpacing[0], ConstPixelSpacing[1], ConstPixelSpacing[2]

print 'x FOV = ', round( ConstPixelSpacing[0] * ConstPixelDims[0], 2 ), ' mm'
print 'y FOV = ', round( ConstPixelSpacing[1] * ConstPixelDims[1], 2 ), ' mm'
print 'axial FOV = ', round( ConstPixelSpacing[2] * ConstPixelDims[2], 2 ), ' mm'
    
# The array is sized based on 'ConstPixelDims'
ArrayDicom = numpy.zeros(ConstPixelDims, dtype=float)

# loop through all the DICOM files and copy yo numpy array
for filenameDCM in lstFilesDCM:
    ds = dicom.read_file(filenameDCM)
    RescaleIntercept = float( ds[0x28,0x1052].value ) ## (0028, 1052) Rescale Intercept DS: '0'
    RescaleSlope = float( ds[0x28,0x1053].value ) ## (0028, 1053) Rescale Slope DS: '2.97373'
    ArrayDicom[:, :, ds[0x20,0x13].value - 1] = ds.pixel_array * RescaleSlope + RescaleIntercept ## [0x20,0x13] is the 'Instance Number'.  This will order the image correctly in the array.
        
fig = pyplot.figure()
ax = fig.add_subplot(2, 2, 1)
ax.imshow(numpy.sum(ArrayDicom, axis=2), interpolation = 'none')
        
## Loop through for up to 6 points
for points in numpy.arange(0, 3):
    print
    print 'Point found number', points+1      
    ## http://stackoverflow.com/questions/3584243/python-get-the-position-of-the-biggest-item-in-a-numpy-array
    MaxIndices = unravel_index(ArrayDicom.argmax(), ArrayDicom.shape) 
    
    ## calc 30mm cube size in pixels around point 
    pointx = int(round(30/ConstPixelSpacing[0]))
    pointy = int(round(30/ConstPixelSpacing[1])) 
    pointz = int(round(30/ConstPixelSpacing[2]))
    
    ## extract cube around point
    PointArray = ArrayDicom[MaxIndices[0]-int(pointx/2): pointx+MaxIndices[0]-int(pointx/2), MaxIndices[1]-int(pointy/2): pointy+MaxIndices[1]-int(pointy/2), MaxIndices[2]-int(pointz/2): pointz+MaxIndices[2]-int(pointz/2)]
    
    print 'Line response counts = ', numpy.sum(PointArray), '. Must be atleast 100,000 counts.'
    
    if numpy.sum(PointArray) > 100000:
        
        ## Sum cube in to square
        FlatPointArray1 = numpy.sum(PointArray, axis=0)
        FlatPointArray2 = numpy.sum(PointArray, axis=1)
        
        ## Sum squares in to line response function
        xLineResponse = numpy.sum(FlatPointArray2, axis=1)
        yLineResponse = numpy.sum(FlatPointArray1, axis=1)
        zLineResponse = numpy.sum(FlatPointArray1, axis=0)

        ## Caclulate the FWHM of the line response function
        x_info = NEMA_Resolution_lib.Calculate_x_Resolution(MaxIndices, pointx, xLineResponse, ConstPixelDims, ConstPixelSpacing)
        y_info = NEMA_Resolution_lib.Calculate_y_Resolution(MaxIndices, pointy, yLineResponse, ConstPixelDims, ConstPixelSpacing)
        z_info = NEMA_Resolution_lib.Calculate_z_Resolution(MaxIndices, pointz, zLineResponse, ConstPixelDims, ConstPixelSpacing)
           
        Title = 'Location (' + str(x_info[2]) + ',' + str(y_info[2]) + ',' + str(z_info[2]) + '%) \n ' + \
                'FWHM (' + str(x_info[0]) + ',' + str(y_info[0]) + ',' + str(z_info[0]) + ') mm ' + \
                '(' + str(x_info[1]) + ',' + str(y_info[1]) + ',' + str(z_info[1]) + ') pixels'
                
        # Plot the grid
        ax = fig.add_subplot(2, 2, points+2)
        ax.set_title(Title)
        ax.set_ylabel('Line response cnts = \n'+str(numpy.sum(PointArray)))
        ax.imshow(numpy.sum(PointArray, axis=2), interpolation = 'none')
        
        ## Set that point source to zero in the image so it wont be selected again in the analysis.
        ArrayDicom[MaxIndices[0]-int(pointx/2): pointx+MaxIndices[0]-int(pointx/2), MaxIndices[1]-int(pointy/2): pointy+MaxIndices[1]-int(pointy/2), MaxIndices[2]-int(pointz/2): pointz+MaxIndices[2]-int(pointz/2)] = 0
        
