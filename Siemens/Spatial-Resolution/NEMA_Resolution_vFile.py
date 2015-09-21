# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 18:28:58 2015

@author: Jye Smith

NEMA NU 2-2007

Set 'PathDicom' to dir with dicom files. Can calculate FWHM of up to 3 points sources in a image.
"""

import numpy as np
from numpy import unravel_index
import matplotlib.pyplot as pyplot
import NEMA_Resolution_lib as lib

####################################################
Dir = "mMR vFile/" 

## File naming convention below can be either 'nemaresolution_1_2x_000_000.v.hdr' 
##                                         or '2007resolution_1_2x_000_000.v.hdr'.

j = 6 # Number of acquisitions
####################################################

RefImageHeader = lib.InterfileHeader(Dir+'2007resolution_1_2x_000_000.v.hdr')

MatrixSizex = int(RefImageHeader.get('matrix size[1]:'))
MatrixSizey = int(RefImageHeader.get('matrix size[2]:'))
i = int(RefImageHeader.get('matrix size[3]:'))
ConstPixelDims = [MatrixSizex, MatrixSizey, i]

x_scale = float(RefImageHeader.get('scale factor (mm/pixel) [1]:'))
y_scale = float(RefImageHeader.get('scale factor (mm/pixel) [2]:'))
slice_thickness = float(RefImageHeader.get('scale factor (mm/pixel) [3]:'))
ConstPixelSpacing = [x_scale, y_scale, slice_thickness]

Image = np.zeros([MatrixSizex, MatrixSizey, i], dtype='<f4')

for acquisition in range(j):
    print 'Reading file = ', Dir+'2007resolution_'+str(acquisition+1)+'_2x_000_000.v.hdr'
    ImageHeader = lib.InterfileHeader(Dir+'2007resolution_'+str(acquisition+1)+'_2x_000_000.v.hdr')
    DataFile = ImageHeader.get('!name of data file:')
    ImageVolume = np.fromfile(Dir+DataFile, dtype=np.dtype('<f4'), count=(MatrixSizex*MatrixSizey*i))
    ImageVolume.resize(i, MatrixSizex, MatrixSizey)
    for images in range(i):
        Image[:,:,images] += ImageVolume[images,:,:]

fig = pyplot.figure(0)
        
## Loop through for up to 3 points
for points in range(6):
    print 'Point found number', points+1      
    ## http://stackoverflow.com/questions/3584243/python-get-the-position-of-the-biggest-item-in-a-numpy-array
    MaxIndices = unravel_index(Image.argmax(), Image.shape) 
    
    ## calc 30mm cube size in pixels around point 
    pointx = int(round(30/ConstPixelSpacing[0]))
    pointy = int(round(30/ConstPixelSpacing[1])) 
    pointz = int(round(30/ConstPixelSpacing[2]))
    
    ## extract cube around point
    PointArray = Image[MaxIndices[0]-int(pointx/2): pointx+MaxIndices[0]-int(pointx/2), MaxIndices[1]-int(pointy/2): pointy+MaxIndices[1]-int(pointy/2), MaxIndices[2]-int(pointz/2): pointz+MaxIndices[2]-int(pointz/2)]
    
    print 'Line response counts = ', np.sum(PointArray), '. Must be atleast 100,000 counts.'
    
    if np.sum(PointArray) > 100000:
        
        ## Sum cube in to square
        FlatPointArray1 = np.sum(PointArray, axis=0)
        FlatPointArray2 = np.sum(PointArray, axis=1)
        
        ## Sum squares in to line response function
        xLineResponse = np.sum(FlatPointArray2, axis=1)
        yLineResponse = np.sum(FlatPointArray1, axis=1)
        zLineResponse = np.sum(FlatPointArray1, axis=0)

        ## Caclulate the FWHM of the line response function
        x_info = lib.Calculate_x_Resolution(MaxIndices, pointx, xLineResponse, ConstPixelDims, ConstPixelSpacing)
        y_info = lib.Calculate_y_Resolution(MaxIndices, pointy, yLineResponse, ConstPixelDims, ConstPixelSpacing)
        z_info = lib.Calculate_z_Resolution(MaxIndices, pointz, zLineResponse, ConstPixelDims, ConstPixelSpacing)
           
        Title = 'Location (' + str(x_info[2]) + ',' + str(y_info[2]) + ',' + str(z_info[2]) + '%) \n ' + \
                'FWHM (' + str(x_info[0]) + ',' + str(y_info[0]) + ',' + str(z_info[0]) + ') mm ' + \
                '(' + str(x_info[1]) + ',' + str(y_info[1]) + ',' + str(z_info[1]) + ') pixels'
                
        # Plot the grid
        ax = fig.add_subplot(2, 3, points+1)
        ax.set_title(Title, fontsize=10)
        ax.set_ylabel('Line response cnts = \n'+str(np.sum(PointArray)), fontsize=10)
        ax.imshow(np.sum(PointArray, axis=2), interpolation = 'none')
        
        ## Set that point source to zero in the image so it wont be selected again in the analysis.
        Image[MaxIndices[0]-int(pointx/2): pointx+MaxIndices[0]-int(pointx/2), MaxIndices[1]-int(pointy/2): pointy+MaxIndices[1]-int(pointy/2), MaxIndices[2]-int(pointz/2): pointz+MaxIndices[2]-int(pointz/2)] = 0

