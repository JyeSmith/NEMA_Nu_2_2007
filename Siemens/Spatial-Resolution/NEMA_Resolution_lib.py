# -*- coding: utf-8 -*-
"""
Created on Fri Sep 04 19:17:49 2015

@author: Jye Smith
"""

import numpy

def Calculate_x_Resolution(MaxIndices, point, LineResponse, ConstPixelDims, ConstPixelSpacing):
    
    MaximumIndices = LineResponse.argmax()
    PolyCoeff = numpy.polyfit([MaximumIndices-1, MaximumIndices, MaximumIndices+1], [LineResponse[MaximumIndices-1], LineResponse[MaximumIndices], LineResponse[MaximumIndices+1]], 2)
    p = numpy.poly1d(PolyCoeff)
    
    ## Differentiate the polynomial and solve for zero. This is the maximum index.
    MaxIndex = -PolyCoeff[1] / (2 * PolyCoeff[0])
    LineResponseMax = p(MaxIndex)
    
    FWHM1 = -1
    FWHM2 = -1
    ## find indices of FWHM
    for i in range(0, point):
        if LineResponse[i] >  LineResponseMax/2 and FWHM1 == -1:
            FWHM1Coeff = numpy.polyfit([i-1, i], [LineResponse[i-1], LineResponse[i]], 1)
            FWHM1 = ( LineResponseMax/2 - FWHM1Coeff[1] ) / FWHM1Coeff[0]
        if LineResponse[i] <  LineResponseMax/2 and FWHM2 == -1 and FWHM1 != -1:
            FWHM1Coeff = numpy.polyfit([i-1, i], [LineResponse[i-1], LineResponse[i]], 1)
            FWHM2 = ( LineResponseMax/2 - FWHM1Coeff[1] ) / FWHM1Coeff[0]
    
    FWHM = round( ConstPixelSpacing[0] * (FWHM2 - FWHM1), 2 )
    FWHM_in_pixels = round( FWHM2 - FWHM1, 2 )
    location = round( ConstPixelSpacing[1] * ( MaxIndices[1] - ConstPixelDims[1]/2. - int(point/2) + MaxIndex ), 1)
 
    ## Print points z location with respect to z direction
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print 'x location = ', location, ' mm'
    
    ## Print relevant NEMA information about acquisition
    print 'FWHM pixel width = ', FWHM_in_pixels, ". Must be atleast 3 e.g. 3 pixels accross the FWHM, if not acquire at a smaller resolution."
    print 'FWHM = ', FWHM, ' mm'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    
    return FWHM, FWHM_in_pixels, location
    
def Calculate_y_Resolution(MaxIndices, point, LineResponse, ConstPixelDims, ConstPixelSpacing):
    
    MaximumIndices = LineResponse.argmax()
    PolyCoeff = numpy.polyfit([MaximumIndices-1, MaximumIndices, MaximumIndices+1], [LineResponse[MaximumIndices-1], LineResponse[MaximumIndices], LineResponse[MaximumIndices+1]], 2)
    p = numpy.poly1d(PolyCoeff)
    
    ## Differentiate the polynomial and solve for zero. This is the maximum index.
    MaxIndex = -PolyCoeff[1] / (2 * PolyCoeff[0])
    LineResponseMax = p(MaxIndex)
    
    FWHM1 = -1
    FWHM2 = -1
    ## find indices of FWHM
    for i in range(0, point):
        if LineResponse[i] >  LineResponseMax/2 and FWHM1 == -1:
            FWHM1Coeff = numpy.polyfit([i-1, i], [LineResponse[i-1], LineResponse[i]], 1)
            FWHM1 = ( LineResponseMax/2 - FWHM1Coeff[1] ) / FWHM1Coeff[0]
        if LineResponse[i] <  LineResponseMax/2 and FWHM2 == -1 and FWHM1 != -1:
            FWHM1Coeff = numpy.polyfit([i-1, i], [LineResponse[i-1], LineResponse[i]], 1)
            FWHM2 = ( LineResponseMax/2 - FWHM1Coeff[1] ) / FWHM1Coeff[0]
        
    FWHM = round( ConstPixelSpacing[1] * (FWHM2 - FWHM1), 2 )
    FWHM_in_pixels = round( FWHM2 - FWHM1, 2 )
    location = round( ConstPixelSpacing[0] * ( ConstPixelDims[0]/2. - MaxIndices[0] + int(point/2) - MaxIndex ), 1)
    
    ## Print points z location with respect to z direction
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print 'y location = ', location, ' mm'
    
    ## Print relevant NEMA information about acquisition
    print 'FWHM pixel width = ', round( FWHM2 - FWHM1, 2 ), ". Must be atleast 3 e.g. 3 pixels accross the FWHM, if not acquire at a smaller resolution."
    print 'FWHM = ', round( ConstPixelSpacing[1] * (FWHM2 - FWHM1), 2 ), ' mm'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    return FWHM, FWHM_in_pixels, location

def Calculate_z_Resolution(MaxIndices, point, LineResponse, ConstPixelDims, ConstPixelSpacing):
    
    MaximumIndices = LineResponse.argmax()
    PolyCoeff = numpy.polyfit([MaximumIndices-1, MaximumIndices, MaximumIndices+1], [LineResponse[MaximumIndices-1], LineResponse[MaximumIndices], LineResponse[MaximumIndices+1]], 2)
    p = numpy.poly1d(PolyCoeff)
    
    ## Differentiate the polynomial and solve for zero. This is the maximum index.
    MaxIndex = -PolyCoeff[1] / (2 * PolyCoeff[0])
    LineResponseMax = p(MaxIndex)
    
    FWHM1 = -1
    FWHM2 = -1
    ## find indices of FWHM
    for i in range(0, point):
        if LineResponse[i] >  LineResponseMax/2 and FWHM1 == -1:
            FWHM1Coeff = numpy.polyfit([i-1, i], [LineResponse[i-1], LineResponse[i]], 1)
            FWHM1 = ( LineResponseMax/2 - FWHM1Coeff[1] ) / FWHM1Coeff[0]
        if LineResponse[i] <  LineResponseMax/2 and FWHM2 == -1 and FWHM1 != -1:
            FWHM1Coeff = numpy.polyfit([i-1, i], [LineResponse[i-1], LineResponse[i]], 1)
            FWHM2 = ( LineResponseMax/2 - FWHM1Coeff[1] ) / FWHM1Coeff[0]
           
    FWHM = round( ConstPixelSpacing[2] * (FWHM2 - FWHM1), 2 )
    FWHM_in_pixels = round( FWHM2 - FWHM1, 2 )
    FOV = ConstPixelSpacing[2] * ConstPixelDims[2]
    location = ConstPixelSpacing[2] * ( (MaxIndices[2]-int(point/2) + MaxIndex) )
    location2 = round( ConstPixelSpacing[2] * ( (MaxIndices[2]-int(point/2) + MaxIndex) - (ConstPixelDims[2]/2.) ), 1)
    location4 = round( ConstPixelSpacing[2] * ( (MaxIndices[2]-int(point/2) + MaxIndex) - (ConstPixelDims[2]/4.) ), 1)
    locationPercentage = round( 100. * location / FOV, 1)
    
    ## Print points z location with respect to z direction
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print 'z location from 1/2 axial FOV = ', location2, ' mm'
    print 'z location from 1/4 axial FOV = ', location4, ' mm'
    
    ## Print relevant NEMA information about acquisition
    print 'FWHM pixel width = ', FWHM_in_pixels, ". Must be atleast 3 e.g. 3 pixels accross the FWHM, if not acquire at a smaller resolution."
    print 'FWHM = ', FWHM, ' mm'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    return FWHM, FWHM_in_pixels, locationPercentage
