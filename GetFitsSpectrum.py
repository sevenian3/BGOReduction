# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 13:52:11 2022

Based on tutrial at
https://docs.astropy.org/en/stable/io/fits/index.html

Order of processing BGO spectra:
    
    1) Get the spectrum from the FITS file:
       spectrum = getFitsSpectrum()
       - automatically bias subtracted and flat-field corrected
    2) Calibrate the wavelength scale with ar-ne calib.ini data:
        wave = waveCalib(spectrum)
    3) Correct for instrumental response and rectify:
        [wave, spectrum] = specNorm2(wave, spectrum)

@author: Ian Short
"""

import numpy as np
import math
from astropy.io import fits
import matplotlib.pyplot as plt

def getFitsSpectrum(subPath, file):
    
    #file = "LIVEOBS-ID18104-HR3309-OC159232-GR5217-ALPY.FIT"
    path = "C:/Users/s7970488/Documents/BGO/" + subPath + "/"
    
    with fits.open(path+file) as hdul:
    
        hdul.info()
        hdr = hdul[0].header

        hdr['OBJECT']
        hdr['DATE-OBS']
        
        #data is a numpy ndarray object
        data = np.array(hdul[0].data)
        data.shape
        data.dtype.name
        
        #number of pixels in dispersion direction:
        numCols = hdr['NAXIS1']
        #number of pixels perpendicular to disperion direction
        numRows = hdr['NAXIS2']
        print("numCols ", numCols, " numRows ", numRows)

# FIRST!  Subtract off zero-based ADU correction (FITS header 'PEDESTAL' field

    pedestal = int(hdr['PEDESTAL'])
    buffer = 10
    #Want a +ve value that we explicitly subtract:
    if (pedestal < 0):
        pedestal = -1 * pedestal
    pedestal = pedestal - buffer #safety

    #Subtract off pedestal:
        
    #For testing:
    #pedestal = 0
    #buffer = 0
    print("pedestal ", pedestal)   
    for col in range(0, numCols):
        for row in range(0, numRows):
            if ((col == 500) and (row == 10)):
                print('before: ', data[10][500])
            if (data[row][col] <= pedestal):
                data[row][col] = buffer
            else:
                data[row][col] = data[row][col] - pedestal
            if ((col == 500) and (row == 10)):
                print('after: ', data[10][500])
    

#THen:
# Fit (0th, 0th) order ("flat") scattered light and other background and 
# subtract from all pixels:
# Fit to two long narrow areas just inside edges of CCD:
    
    count = 0 #pixel number accumulator
    total = 0.0 #total counts accumulator
    for col in range(10, numCols-10):
        for row in range(10, 50):
            total += data[row][col]
            count+=1
        for row in range(numRows-50, numRows-10):
            total += data[row][col]
            count+=1    
            
    mean = float(total) / float(count)
    background0 = int(mean) - 1
    
    print("count ", count, " total ", total)
    print("background0 ", background0)

        
    #Subtract off background:
    for col in range(0, numCols):
        for row in range(0, numRows):
            if (data[row][col] <= background0):
                data[row][col] = 0
            else:
                data[row][col] = data[row][col] - background0
                
            
    #Plot up three representative crosss-sections perpendicular to dispersion to find out where spectrum is on the chip:
        
    refCols = [100, int(round(numCols/2)), numCols-100]
    #print(refPoints)
    cross = np.array([ [0.0 for i in range(numRows)] for j in range(3) ])
    
    for row in range(numRows):
        cross[0][row] = data[row][refCols[0]]
        cross[1][row] = data[row][refCols[1]]
        cross[2][row] = data[row][refCols[2]]
        
    maxCounts = np.amax(data)
    print("maxCounts", maxCounts)
    maxCountsMid = np.amax(cross[1])
    #Exact value of cross-dispersion profile weighting doesn't really matter:
    maxWeight = math.sqrt(maxCountsMid)
        
    #Find the approximate centre of the cross-dispersion profile at the 3
    #reference columns by just locating the pixel with the most counts:
        
    max0 = cross[0].argmax()
    max1 = cross[1].argmax()   
    max2 = cross[2].argmax()
    
    #Find the approximate half-power widths of the 3 reference columns:
    
    halfVal0 = int(round(cross[0][max0]/2.0))
    half0 = np.where(cross[0] >= halfVal0)
    halfVal1 = int(round(cross[1][max1]/2.0))
    half1 = np.where(cross[1] >= halfVal1)
    halfVal2 = int(round(cross[2][max2]/2.0))
    half2 = np.where(cross[2] >= halfVal2)
    
    print("initial values: max0 ", max0, " max1 ", max1, " max2 ", max2)
      
    #print("int(cross[0][max0]) ", int(cross[0][max0]) , "int(cross[1][max1]) ", int(cross[1][max1]), "int(cross[2][max2]) ") #, int(cross[2][max2]))                                                                      )
    #print("halfVal0 ", halfVal0, " halfVal1 ", halfVal1, " halfVal2 ", halfVal2)
    #print("half0 ", half0)
    #print("half1 ", half1)
    #print("half2 ", half2)
    
    #np.where returns a 2-element tuple where the first element [0] is the 1D array
    # we need to access
    
    max0 = int(round(((half0[0][0]+half0[0][len(half0[0])-1])/2)))
    max1 = int(round(((half1[0][0]+half1[0][len(half1[0])-1])/2)))
    max2 = int(round(((half2[0][0]+half2[0][len(half2[0])-1])/2)))
    
    print("max0 ", max0, " max1 ", max1, " max2 ", max2)
    
    plt.figure()
    plt.subplot(1, 1, 1)    
    plt.title = ""
    plt.ylabel('Counts')
    plt.xlabel('Pixel row (x)')
    plotLabel = file[0:10]
        
    #plt.plot(cross[0], color='blue')
    #plt.plot(cross[1], color='orange')
    #plt.plot(cross[2], color='green')
    
    #plt.plot([max0, max0], [0, maxCounts], color='blue')
    #plt.plot([max1, max1], [0, maxCounts], color='orange')
    #plt.plot([max2, max2], [0, maxCounts], color='green')

#cross-dispersion profile centroid is a function of column number due to spectrum tilt:
# "y" is pixel row and "x" is pixel column
    slope = (max2-max0)/(refCols[2]-refCols[0])
    #2 y-intercept estimates:
    b1 = max1 - (slope*refCols[1])
    b2 = max2 - (slope*refCols[2])
    #mean y-interecept
    b = (b1 + b2)/2.0
    #print("b ", b)
    #Calculate "tilt" of spectrum centroid with respect to pixel rows, in RAD:
    theta = math.atan(slope) # Radians!
    print("Spectrum is tilted by ", theta*180.0/math.pi, " deg wrt pixel rows")
    
    #Co-add over subset of rows, deltaRows, centered on centroid
    deltaRows = 30
    halfDelta = int(deltaRows/2)
    #Number of columns in co-adding region correspond to tilt of spectrum:
    driftCols = deltaRows * math.tan(theta) #expects Radians!
    print("Due to tilt, one wavelength element corresponds to ", driftCols, " pixel columns") 

#Calculate centroid pixel (y) as a function of pixel column (x) across chip: 
    centroid0 = ( max0 + max1 + max2 ) / 3.0
    centroid0 = int(round(centroid0))
    #print("centroid0 ", centroid0) # for checking
    #centroid = [0.0 for x in range(numCols)]
    centroid = [(slope*x)+b for x in range(numCols)]
    centroid = [int(round(x)) for x in centroid]
    print("centroid[550] ", centroid[550])

#Get row-wise root-N weights from average cross-dispersion profile of 20 central 
#columns:   
    numWeightCols = 20
    numWeightCols2 = int(numWeightCols/2)
    #default value
    weight0 = 1.0/maxWeight
    #initialize with default values
    colWeight = [weight0 for y in range(numWeightCols)]
    weight = [0.0 for y in range(deltaRows)]
    #pixel-wise weight
    weightRow = 0
    for row in range(centroid[refCols[1]]-halfDelta, centroid[refCols[1]]+halfDelta):
        #print("row ", row, " weightRow ", weightRow)
        weightCol = 0
        for col in range(refCols[1]-numWeightCols2, refCols[1]+numWeightCols2):
            if (data[row][col] < 1):
                colWeight[weightCol] = weight0
            else:
                colWeight[weightCol] = math.sqrt(data[row][col])
            weight[weightRow] += colWeight[weightCol]
            weightCol += 1
        weight[weightRow] = weight[weightRow] / maxWeight
        weight[weightRow] = weight[weightRow] / numWeightCols
        weightRow += 1
        
    plt.plot(cross[1], linewidth=2, color='blue')
    plt.plot(range(centroid[refCols[1]]-halfDelta, centroid[refCols[1]]+halfDelta), [x*4000 for x in weight], linewidth=2, color='red')
    plt.annotate(plotLabel, xy=(250, 100))
    plt.annotate("Mid-column X-dispersion profile", xy=(250, 3000), color='blue')
    plt.annotate("rootN weight", xy=(250, 2500), color='red')
    
    #Show zero level:
    plt.plot([-20, 450], [0.0, 0.0], '--', color='gray')
    
    #Show background subtraction regions:
        
    plt.plot([10, 10], [0, 500], color='black')
    plt.plot([50, 50], [0, 500], color='black')
    
    plt.plot([380, 380], [0, 500], color='black')
    plt.plot([420, 420], [0, 500], color='black')
    
    
    #Now produce the row-co-added spectrum!   
    spectrum = [0.0 for i in range(numCols)]
    for col in range(numCols):
        #for row in range(numRows):
        #row range being co-added slowly drifts with column number along chip
        #BUT 1D weight array is always the same:
        weightRow = 0
        for row in range(centroid[col]-halfDelta, centroid[col]+halfDelta):
            spectrum[col] = spectrum[col] + weight[weightRow]*data[row][col]
            weightRow += 1
            ##For testing:
            #spectrum[col] = spectrum[col] + data[row][col]
        
    #plt.plot(data[centroid])
    #plt.plot(spectrum)
    
    plt.show()

    #Save as encapsulated postscript (eps) for LaTex
    fileStem = "CrossDispersionWeight"+plotLabel
    epsName = fileStem + '.eps'
    plt.savefig(epsName, format='eps', dpi=1000)
    #Save as jpg for PowerPoint
    jpgName = fileStem + '.jpg'
    plt.savefig(jpgName, format='jpg')    
     
    
    #hdul.close()
    #return hdr, data
    return hdr, spectrum, cross, data, weight