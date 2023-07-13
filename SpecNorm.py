# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 13:35:00 2022

@author: s7970488
"""

import numpy as np
import matplotlib.pyplot as plt

def specNorm(spectrum):
    
    
    numPix = len(spectrum)
    spectrum = np.array(spectrum)

# 1st order normalization:
    
    minPix = spectrum.argmin()
    maxPix = spectrum.argmax()
    #Find mina nd max Y this way to Ensure consistenty:
    yMin = spectrum[minPix]
    yMax = spectrum[maxPix]

    deltaY = yMax - yMin
    deltaX = maxPix - minPix
    
    m = deltaY / deltaX # slope
    b = yMin - (m * minPix)
    
    # Generate 1st order continuum spectrum:
    cont1 = [ (m*x) + b for x in range(numPix)]

    # 1st order rectification:
    spectrum1 = [spectrum[x]/cont1[x] for x in range(numPix)]
    #spectrum1 = [0.0 for x in range(numPix)]
    #for i in range(numPix):
    #    spectrum1[i] = spectrum[i]/cont1[i]
    
    #plt.plot(spectrum)
    #plt.plot(cont1)
    #plt.plot(spectrum1)

    """    
    # 2nd order normalization:
    # Assume spectrum 1 is "smiling"
    # Fit a symmetric offset parabola:  y = c2*(DeltaX)^2 + C0
        
    spectrum1 = np.array(spectrum1)
        
    DeltaX = spectrum1.argmin()
    xPrime = [x-DeltaX for x in range(numPix)]
    c0 = spectrum1[DeltaX]  #minimum y value
    
    maxPix = spectrum1.argmax()
    maxPixPrime = maxPix - DeltaX
    yMax = spectrum1[maxPix]
    
    c2 = (yMax - c0) / (maxPixPrime**2)
    
    # Generate 2nd order continuum spectrum:   
    cont2 = [ (c2*(xPrime[x]**2)) + c0 for x in range(numPix)]
    
    # 2nd order rectification:
    spectrum2 = [spectrum1[x]/cont2[x] for x in range(numPix)]

    
    plt.plot(spectrum1)
    #plt.plot(cont2)
    plt.plot(spectrum2)   
    """
    
    """
    # 4th order normalization:
    #Assume spectrum2 is 'frowning'
    # Fit a symmetric offset quartic:  y = c4*(DeltaX)^4 + C0    
        
    #spectrum2= np.array(spectrum2)
    #Trick to jump from 1st order to 4th order normalization. skipping 2nd order
    # IF 2nd order block commented out above
    spectrum2 = np.array(spectrum1)
        
    DeltaX = spectrum2.argmax()
    xPrime = [x-DeltaX for x in range(numPix)]
    c0 = spectrum2[DeltaX]   #maximum y value
    
    minPix = spectrum2.argmin()
    minPixPrime = minPix - DeltaX
    yMin = spectrum2[minPix]
    
    c4 = (yMin - c0) / (minPixPrime**4)
    
    # Generate 4th order continuum spectrum:   
    cont4 = [ (c4*(xPrime[x]**4)) + c0 for x in range(numPix)]
    
    # 4th order rectification:
    spectrum4 = [spectrum2[x]/cont4[x] for x in range(numPix)]

    
    plt.plot(spectrum2)
    plt.plot(cont4)
    plt.plot(spectrum4) 
    """
    
    # 4th order normalization:
    #Assume spectrum2 is "smiling"
    # Fit a symmetric offset quartic:  y = c4*(DeltaX)^4 + C0    
        
    #spectrum2= np.array(spectrum2)
    #Trick to jump from 1st order to 4th order normalization. skipping 2nd order
    # IF 2nd order block commented out above
    spectrum2 = np.array(spectrum1)
        
    DeltaX = spectrum2.argmin()
    xPrime = [x-DeltaX for x in range(numPix)]
    c0 = spectrum2[DeltaX]   #maximum y value
    
    maxPix = spectrum2.argmax()
    maxPixPrime = maxPix - DeltaX
    yMax = spectrum2[maxPix]
    
    c4 = (yMax - c0) / (maxPixPrime**4)
    
    # Generate 4th order continuum spectrum:   
    cont4 = [ (c4*(xPrime[x]**4)) + c0 for x in range(numPix)]
    
    # 4th order rectification:
    spectrum4 = [spectrum2[x]/cont4[x] for x in range(numPix)]

    
    plt.plot(spectrum2)
    plt.plot(cont4)
    plt.plot(spectrum4)   
    
    
    
    
    
    
    
    
        
    
        
    
    