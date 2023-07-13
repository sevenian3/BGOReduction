# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 14:06:39 2022

@author: Ian Short

Order of processing BGO spectra:
    
    1) Get the spectrum from the FITS file:
       spectrum = getFitsSpectrum()
       - automatically bias subtracted and flat-field corrected
    2) Calibrate the wavelength scale with ar-ne calib.ini data:
        wave = waveCalib(spectrum)
    3) Correct for instrumental response ONLY - no rectification:
        [wave, spectrum] = specInstrum(wave, spectrum)

@author: Ian Short
"""
import numpy as np
import matplotlib.pyplot as plt

def specInstrum(wave, spectrum):
    
    #Read in DL's "instrument response" - made by oberving an A0V star:
    dataPath = "C:/Users/s7970488/Documents/BGO/"
    #outPath = absPath + "Outputs/"

    numPix = len(spectrum)
    wavStr = ""
    flxStr = ""
    inLine = ""
    fields = [" " for i in range(2)]         
    numResp = 2398
    waveResp = [0.0 for i in range(numResp)]
    fluxResp = [0.0 for i in range(numResp)]
    
    inFile = dataPath + "A0V-instrument-response.dat"
    
    #inFile = dataPath + "A0V-smoothed.dat"
    with open(inFile, 'r') as inputHandle:
        
        #no header
        for i in range(numResp):
            inLine = inputHandle.readline()  
            fields = inLine.split()
            wavStr = fields[0].strip(); flxStr = fields[1].strip()
            waveResp[i] = float(wavStr); fluxResp[i] = float(flxStr)    
            
    #plt.plot(waveResp, fluxResp)
    
    #Data spectrum extends beyong range of instrumental response spectrum - must
    # be truncated before interpolation:
    
    print("numResp ", numResp, " waveResp[numResp-1] ", waveResp[numResp-1])
    print("wave[numPix-1] ", wave[numPix-1])
    for i in range(numPix):
        if (wave[i] >= waveResp[numResp-1]):
            break
        i+=1
    numT = i
    print("numPix ", numPix, " numT ", numT)
    
    waveT = np.array([x for x in wave[0:numT]])
    spectrumT = np.array([x for x in spectrum[0:numT]])
    
    #Interpolate instrument response onto 
    # Need wavelength calibration first:
    fluxResp2 = np.interp(waveT, waveResp, fluxResp)
    #plt.plot(waveT, fluxResp2)
    

    spectrumT = np.array(spectrumT)

    #Divide out instrumental response:
    
    spectrumInst = [spectrumT[x]/fluxResp2[x] for x in range(numT)]
            
    spectrumInst = np.array(spectrumInst)
    
    #plt.plot(wave, spectrum)
    #plt.plot(waveT, fluxResp2)
    #plt.plot(waveT, spectrumInst)

    #0th order normailzation:
    maxPix = spectrumInst.argmax()
    yMax = spectrumInst[maxPix]
    
    # Generate 0th order continuum spectrum:
    cont0 = [yMax for i in range(numT)]

    # 0th order rectification:
    spectrum0 = [spectrumInst[x]/cont0[x] for x in range(numT)]

    
    return waveT, spectrum0