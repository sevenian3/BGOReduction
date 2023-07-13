# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 11:23:12 2022

Order of processing BGO spectra:
    
    1) Get the spectrum from the FITS file:
       spectrum = getFitsSpectrum()
       - automatically bias subtracted and flat-field corrected
    2) Calibrate the wavelength scale with ar-ne calib.ini data:
        wave = waveCalib(spectrum)
    3) Correct for instrumental response and rectify:
        [wave, spectrum] = specNorm2(wave, spectrum)

Input:  
    spectrum: needed only to get the number of pixels in dispersion direction of 
    spectrum being calibrated
    
@author: Ian Short
"""
import numpy as np

def waveCalib(spectrum):
    
    #Read in DL's "ar-ne calib.ini" - made by oberving an Ar-Ne lamp?:
    #Plain text ascii
    dataPath = "C:/Users/s7970488/Documents/BGO/BGORed/"
    #outPath = absPath + "Outputs/"

    numDisp = len(spectrum)
    #print("numDisp ", numDisp)
    wavStr = ""
    pixStr = ""
    inLine = ""
    fields = [" " for i in range(2)]         
    numCalibPoints = 9  #hard-wire for now
    
    wav = [0.0 for i in range(numCalibPoints)]
    pix = [0.0 for i in range(numCalibPoints)]
    
    inFile = dataPath + "ar-ne calib.ini"
    with open(inFile, 'r') as inputHandle:
        
        #one line of header:
        inLine = inputHandle.readline()
        #print(inLine)
        for i in range(numCalibPoints):
            #Two lines per calibration point:
            inLine = inputHandle.readline()  
            fields = inLine.split("=")
            wavStr = fields[1].strip()
            inLine = inputHandle.readline()
            fields = inLine.split("=")
            pixStr = fields[1].strip()
            wav[i] = float(wavStr); pix[i] = float(pixStr)
 
    #pix2 = [x for x in range(numDisp)]
    #WRONG! #Interpolate wav vs pix onto the number of pixels along dispersion axis of detector:
    #    
    #wave2 = np.interp(pix2, pix, wav)
    
    #Fit parabola to lambda vs pixel relation:
    c = np.polyfit(pix, wav, 2)
    #print("Wavelength calibration coeffs: c[0] ", c[0], " c[1] ", c[1], " c[2] ", c[2])
    
    wave2 = [(c[0]*(x**2))+(c[1]*x)+c[2] for x in range(numDisp)]
       
    return wave2