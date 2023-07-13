# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 15:00:14 2022

Inputs:
    lambdaStart, lambdaStop:  (Angstroms) Restrict wavelength range to be fitted 
        with a continuum to this range

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
import matplotlib.pyplot as plt

def specNorm2(file, wave, spectrum, lambdaStart, lambdaStop):
    
    label = file[0:10]
    
    #Read in DL's "instrument response" - made by oberving an A0V star:
    dataPath = "C:/Users/s7970488/Documents/BGO/BGORed/"
    #outPath = absPath + "Outputs/"

    numPix = len(spectrum)
    wavStr = ""
    flxStr = ""
    inLine = ""
    fields = [" " for i in range(2)]         
    numResp = 2398
    waveResp = [0.0 for i in range(numResp)]
    fluxResp = [0.0 for i in range(numResp)]
    
    #inFile = dataPath + "A0V-smoothed.dat"
    inFile = dataPath + "A0V-instrument-response.dat"
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
    
    #print("numResp ", numResp, " waveResp[numResp-1] ", waveResp[numResp-1])
    #print("wave[numPix-1] ", wave[numPix-1])
    #print("wave[0:100] ", wave[0:100])
    for i in range(numPix):
        if (wave[i] >= waveResp[numResp-1]):
            break
        i+=1
    numT = i
    #print("numPix ", numPix, " numT ", numT)
    
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


    """
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
    #plt.plot(spectrum)
    #plt.plot(cont1)
    """
    
    ####  A trick to skip the 1st order corrections above:
        
    spectrum1 = np.array(spectrum0) ###  Comment this out if you want to include the 1st order correction first!


    #Normalize in TWO stages:
    # 1 - a low order fit to the MEAN flux in each lambda bin to remove global curvature and slope
    # 2 - a highe rorder fit to MAXIMUM flux in each bin to continuum rectify

    #Truncate spectrum to wavelength range to which we should try to fit the continuum.   
    if (lambdaStart < waveT[0]):
        lambdaStart = waveT[0]
    if (lambdaStop > waveT[numT-1]):
        lambdaStop = waveT[numT-1]
    
    #Find 1st index of truncated spectrum:
    for iStart in range(numT):
        if (waveT[iStart] >= lambdaStart):
            break
    #Find last index of truncated spectrum:
    for iStop in range(numT):
        if (waveT[iStop] >= lambdaStop):
            break
        iStop-=1
        
    waveTtrunc = waveT[iStart: iStop]
    spectrum1trunc = spectrum1[iStart: iStop]
    numTtrunc = iStop - iStart


        
    #######################
    
    ######## 1:
    
    ########################    
    
    # Fit a polynomial to binned MEAN flux values

    # Divide the spectrum into coarse bins and caluclate flux means in each bin
    # these are the ordinate "y" values we are fitting:
    numBins = 13 #200 A bins
    #Order of polynomial fit:
    order = 5
    
    #Avoid edge effects by omitting 'numEdge' pixels on either end:
    
    numEdge = 1 #initial guess
    numT2 = numTtrunc - (2*numEdge)
    numPerBin = int(round((numT2/numBins)-1.0))
    halfBin = int(round(numPerBin/2))
    
    #Fit to MEAN flux in each bin:
    yBin = [sum(spectrum1trunc[numEdge+(i*numPerBin):numEdge+numPerBin+(i*numPerBin)])/numPerBin for i in range(numBins)]
    ## Fit to MAXIMUM flux in each bin:
    #yBin = [max(spectrum1[numEdge+(i*numPerBin):numEdge+numPerBin+(i*numPerBin)]) for i in range(numBins)]
    
    #i = 0
    #print("1st bin: pixel range: ", [numEdge+(i*numPerBin), numEdge+numPerBin+(i*numPerBin)])
    #print("fluxes: ", spectrum1[numEdge+(i*numPerBin):numEdge+numPerBin+(i*numPerBin)])
    #print("Max (y) ", max(spectrum1[numEdge+(i*numPerBin):numEdge+numPerBin+(i*numPerBin)]))
   
    #Corresponding abscissae at bin mid-pixel - just for visualization:
    xBin = [numEdge+halfBin+(i*numPerBin) for i in range(numBins)]
    
    #print("halfBin ", halfBin)
    #print("xBin ", xBin, " waveT[0] ", waveT[0], " xBin[0] ", xBin[0], " waveT[x[0]] ", waveT[x[0]])
    #print("waveT[0:100] ", waveT[0:100])

    
    plt.figure()
    plt.subplot(1, 1, 1)    
    plt.title = ""
    plt.ylabel('$F_\lambda/F^C_\lambda$')
    plt.xlabel('$\lambda$ (A)')
    plt.ylim([0.0, 1.5])
    
    plt.annotate(label, xy=(4500, 0.1))
    
    plt.plot(waveT, spectrum1, color='green')
    
    plt.annotate("Step 1", xy=(6000, 0.1), color='red')
    plt.plot(waveTtrunc[xBin], yBin, 'o', color='red')
    
    
    #Solve for polynomial coefficients, vector c  
    c = np.polyfit(xBin, yBin, order)
    #print(c)
    
    
    # Generate nth order continuum spectrum:
    #contN = [ (c[0]*(x**4)) + (c[1]*(x**3)) + (c[2]*(x**2)) + (c[3]*x) + c[4] for x in range(numPix) ]
    contN = [0.0 for i in range(numTtrunc)]
    for x in range(numTtrunc):
        for n in range(order):
            polyI = 0.0
            exponent = order - n
            polyI = c[n]*(x**exponent)
            contN[x] = contN[x] + polyI
        contN[x] = contN[x] + c[order]
        

    # Nth order rectification:
    spectrumN = [spectrum1trunc[x]/contN[x] for x in range(numTtrunc)]

    
    plt.plot(waveTtrunc, spectrumN, color='blue', linewidth=2)
    plt.plot(waveTtrunc, contN, color='red')
    ##print("len(waveTtrunc) ", len(waveTtrunc), " numTtrunc ", numTtrunc)
    plt.plot(waveTtrunc, [1.0 for i in range(numTtrunc)], color='green')
    
    """
    plt.show()

    #Save as encapsulated postscript (eps) for LaTex
    fileStem = "SpecNorm2Step1"+label
    epsName = fileStem + '.eps'
    plt.savefig(epsName, format='eps', dpi=1000)
    #Save as jpg for PowerPoint
    jpgName = fileStem + '.jpg'
    plt.savefig(jpgName, format='jpg')
    """
    
    ############################
    
    # 2. 
    
    ############################
    
    # Fit a polynomial to binned MAXIMUM flux values

    # Divide the spectrum into coarse bins and caluclate flux means in each bin
    # these are the ordinate "y" values we are fitting:
    numBins2 = 13  # 200 A bins
    #Order of polynomial fit:
    order2 = 5
    
    #Avoid edge effects by omitting 'numEdge' pixels on either end:
    
    numEdge2 = 20 #initial guess
    numT2 = numTtrunc - (2*numEdge2)
    numPerBin2 = int(round((numT2/numBins2)-1.0))
    halfBin2 = int(round(numPerBin2/2))
    
    ##Fit to MEAN flux in each bin:
    #y2 = [sum(spectrumN[numEdge2+(i*numPerBin2):numEdge2+numPerBin2+(i*numPerBin2)])/numPerBin2 for i in range(numBins2)]
    # Fit to MAXIMUM flux in each bin:
    yBin2 = [max(spectrumN[numEdge2+(i*numPerBin2):numEdge2+numPerBin2+(i*numPerBin2)]) for i in range(numBins2)]
    
   
    #Corresponding abscissae at bin mid-pixel:
    xBin2 = [numEdge2+halfBin2+(i*numPerBin2) for i in range(numBins2)]
    
    #print("halfBin2 ", halfBin2)
    #print("x ", x, " waveT[0] ", waveT[0], " x[0] ", x[0], " waveT[x[0]] ", waveT[x[0]])
    #print("waveT[0:100] ", waveT[0:100])
    
    """
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.title = ""
    plt.ylabel('$F_\lambda/F^C_\lambda$')
    plt.xlabel('$\lambda$ (A)')
    plt.ylim([0.0, 1.5])
    
    plt.annotate(label, xy=(4500, 1.35))
    
    plt.annotate("Step 2", xy=(6000, 1.35), color='red')
    plt.plot(waveTtrunc, spectrumN, color='blue')  # From Step 1 above  
    plt.plot(waveTtrunc[xBin2], yBin2, 'o', color='red')
    """
    
    #Solve for polynomial coefficients, vector c  
    c2 = np.polyfit(xBin2, yBin2, order2)
    #print(c2)
    
    
    # Generate nth order continuum spectrum:
    #contN = [ (c[0]*(x**4)) + (c[1]*(x**3)) + (c[2]*(x**2)) + (c[3]*x) + c[4] for x in range(numPix) ]
    contN2 = [0.0 for x in range(numTtrunc)]
    for x in range(numTtrunc):
        for n in range(order2):
            polyI = 0.0
            exponent = order2 - n
            polyI = c2[n]*(x**exponent)
            contN2[x] = contN2[x] + polyI
        contN2[x] = contN2[x] + c2[order2]
        
    #Fix the "wings"
    for x in range(numTtrunc):
        if (waveTtrunc[x] < waveTtrunc[xBin2[0]]):
            contN2[x] = contN2[xBin2[0]]
        if (waveTtrunc[x] > waveTtrunc[xBin2[numBins2-1]]):
            contN2[x] = contN2[xBin2[numBins2-1]]
        

    # Nth order rectification:
    spectrumN2 = [spectrumN[x]/contN2[x] for x in range(numTtrunc)]

    
    plt.plot(waveTtrunc, spectrumN2, color='black', linewidth=2)
    plt.plot(waveTtrunc, contN2, color='red')
    plt.plot(waveTtrunc, [1.0 for i in range(numTtrunc)], color='green')
    
    """
    plt.show()

    #Save as encapsulated postscript (eps) for LaTex
    fileStem = "SpecNorm2Step2"+label
    epsName = fileStem + '.eps'
    plt.savefig(epsName, format='eps', dpi=1000)
    #Save as jpg for PowerPoint
    jpgName = fileStem + '.jpg'
    plt.savefig(jpgName, format='jpg')
    """
 
    return waveTtrunc, spectrumN2
        
    