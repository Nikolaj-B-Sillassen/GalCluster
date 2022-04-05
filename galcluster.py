# -*- coding: utf-8 -*-
"""
# Automated Identification of Galaxy Clusters
# galcluster v0.1
@author: Nikolaj Bjerregaard Sillassen
"""
from astropy.io import fits
from astropy.modeling import models, fitting
import numpy as np
import time
import multiprocessing as mp
import sys
from sklearn import cluster,neighbors


def sigmaCalcFun(data):
    """
    Function used to calculate sigmas for galcluster
    
    Parameters
    ----------
    data : ndarray
        ndarray of objects containing: [RA,DEC,square_size,sigmas]
        RA : ndarray
            ndarray of length N containing right ascension for all sources
        DEC : ndarray
            ndarray of length N containing declination for all sources
        square_size : float
            size of the moving window in which to calculate distances 
        sigmas : list
            list of length M containing sigmas used for sigma calculation

    Returns
    -------
    sigmaArray : ndarray
        ndarray of size M by N containing calculated sigma values for all 
        sigmas and all sources

    """
    RA = data[0]
    DEC = data[1]
    square_size = data[2]
    sigmas = []
    for i in range(len(data)-3):
        sigmas.append(data[3+i])
    i = 0
    sigmaArray = np.zeros([len(sigmas[0]), len(RA)])
    total_time = 0
    while i <= len(RA)-1:
        start = time.time()
        [squareCenterRA, squareCenterDEC] = [
            RA[i], DEC[i]]  # Find center of square window
        [lowRAID, highRAID, lowDECID, highDECID] = [RA >= squareCenterRA - square_size, RA <= squareCenterRA +
                                                    square_size, DEC >= squareCenterDEC - square_size, DEC <= squareCenterDEC + square_size]
        SourcesInSquare = lowRAID & highRAID & lowDECID & highDECID
        distances = np.sqrt((RA[SourcesInSquare]-squareCenterRA)
                            ** 2+(DEC[SourcesInSquare]-squareCenterDEC)**2)
        sortedDistances = np.sort(distances)
        for j in range(len(sigmas[0])):
            sigma = sigmas[0][j]
            if sum(SourcesInSquare) >= sigma+1:
                sigmaArray[j][i] = sigma/(np.pi*sortedDistances[int(sigma)]**2)
            else:
                sigmaArray[j][i] = float('NaN')
        stop = time.time()
        total_time += stop-start
        if i % 100 == 0:
            print("Progress: {0:.2f}%".format(i/(len(RA)-1)*100))
            print("Elapsed time: {0:.1f}s".format(total_time))
        i += 1
    return sigmaArray


def clusterFun(RA, DEC, sigma, SNsigma, SNlimit, min_c_size, ODSigma, Z):
    """
    Function used to generate overdensity regions and cluster candidates
    
    Parameters
    ----------
    RA : ndarray
        ndarray of length N containing right ascension for all sources
    DEC : ndarray
        ndarray of length N containing declination for all sources
    sigma : ndarray
        ndarray of length N containing sigma values for all sources, as 
        calculated by sigmaCalcFun
    SNsigma : ndarray
        ndarray of length N containing signal to noise of sigma values, for 
        all sources
    SNlimit : float
        cutoff in signal to noise of sigma, for a source to be significantly
        overdense
    min_c_size : int
        minimum number of sources in an overdensity region to be considered
        a cluster candidate
    ODSigma : int
        value of sigma used in sigma calculations
    Z : ndarray (optional)
        ndarray of length N, containing redshift for all sources, if redshift
        is not known input an array of length N containing nans

    Returns
    -------
    clusterCands : list
        list containing all cluster candidates, will be empty if redshift is
        not passed
        each cluster candidate is a list of the following values:
        [RA,DEC,z,Sigma,SignalToNoiseSigma] for all sources in the
        candidates
    ODRs : list
        list containing all overdensity regions
        each overdensity region is a list of the following values:
        [RA,DEC,z,Sigma,SignalToNoiseSigma] for all sources in the
        regions

    """
    print("Clustering overdensities in Sigma{0:d}".format(ODSigma))
    # Clustering of sigma overdensities
    largestSourceIDs = SNsigma >= SNlimit
    largestSources = SNsigma[largestSourceIDs]
    largestOD = SNsigma[largestSourceIDs]
    largestSigma = sigma[largestSourceIDs]
    LRA, LDEC = RA[largestSourceIDs], DEC[largestSourceIDs]
    LZ = Z[largestSourceIDs]
    RAsortedIDs = np.argsort(LRA)
    sLRA, sLDEC, sLZ, ssigma, sSNS = LRA[RAsortedIDs], LDEC[RAsortedIDs], LZ[
        RAsortedIDs], largestSigma[RAsortedIDs], largestOD[RAsortedIDs]
    # & (abs(sLDEC5[:-1]-sLDEC5[1:])>=0.01)# & (abs(sLZ5[:-1]-sLZ5[1:])>=0.1)
    cluster = (abs(sLRA[:-1]-sLRA[1:]) >= ODSigma/500)
    clusterLim = np.where(cluster == True)
    N_clustersRA = sum(cluster)+1
    RAClusters = np.empty(N_clustersRA, dtype=object)
    cRA, cDEC, cZ, cS, cSNS = np.split(sLRA, clusterLim[0]+1), np.split(sLDEC, clusterLim[0]+1), np.split(
        sLZ, clusterLim[0]+1), np.split(ssigma, clusterLim[0]+1), np.split(sSNS, clusterLim[0]+1)

    tc = list()
    od = list()
    # Cluster Determination
    for i in range(N_clustersRA):
        RAClusters[i] = np.array([cRA[i], cDEC[i], cZ[i], cS[i], cSNS[i]])
        sDecS = np.argsort(RAClusters[i][1])
        RAClusters[i] = np.array(
            [cRA[i][sDecS], cDEC[i][sDecS], cZ[i][sDecS], cS[i][sDecS], cSNS[i][sDecS]])
        clusterDec = abs(RAClusters[i][1][:-1] -
                         RAClusters[i][1][1:]) >= ODSigma/1000
        N_clustersDEC = sum(clusterDec)+1
        if N_clustersDEC != 1:
            climDec = np.where(clusterDec == True)
            tcRA, tcDEC, tcZ, tcS, tcSN = np.split(RAClusters[i][0], climDec[0]+1), np.split(RAClusters[i][1], climDec[0]+1), np.split(
                RAClusters[i][2], climDec[0]+1), np.split(RAClusters[i][3], climDec[0]+1), np.split(RAClusters[i][4], climDec[0]+1)
            for o in range(len(tcRA)):
                if len(tcRA[o]) >= min_c_size:
                    od.append(
                        np.array([tcRA[o], tcDEC[o], tcZ[o], tcS[o], tcSN[o]], dtype=object))
        else:
            tcRA, tcDEC, tcZ, tcS, tcSN = RAClusters[i][0], RAClusters[i][
                1], RAClusters[i][2], RAClusters[i][3], RAClusters[i][4]
            # if len(tcRA) >= min_c_size:
            od.append(np.array([tcRA, tcDEC, tcZ, tcS, tcSN], dtype=object))
        tempC = np.empty(N_clustersDEC, dtype=object)
        countdec = 0
        for j in range(N_clustersDEC):
            if len(clusterDec) == 0:
                countdec += 1
                countz = 0
                break
            else:
                try:
                    if N_clustersDEC != 1:
                        tempC[countdec] = np.array(
                            [tcRA[countdec], tcDEC[countdec], tcZ[countdec], tcS[countdec], tcSN[countdec]])
                    else:
                        tempC[countdec] = np.array(
                            [tcRA, tcDEC, tcZ, tcS, tcSN])
                    if len(tempC[countdec][0]) > 1:
                        sZS = np.argsort(tempC[countdec][2])
                        tempC[countdec] = np.array([tempC[countdec][0][sZS], tempC[countdec][1][sZS],
                                                   tempC[countdec][2][sZS], tempC[countdec][3][sZS], tempC[countdec][4][sZS]])
                        clusterz = abs(
                            (tempC[countdec][2][:-1]-tempC[countdec][2][1:])) >= 0.1*(1+tempC[countdec][2][:-1])
                        N_clustersz = sum(clusterz)+1
                    if N_clustersz != 1:
                        climz = np.where(clusterz == True)
                        t2cRA, t2cDEC, t2cZ, t2cS, t2cSN = np.split(tempC[countdec][0], climz[0]+1), np.split(tempC[countdec][1], climz[0]+1), np.split(
                            tempC[countdec][2], climz[0]+1), np.split(tempC[countdec][3], climz[0]+1), np.split(tempC[countdec][4], climz[0]+1)
                    else:
                        t2cRA, t2cDEC, t2cZ, t2cS, t2cSN = tempC[countdec][0], tempC[countdec][
                            1], tempC[countdec][2], tempC[countdec][3], tempC[countdec][4]
                    countdec += 1
                    countz = 0
                except:
                    N_clustersz = 0
                    countdec += 1
                    countz = 0
                    pass

            countz = 0
            for k in range(N_clustersz):
                if len(clusterz) == 0:
                    countz += 1
                    break
                if N_clustersz != 1:
                    if t2cRA[countz].size >= min_c_size:
                        tc.append(np.array(
                            [t2cRA[countz], t2cDEC[countz], t2cZ[countz], t2cS[countz], t2cSN[countz]]))
                        countz += 1
                    else:
                        countz += 1
                else:
                    if t2cRA.size >= min_c_size:
                        tc.append(np.array([t2cRA, t2cDEC, t2cZ, t2cS, t2cSN]))
                        countz += 1
                    else:
                        countz += 1
    clusterCands = tc
    ODRs = od
    return clusterCands, ODRs


def galcluster(filename=None, RAKey="RA", DECKey="DEC", method="sigma", sigmas=[5, 10], multi_processing=True, mute_plots=False, SNlimit=[5, 7], min_c_size=3, zKey=None):
    """
    Automated Identification of Galaxy Clusters
    galcluster v0.1

    Program created during Synthesis Project at DTU Space, in spring 2022

    Function description:

    galcluster(filename=None, RAKey="RA", DECKey="DEC", method="sigma", sigmas=[5,10], multi_processing=True, mute_plots=False, SNlimit=[5,7],min_c_size=3,zKey=None)

    ### Input Parameters

      filename : str, optional
          Name of the file to run clustering algorithm on, must be .fits format 
      
      RAKey : str, optional
          Keyword of the column containing right ascension
        
      DECKey : str, optional
          Keyword of the column containing declination

      method : str, optional 
          "sigma" (default) will use sigma5 and sigma10 overdensities,
          "dbscan" will use Density Based Spacial Clustering for Applications with Noise, to produce galaxy cluster candidates

      sigmas : list, optional
          list of sigmas, to use for overdensity calculations, default is 
          [5,10]

      multi_processing : bool, optional 
          If true (default) enables parallelization using multiprocessing, 
          if not no parallelization will be used, which is much slower 
          must be run under if __name__ == "__main__": protection when using 
          multiprocessing 

      mute_plots : bool, optional 
          If false (default), plots are created, otherwise plots are muted. 

      SN_limit : list, optional
          limit the signal to noise of overdensity for clustering for each SigmaN, [5, 7] is default

      min_c_size : int, optional 
          minumum number of sources in a cluster candidate, 3 is default

      zKey : str, optional
          keyword for column containing redshift

      z_cut : positive float, optional 
          minimum redshift of the sources that are being clustered, 1.5 is 
          default. 

     ### Returns

      out : [RA,DEC,Z,Sigmas,SignalToNoiseSigmas], clusterList, ODRlist

      First list:
          RA : ndarray 
              array containing right ascension for all sources
          DEC : ndarray
              array containing declination for all sources
          Z : ndarray
              array containing redshifts for all sources, not present if 
              redshift keyword is not given
          Sigmas : ndarray
              array containing all sigma values, for all pre-selected sources
          SignalToNoiseSigmas : ndarray
              array containing all signal to noise values of all sigmas, for 
              all pre-selected sources
      clusterList : list
          List is not present if keyword for redshift is not given
          First Dimension : list
              List of cluster candidate lists, for all values of input sigmas
          Second Dimension : list
              [RA,DEC,Z,Sigma,SignalToNoiseSigma] for all sources in the
              candidates
      ODRlist : list
          First Dimension : list
              List of overdensity region lists, for all input sigmas
          Second Dimension : list
              [RA,DEC,Z,Sigma,SignalToNoiseSigma] for all sources in the
              overdensity regions
    """
    while True:
        # Setup all parameters
        from astropy.io import fits
        if filename == None:
            while True:
                filename = input("Please enter the filename: ")
                try:
                    dataTable = fits.open(filename)
                    dataset = dataTable[1].data
                except:
                    sys.exit("Error: {0} not found".format(filename))
                RAKey = input("Please enter the keyword for right ascension: ")
                try:
                    test = dataset[RAKey]
                except:
                    sys.exit("Error: {0} not found in {1}".format(RAKey,filename))
                DECKey = input("Please enter the keyword for declination: ")
                try:
                    test = dataset[DECKey]
                except:
                    sys.exit("Error: {0} not found in {1}".format(DECKey,filename))
                zpresel = input("Does the catalog contain redshift? (y/n): ")
                if zpresel.lower() == "y" or zpresel.lower() == "yes":
                    zKey = input(
                        "Please enter the keyword for redshift: ")
                    try:
                        test = dataset[zKey]
                    except:
                        sys.exit("Error: key {0} not found in file {1}".format(
                            zKey, filename))
                defaultSettings = input(
                    "Do you wish to use default settings? (y/n): ")
                if defaultSettings.lower() == "n" or defaultSettings.lower() == "no":
                    typ = input(
                        "Please select the clustering method; 1 for sigma, 2 for DBSCAN: ")
                    if int(typ) == 1 or typ.lower() == "sigma":
                        method = "sigma"
                    elif int(typ) == 2 or typ.lower() == "dbscan":
                        method = "dbscan"
                    else:
                        sys.exit("Error: not a valid selection")
                    if method == "sigma":
                        Nsigmas = input(
                            "How many sigma values would you like to calculate overdensities for? ")
                        Nsigmas = int(Nsigmas)
                        sigmas = []
                        SNlimit = []
                        for i in range(Nsigmas):
                            sigm = input(
                                "Please input sigma {0}: ".format(i+1))
                            SNsigm = input(
                                "Please input the signal to noise minimum for overdensities in Sigma{}".format(int(sigm)))
                            sigmas.append(int(sigm))
                            SNlimit.append(float(SNsigm))
                        mult = input(
                            "Do you wish to use multiprocessing? (y/n): ")
                        if mult.lower() == "y" or mult.lower() == "yes":
                            multi_processing = True
                        elif mult.lower() == "n" or mult.lower() == "no":
                            multi_processing = False
                        else:
                            sys.exit("Error: invalid input")
                    plts = input("Do you wish to produce plots? (y/n): ")
                    if plts.lower() == "y" or plts.lower() == "yes":
                        mute_plots = False
                    elif plts.lower() == "n" or plts.lower() == "no":
                        mute_plots = True
                    else:
                        sys.exit("Error: invalid input")
                    try:
                        min_c_size = int(
                            input("Please enter the minimum size of a cluster for candidate detection: "))
                        break
                    except:
                        sys.exit("Error: not a number")
                else:
                    method = "sigma"
                    sigmas = [5, 10]
                    multi_processing = True
                    min_c_size = 3
                    SNlimit = [5, 7]
                    mute_plots = False
                    break
        import numpy as np
        dataTable = fits.open(filename)
        dataset = dataTable[1].data
        N_sources = len(dataset)
        print("Number of sources: {}".format(N_sources))
        if N_sources >= 100000:
            print("Warning: Dataset very Large")
        RA = dataset[RAKey]
        DEC = dataset[DECKey]
        if zKey != None:
            Z = dataset[zKey]
        else:
            Z = np.nan*np.ones(N_sources)
        avg_density = N_sources/((max(RA)-min(RA))*(max(DEC)-min(DEC)))
        start1 = time.time()

        square_size = 5*np.sqrt(20/avg_density)
        if method == "sigma":
            sigmaList = np.empty(len(sigmas), dtype=object)
            if multi_processing:
                if N_sources/mp.cpu_count() > 10000:
                    N_iters = round(N_sources/10000)
                    sigmaTotal = np.empty(len(sigmas), dtype=object)
                    RASplittemp, DECSplittemp = np.array_split(
                        RA, N_iters), np.array_split(DEC, N_iters)
                    print("Starting overdensity calculation")
                    for j in range(N_iters):
                        RASplit, DECSplit = np.array_split(
                            RASplittemp[j], mp.cpu_count()), np.array_split(DECSplittemp[j], mp.cpu_count())
                        sigmastruct = list(np.asarray(
                            [sigma*np.ones(mp.cpu_count()) for sigma in sigmas]).T)
                        dataStruct = np.array(np.transpose(
                            [RASplit, DECSplit, square_size*np.ones(mp.cpu_count()), sigmastruct]), dtype=object)
                        dataStruct = np.hstack([dataStruct, sigmastruct])
                        with mp.Pool(mp.cpu_count()) as p:
                            sigmastemp = p.map(sigmaCalcFun, dataStruct)
                            for k in range(len(sigmastemp)):
                                if k == 0:
                                    for p in range(len(sigmas)):
                                        sigmaTotal[p] = np.array(sigmastemp[0][p],dtype=float)
                                else:
                                    for p in range(len(sigmas)):
                                        sigmaTotal[p] = np.concatenate(
                                            (sigmaTotal[p], sigmastemp[k][p]))
                            stop = time.time()
                        if j == 0:
                            for k in range(len(sigmas)):
                                sigmaList[k] = sigmaTotal[k]
                        else:
                            for k in range(len(sigmas)):
                                sigmaList[k] = np.concatenate((sigmaList[k],sigmaTotal[k]))
                            #sigmaList = np.vstack((sigmaList, sigmaTotal))
                        print("Progress: {0:.2f}%, time: {1:.2f}s".format(
                            (j+1)/(N_sources/10000)*100, stop-start1))

                else:
                    RASplit, DECSplit = np.array_split(
                        RA, mp.cpu_count()), np.array_split(DEC, mp.cpu_count())
                    sigmastruct = list(np.asarray(
                        [sigma*np.ones(mp.cpu_count()) for sigma in sigmas]).T)
                    dataStruct = np.array(np.transpose(
                        [RASplit, DECSplit, square_size*np.ones(mp.cpu_count()), sigmastruct]), dtype=object)
                    print("Starting overdensity calculations")
                    with mp.Pool(mp.cpu_count()) as p:
                        sigmastemp = p.map(sigmaCalcFun, dataStruct)
                        for j in range(len(sigmastemp)):
                            if j == 0:
                                for k in range(len(sigmas)):
                                    sigmaList[k] = sigmastemp[0][k]
                            else:
                                for k in range(len(sigmas)):
                                    sigmaList[k] = np.concatenate(
                                        (sigmaList[k], sigmastemp[j][k]))
                        stop = time.time()
                        print("Done calculating densities for Sigma")
                        print("Total calculation time: {0:.2f}s".format(
                            stop-start1))

            else:
                dataStruct = [RA, DEC, square_size]
                sigmasTemp = sigmaCalcFun(dataStruct)
                sigmaList.append(sigmasTemp)
                stop = time.time()
                print("Calculation time: {0:.2f}s".format(stop-start1))

            # Clustering of sigma overdensities
            histList = []
            binCenterList = []
            meanSigmaList = []
            stdSigmaList = []
            SignalToNoiseList = []
            clusterList = []
            ODRlist = []
            gList = []
            LRAlist, LDEClist, LSNList = [], [], []
            LSIDlist = []
            z_cut = min(Z)
            for i in range(len(sigmaList)):
                sigmahist = np.histogram(sigmaList[i], np.logspace(
                    min(np.log10(sigmaList[i])), max(np.log10(sigmaList[i])), 25))
                binCenters = np.mean(
                    np.vstack([sigmahist[1][0:-1], sigmahist[1][1:]]), axis=0)
                g_init = models.Gaussian1D(amplitude=np.nanmax(sigmahist[0]), mean=np.nanmean(
                    np.log10(sigmaList[i])), stddev=np.nanstd(np.log10(sigmaList[i])))
                fit_g = fitting.LevMarLSQFitter()
                g = fit_g(g_init, np.log10(binCenters), sigmahist[0])
                gList.append(g)
                sigmamean, sigmastd = g.mean.value, g.stddev.value
                SignalToNoiseSigma = (
                    np.log10(sigmaList[i])-sigmamean)/sigmastd
                histList.append(sigmahist)
                binCenterList.append(binCenters)
                meanSigmaList.append(sigmamean)
                stdSigmaList.append(sigmastd)
                SignalToNoiseList.append(SignalToNoiseSigma)
                tc, od = clusterFun(
                    RA, DEC, sigmaList[i], SignalToNoiseSigma, SNlimit[i], min_c_size, sigmas[i], Z=Z)
                clusterList.append(tc)
                ODRlist.append(od)
                largestSourceID = SignalToNoiseSigma >= SNlimit[i]
                LRAlist.append(RA[largestSourceID])
                LDEClist.append(DEC[largestSourceID])
                LSNList.append(SignalToNoiseSigma[largestSourceID])
                LSIDlist.append(largestSourceID)

            # Plotting and Visualisation
            if mute_plots == False:
                for i in range(len(sigmas)):
                    import matplotlib.pyplot as plt
                    import matplotlib
                    import random
                    from scipy import interpolate

                    gist_heat_r = matplotlib.cm.get_cmap('gist_heat_r', 256)
                    new_cmp = matplotlib.colors.ListedColormap(np.vstack([gist_heat_r(
                        np.linspace(0, 10**(-.75), 50)), gist_heat_r(np.logspace(-.75, 0, 206))]))

                    def color_generator(n): return list(
                        map(lambda i: "#" + "%06x" % random.randint(0x333333, 0x555555), range(n)))
                    colors5 = color_generator(len(clusterList[i]))

                    meshRA, meshDEC = np.meshgrid(np.linspace(min(RA), max(RA), 2*round(
                        np.sqrt(len(RA)))), np.linspace(min(DEC), max(DEC), 2*round(np.sqrt(len(DEC)))))
                    interp = interpolate.NearestNDInterpolator(
                        list(zip(RA, DEC)), SignalToNoiseList[i])
                    intSN = interp(meshRA, meshDEC)
                    
                    plt.figure()
                    ax = plt.axes()
                    ax.set_facecolor([1, 1, 1])
                    ax.set_xlim([max(RA), min(RA)])
                    ax.set_ylim([min(DEC), max(DEC)])
                    M = ax.transData.get_matrix()
                    xscale = -M[0, 0]
                    yscale = M[1, 1]
                    maxSN = np.floor(np.nanmax(SignalToNoiseList[i]))
                    plt.pcolormesh(meshRA, meshDEC, intSN, cmap=new_cmp,
                                   shading="gouraud", vmin=-2, vmax=maxSN)  # ,antialiased=True)
                    cbar = plt.colorbar()

                    for j in range(len(ODRlist[i])):
                        if j == 0:
                            plt.scatter(np.mean(ODRlist[i][j][0]), np.mean(ODRlist[i][j][1]), s=xscale*(np.sqrt(sigmas[i]/(np.mean(
                                ODRlist[i][j][3])*np.pi)))*10, label="S/N $\geq$ {0}".format(SNlimit[i]), edgecolors='b', facecolors="none")
                        else:
                            plt.scatter(np.mean(ODRlist[i][j][0]), np.mean(ODRlist[i][j][1]), s=xscale*(np.sqrt(
                                sigmas[i]/(np.mean(ODRlist[i][j][3])*np.pi)))*10, edgecolors='b', facecolors="none")
                    plt.legend()
                    if zKey != None:
                        if i == 0:
                            ccandp = input(
                                "Do you wish to plot cluster candidates? (y/n): ")
                        if ccandp.lower() == "y" or ccandp.lower() == "yes":
                            for j in range(len(clusterList[i])):
                                plt.scatter(clusterList[i][j][0], clusterList[i][j][1], 40, edgecolors=colors5[j], facecolors="none")
                                plt.annotate("#{0:d}, z~ = {1:.2f}".format(j+1, np.median(clusterList[i][j][2])), [
                                             np.amax(clusterList[i][j][0])+0.1*(max(RA)-min(RA)), np.amin(clusterList[i][j][1])-0.05*(max(DEC)-min(DEC))], color=colors5[j])
                    plt.xlabel("ra [deg]", fontsize=14)
                    plt.ylabel("dec [deg]", fontsize=14)
                    cbar.ax.get_yaxis().labelpad = 15
                    cbar.ax.set_ylabel("S/N", rotation=90, fontsize=14)
                    plt.title("Signal to Noise ratio of $\Sigma_{%d}$" % (
                        sigmas[i]), fontsize=18)
                    plt.figure()
                    plt.plot(np.log10(
                        binCenterList[i]), histList[i][0], 'r', label="$\Sigma_{%d}$ histogram" % (sigmas[i]))
                    plt.plot(np.linspace(min(np.log10(binCenterList[i])), max(np.log10(binCenterList[i])), 100), gList[i](np.linspace(min(
                        np.log10(binCenterList[i])), max(np.log10(binCenterList[i])), 100)), 'k', label="$\Sigma_{%d}$ gaussian fit" % (sigmas[i]))
                    plt.legend(loc="upper left")
                    sSigma = np.sort(sigmaList[i][LSIDlist[i]])
                    sSigma = sSigma[::-1]
                    try:
                        for j in range(10):
                            plt.arrow(np.log10(sSigma[j]), max(histList[i][0]), 0, -.1*max(
                                histList[i][0]), width=.01, fill=0, head_length=.01*max(histList[i][0]))
                    except:
                        pass
                    try:
                        plt.annotate("Largest overdensities", [
                                     np.log10(sSigma[-1])-.5, max(histList[i][0])+5])
                    except:
                        pass
                    plt.title(
                        "Histogram of $\mathrm{log}_{10}\Sigma_{%d}$" % (sigmas[i]), fontsize=18)
                    plt.xlabel(
                        "$\mathrm{log}_{10}\Sigma_{%d}/\mathrm{deg}^{-2}$" % (sigmas[i]), fontsize=14)
                    plt.ylabel("N", fontsize=14)
                    if zKey != None:
                        if i == 0:
                            zPhotHistPlot = input(
                                "Do you wish to plot histograms of {0} for all overdensity regions with more than {1} sources? (y/n): ".format(zKey,min_c_size))
                        if zPhotHistPlot.lower() == "y" or zPhotHistPlot.lower() == "yes":
                            # z_phot histograms for all significant overdensity regions
                            for j in range(len(ODRlist[i])):
                                if len(ODRlist[i][j][2]) >= min_c_size:
                                    plt.figure()
                                    hist = np.histogram(ODRlist[i][j][2], bins=np.linspace(
                                        z_cut, 6, round((6-z_cut)/0.2)))
                                    binCenters = np.mean(
                                        np.vstack([hist[1][0:-1], hist[1][1:]]), axis=0)
                                    plt.plot(binCenters, hist[0], 'o-')
                                    plt.xlabel("{}".format(zKey),fontsize=14)
                                    plt.ylabel("N",fontsize=14)
                                    plt.title("ODR Histogram of redshift Sigma{0:d}, RA: {1:.4f} DEC: {2:.4f}".format(
                                        sigmas[i], np.mean(ODRlist[i][j][0]), np.mean(ODRlist[i][j][1])))
                    plt.show()
            save = input(
                "Do you wish to save the output in a .fits file? (y/n): ")
            if save.lower() == "y" or save.lower() == "yes":
                saveName = input(
                    "Input name of file you would like to save to: ")
                from astropy.io import fits
                import astropy.table as tb
                t = tb.Table([RA, DEC, Z], names=["RA", "DEC", "z"])
                for i in range(len(sigmas)):
                    sigmatemp = tb.Column(
                        sigmaList[i], name="Sigma{0}".format(sigmas[i]))
                    SNsigmatemp = tb.Column(
                        SignalToNoiseList[i], name="SNR_Sigma{0}".format(sigmas[i]))
                    t = tb.hstack([t, sigmatemp])
                    t = tb.hstack([t, SNsigmatemp])
                t.write(saveName, format='fits')
            saveODR = input(
                "Do you wish to save the overdensity regions in a .fits file? (y/n): ")
            if saveODR.lower() == "y" or saveODR.lower() == "yes":
                saveODRName = input(
                    "Input name of file you would like to save overdensity regions to: ")
                from astropy.io import fits
                import astropy.table as tb
                t = tb.Table()
                prevN = 1
                for i in range(len(sigmas)):
                    ttemp = tb.Table()
                    for j in range(len(ODRlist[i])):
                        tabletemp = tb.Table(data=[np.ones(len(ODRlist[i][j][0]))*sigmas[i], np.ones(len(ODRlist[i][j][0]))*prevN, ODRlist[i][j][0].astype(float), ODRlist[i][j][1].astype(
                            float), ODRlist[i][j][3].astype(float), ODRlist[i][j][4].astype(float), ODRlist[i][j][2].astype(float)], names=("SigmaN", "ODR_num", "RA", "DEC", "Sigma", "SNR_Sigma", "z_phot"))
                        ttemp = tb.vstack([ttemp, tabletemp])
                        prevN += 1
                    t = tb.vstack([t, ttemp])
                t.write(saveODRName, format='fits')
            saveDS9 = input(
                "Do you wish to save the candidates in a DS9 region .reg file? (y/n): ")
            if saveDS9.lower() == "y" or saveDS9.lower() == "yes":
                from astropy.coordinates import SkyCoord,FK5,Longitude,Latitude
                import astropy.units as u
                import regions
                saveDS9name = input("Please input the filename for the DS9 region file: ")
                regionsList = regions.Regions([])
                count = 1
                for j in range(len(sigmas)):
                    if zKey != None:
                        for i in range(len(clusterList[j])):
                            for k in range(len(clusterList[j][i][0])):
                                ra = Longitude(clusterList[j][i][0][k],unit=u.deg)
                                dec = Latitude(clusterList[j][i][1][k],unit=u.deg)
                                coords = SkyCoord(ra,dec,frame=FK5)
                                region = regions.CircleSkyRegion(coords,radius=1.5*u.arcsec,meta={"label":"#{0},z={1}".format(count,clusterList[j][i][2][k])})
                                regionsList.append(region)
                            count+=1
                    else:
                        for i in range(len(clusterList[j])):
                            for j in range(len(clusterList[j][i][0])):
                                ra = Longitude(clusterList[j][i][0][k],unit=u.deg)
                                dec = Latitude(clusterList[j][i][1][k],unit=u.deg)
                                coords = SkyCoord(ra,dec,frame=FK5)
                                region = regions.CircleSkyRegion(coords,radius=1.5*u.arcsec,meta={"label":"#{0}".format(i+1)})
                                regionsList.append(region)
                            count += 1
                regionsList.write(saveDS9name,overwrite=True,format='ds9')
            if zKey != None:
                return [RA, DEC, Z, sigmaList, SignalToNoiseList], clusterList, ODRlist
            else:
                return [RA, DEC, sigmaList, SignalToNoiseList], ODRlist
        
        elif method.lower() == "dbscan":
            import matplotlib.pyplot as plt
            if zKey != None:
                z_phot = Z
            else:
                z_phot = np.nan*np.ones(N_sources)
            """
            RAnnan = ~np.isnan(RA)
            DECnnan = ~np.isnan(DEC)
            Znnan = ~np.isnan(z_phot)
            ids = RAnnan & DECnnan & Znnan
            X = np.array([RA[ids],DEC[ids]])
            z_phot = z_phot[ids]
            """
            
            X = np.array([RA, DEC])
            # For accurate clustering, normalize data:
            meanRA, stdRA = np.nanmean(X[0]), np.nanstd(X[0])
            meanDEC, stdDEC = np.nanmean(X[1]), np.nanstd(X[1])
            
            print("Starting DBSCAN clustering")
            start = time.time()
            # Normalized X
            X = np.array([(X[0]-meanRA)/stdRA,(X[1]-meanDEC)/stdDEC])
            X = np.transpose(X)
            neigh = neighbors.NearestNeighbors(n_neighbors=sigmas[0],n_jobs=-1).fit(X)
            dists, ind = neigh.kneighbors(X)
            distances = [dists[i][min_c_size-1] for i in range(len(dists))]
            distances = 1/np.asarray(distances)
            
            hist = np.histogram(distances,np.logspace(np.log10(min(distances)),np.log10(max(distances)),30))
            binCenters = np.mean(np.vstack([hist[1][0:-1], hist[1][1:]]), axis=0)
            g_init = models.Gaussian1D(amplitude=np.nanmax(hist[0]), mean=np.nanmean(
                np.log10(distances)), stddev=np.nanstd(np.log10(distances)))
            fit_g = fitting.LevMarLSQFitter()
            g = fit_g(g_init, np.log10(binCenters), hist[0])
            mean,std = g.mean.value, g.stddev.value
            
            try:
                selection = SNlimit[0]*std+mean
            except:
                selection = SNlimit*std+mean
            # Histogram of distances to the third neighbor
            plt.figure()
            plt.plot(np.log10(binCenters),hist[0],'r',label="Histogram")
            plt.plot(np.linspace(min(np.log10(binCenters)),max(np.log10(binCenters)),100),g(np.linspace(min(np.log10(binCenters)),max(np.log10(binCenters)),100)),'k',label="Gaussian fit")
            plt.plot([selection, selection],[min(hist[0]),max(hist[0])],'--',label="SNlimit",color="gray")
            plt.xlabel("$\mathrm{log}_{10}$ 1/(%d-dist)" % (sigmas[0]),fontsize=14)
            plt.ylabel("N",fontsize=14)
            plt.title("Histogram of $\mathrm{log}_{10}$ 1/(%d-dist)" % sigmas[0],fontsize=18)
            plt.legend()
            
            # Choice of overdense regions:
            try:
                ids = (np.log10(distances)-mean)/std >= SNlimit[0]
            except:
                ids = (np.log10(distances)-mean)/std >= SNlimit
            epsVal = 1/min(distances[ids])
            print("Eps at chosen significance limit: {0:.3f}".format(epsVal))
            
            # eps could be an input parameter, or interactive selection, or machine selected
            y = np.zeros((X.shape[0], 1))
            model = cluster.DBSCAN(eps=epsVal, min_samples=min_c_size,n_jobs=-1).fit(X)
            centroids = model.components_
            cls = model.labels_
            cRAs = np.empty(max(cls)+1, dtype=object)
            cDECs = np.empty(max(cls)+1, dtype=object)
            cZs = np.empty(max(cls)+1, dtype=object)
            i = 0
            defLabelList = []
            for label in cls:
                if X[i][0] >= (150.3415-meanRA)/stdRA and X[i][0] <= (150.3489-meanRA)/stdRA and X[i][1] >= (2.3314-meanDEC)/stdDEC and X[i][1] <= (2.3368-meanDEC)/stdDEC:
                    test = label
                if X[i][0] >= (150.2245-meanRA)/stdRA and X[i][0] <= (150.2400-meanRA)/stdRA and X[i][1] >= (2.335-meanDEC)/stdDEC and X[i][1] <= (2.357-meanDEC)/stdDEC:
                    test2 = label
                
                if label != -1:
                    if label not in defLabelList:
                        cRAs[label] = X[i][0]*stdRA+meanRA
                        cDECs[label] = X[i][1]*stdDEC+meanDEC
                        cZs[label] = z_phot[i]
                        defLabelList.append(label)
                    else:
                        cRAs[label] = np.append(cRAs[label],X[i][0]*stdRA+meanRA)
                        cDECs[label] = np.append(cDECs[label],X[i][1]*stdDEC+meanDEC)
                        cZs[label] = np.append(cZs[label],z_phot[i])
                i = i + 1
            print("Number of overdense regions found: {}".format(max(defLabelList)+1))
            cs = np.empty((len(cRAs)), dtype=object)
            tcs = np.empty(0, dtype=int)
            for i in range(len(cs)):
                if zKey != None:
                    cs[i] = np.array([cRAs[i], cDECs[i], cZs[i]])
                    try:
                        if len(cRAs[i]) >= min_c_size:
                            tcs = np.append(tcs, i)
                    except:
                        pass
                else:
                    cs[i] = np.array([cRAs[i], cDECs[i]])
                    try:
                        if len(cRAs[i]) >= min_c_size:
                            tcs = np.append(tcs, i)
                    except:
                        pass
            ODRcs = cs
            
            cs = cs[tcs]
            clusterCands = []

            if zKey != None:
                for i in range(len(cs)):
                    temp = cs[i]
                    sZS = np.argsort(temp[2])
                    temp = np.array([temp[0][sZS], temp[1][sZS], temp[2][sZS]])
                    clusterz = abs((temp[2][:-1]-temp[2][1:])
                                   ) >= 0.10*(1+temp[2][:-1])
                    N_clustersz = sum(clusterz)+1
                    if N_clustersz != 1:
                        climz = np.where(clusterz == True)
                        tRA, tDEC, tZ = np.split(
                            temp[0], climz[0]+1), np.split(temp[1], climz[0]+1), np.split(temp[2], climz[0]+1)
                        for i in range(len(tRA)):
                            if len(tRA[i]) >= min_c_size:
                                clusterCands.append(
                                    np.array([tRA[i], tDEC[i], tZ[i]]))
                    else:
                        tRA, tDEC, tZ = temp[0], temp[1], temp[2]
                        if len(tRA) >= min_c_size:
                            clusterCands.append(np.array([tRA, tDEC, tZ]))
            else:
                clusterCands = cs
            stop = time.time()
            print("Time elapsed: {0:.2f}s".format(stop-start))
            Lcs, Lcenters = cs[:], cls[:]
            #Lcenters[0],Lcenters[1] = Lcenters[0]*stdRA+meanRA,Lcenters[1]*stdDEC+meanDEC
            # plot data points color-coded by class, cluster markers and centroids
            # hold(True)
        
            
            import random
            from scipy import interpolate
            
            """
            meshRA, meshDEC = np.meshgrid(np.linspace(min(TcRAs),max(TcRAs),round(np.sqrt(len(TcRAs)))),np.linspace(min(TcDECs),max(TcDECs),round(np.sqrt(len(TcRAs)))))
            interp = interpolate.NearestNDInterpolator(list(zip(TcRAs,TcDECs)),totalVariance[i])#,fill_value=-2)#,kind='cubic')
            intSN = interp(meshRA,meshDEC)
            """
            def get_colors(n): return list(
                map(lambda i: "#" + "%06x" % random.randint(0x333333, 0x555555), range(n)))
            colors = get_colors(len(clusterCands))
            plt.figure()
            
            plt.scatter(RA,DEC,1,'orange',alpha=0.5)
            for i in range(len(clusterCands)):
                plt.scatter(clusterCands[i][0], clusterCands[i][1], 40, facecolors="none", edgecolors=colors[i])
                if zKey != None:
                    plt.annotate("#{0:d},z~={1:.2f}".format(i+1, np.median(clusterCands[i][2])), [
                                 np.mean(clusterCands[i][0])+0.1*(max(RA)-min(RA)), min(clusterCands[i][1])-0.05*(max(DEC)-min(DEC))], color=colors[i])
                else:
                    plt.annotate("#{0:d}".format(i+1), [
                                 np.mean(clusterCands[i][0])+0.1*(max(RA)-min(RA)), np.mean(clusterCands[i][1])-0.05*(max(DEC)-min(DEC))], color=colors[i])
            plt.xlim([max(RA), min(RA)])
            plt.ylim([min(DEC), max(DEC)])
            plt.ylabel("DEC [deg]",fontsize=14)
            plt.xlabel("RA [deg]",fontsize=14)
            plt.title("Cluster candidate map, S/N $\geq$ {0:.2f}".format(SNlimit[0]),fontsize=18)
            
            if zKey != None:
                z_cut = min(Z)
                zPhotHistPlot = input(
                    "Do you wish to plot histograms of z_phot for all overdensity regions with more than {} sources? (y/n): ".format(min_c_size))
                if zPhotHistPlot.lower() == "y" or zPhotHistPlot.lower() == "yes":
                    # z_phot histograms for all significant overdensity regions
                    for i in range(len(ODRcs)):
                        plt.figure()
                        hist = np.histogram(ODRcs[i][2], bins=np.linspace(
                            z_cut, 6, round((6-z_cut)/0.2)))
                        binCenters = np.mean(
                            np.vstack([hist[1][0:-1], hist[1][1:]]), axis=0)
                        # plt.hist(ODRlist[i][j][2],bins=20)
                        plt.plot(binCenters, hist[0], 'o-')
                        plt.xlabel("z_phot")
                        plt.title("Histogram of redshift, RA: {0:.4f} DEC: {1:.4f}".format(
                            np.mean(ODRcs[i][0]), np.mean(ODRcs[i][1])))
            saveODR = input(
                "Do you wish to save the candidates in a .fits file? (y/n): ")
            if saveODR.lower() == "y" or saveODR.lower() == "yes":
                saveODRName = input(
                    "Input name of file you would like to save candidates to: ")
                from astropy.io import fits
                import astropy.table as tb
                t = tb.Table()
                prevN = 1
                if zKey != None:
                    for i in range(len(clusterCands)):
                        tabletemp = tb.Table(data=[np.ones(len(clusterCands[i][0]))*prevN, clusterCands[i][0].astype(float), clusterCands[i][1].astype(
                            float), clusterCands[i][2].astype(float)], names=("Cand_num", "RA", "DEC", "z"))
                        t = tb.vstack([t, tabletemp])
                        prevN += 1
                else:
                    for i in range(len(clusterCands)):
                        tabletemp = tb.Table(data=[np.ones(len(clusterCands[i][0]))*prevN, clusterCands[i][0].astype(float), clusterCands[i][1].astype(
                            float)], names=("Cand_num", "RA", "DEC"))
                        t = tb.vstack([t, tabletemp])
                        prevN += 1
                t.write(saveODRName, format='fits')
            saveDS9 = input(
                "Do you wish to save the candidates in a DS9 region .reg file? (y/n): ")
            if saveDS9.lower() == "y" or saveDS9.lower() == "yes":
                from astropy.coordinates import SkyCoord,FK5,Longitude,Latitude
                import astropy.units as u
                import regions
                saveDS9name = input("Please input the filename for the DS9 region file: ")
                regionsList = regions.Regions([])
                if zKey != None:
                    for i in range(len(clusterCands)):
                        for j in range(len(clusterCands[i][0])):
                            ra = Longitude(clusterCands[i][0][j],unit=u.deg)
                            dec = Latitude(clusterCands[i][1][j],unit=u.deg)
                            coords = SkyCoord(ra,dec,frame=FK5)
                            region = regions.CircleSkyRegion(coords,radius=1.5*u.arcsec,meta={"label":"#{0},z={1}".format(i+1,clusterCands[i][2][j])})
                            regionsList.append(region)
                else:
                    for i in range(len(clusterCands)):
                        for j in range(len(clusterCands[i][0])):
                            ra = Longitude(clusterCands[i][0][j],unit=u.deg)
                            dec = Latitude(clusterCands[i][1][j],unit=u.deg)
                            coords = SkyCoord(ra,dec,frame=FK5)
                            region = regions.CircleSkyRegion(coords,radius=1.5*u.arcsec,meta={"label":"#{0}".format(i+1)})
                            regionsList.append(region)
                regionsList.write(saveDS9name,overwrite=True,format='ds9')
            return Lcs, Lcenters, clusterCands
