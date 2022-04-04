# Synthesis-Project
Synthesis project of Nikolaj Bjerregaard Sillassen, Spring 2022 at DTU

A program to identify cluster candidates given an astrophysical dataset in .fits format.
Can be run either as an interactive program, or as a function by inputting the necessary variables.

The program can run with or without preselection in flux or mag, it is highly recommended to use preselection to avoid large data.
It can also run with or without minimum preselection in redshifts.

Two methods of identification, SigmaN and K_means.
SigmaN uses the distance between each source and the N'th nearest neighbor, following:
SigmaN = N/(Pi d_(Nth neighbor)^2)
Whereas K_means uses K_means hierarchical clustering, to try to identify clusters.

In the end, you are given the possibility to save outputs in a fits file.

Automated Identification of Galaxy Clusters
galcluster v0.1
Program created during Synthesis Project at DTU Space, in spring 2022
    
## Function description:
```     
galcluster(filename=None,FluxPreselection=None,method="sigma",sigmas=[5,10],multi_processing=True,mute_plots=False,SNlimit=5,min_c_size=3,z_cut=1.5)

### Input Parameters

  filename : str, optional
      Name of the file to run clustering algorithm on, must be .fits format 

  method : str, optional 
      "sigma" (default) will use sigma5 and sigma10 overdensities,
      "kmeans" will use K-means spacial clustering,
      "deep" will use convnet to cluster an input image

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

  SN_limit : float, optional
      limit the signal to noise of overdensity for clustering, 5 is default

  min_c_size : int, optional 
      minumum number of sources in a cluster candidate, 3 is default
  
  zPhotKey : str, optional
      keyword for column containing photometric redshift
  
  zSpecKey : str, optional
      keyword for column containing spectroscopic redshift
      
  z_cut : positive float, optional 
      minimum redshift of the sources that are being clustered, 1.5 is 
      default. 

 ### Returns

  out : [RA,DEC,Z,Sigmas,SignalToNoiseSigmas], clusterList, ODRlist

  First list:
      RA : ndarray 
          array containing right ascension for all pre-selected sources
      DEC : ndarray
          array containing declination for all pre-selected sources
      Z : ndarray
          array containing photometric redshifts for all pre-selected 
          sources
      Sigmas : ndarray
          array containing all sigma values, for all pre-selected sources
      SignalToNoiseSigmas : ndarray
          array containing all signal to noise values of all sigmas, for 
          all pre-selected sources
  clusterList : list
      First Dimension : list
          List of cluster candidate lists, for all values of input sigmas
      Second Dimension : list
          [RA,DEC,z_phot,Sigma,SignalToNoiseSigma] for all sources in the
          candidates
  ODRlist : list
      First Dimension : list
          List of overdensity region lists, for all input sigmas
      Second Dimension : list
          [RA,DEC,z_phot,Sigma,SignalToNoiseSigma] for all sources in the
          overdensity regions
```

  @author: Nikolaj Bjerregaard Sillassen
