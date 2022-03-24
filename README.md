# Synthesis-Project
Synthesis project of Nikolaj Bjerregaard Sillassen, Spring 2022 at DTU

A program to identify cluster candidates given a pre-selected astrophysical dataset in .fits format.
Can be run either as an interactive program, or as a function by inputting the necessary variables.

The program can run with or without redshifts, but will only produce cluster candidates with redshift. Without redshift, only overdensity regions are produced.

Two methods of identification, SigmaN and DBSCAN.
SigmaN uses the distance between each source and the N'th nearest neighbor, following:
SigmaN = N/(Pi d_(Nth neighbor)^2).
SigmaN produces a density plot of the entire region. Marking overdensity regions above a given signal to noise, and optionally marking cluster candidates.

Whereas DBSCAN uses Density Based Spacial Clustering for Applications with Noise, to identify cluster candidates. 
DBSCAN only produces a plot of cluster candidates. Without redshift DBSCAN produces overdensity regions.

In the end, you are given the possibility to save outputs in a fits file.

Automated Identification of Galaxy Clusters
galcluster v0.1
Program created during Synthesis Project at DTU Space, in spring 2022

```
## Function description:
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
