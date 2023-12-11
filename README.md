[GalCluster-logo](GalCluster.png)
# GalCluster
Synthesis project of Nikolaj Bjerregaard Sillassen, Spring 2022 at DTU

A program to identify cluster candidates given an astrophysical dataset in .fits format.
Can be run either as an interactive program, or as a function by inputting the necessary variables.

Two methods of identification, SigmaN and DBSCAN.
SigmaN uses the distance between each source and the N'th nearest neighbor, following:
SigmaN = N/(Pi d_(Nth neighbor)^2)
Whereas DBSCAN uses DBSCAN clustering [SKLEARN Documentation](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html).

In the end, you are given the possibility to save outputs in a fits file, and all identified clusters as a DS9 region file.

Automated Identification of Galaxy Clusters
GalCluster v0.1
Program created during Synthesis Project at DTU Space, in spring 2022
    
## Function description:
```     
GalCluster(filename=None, RAKey="RA", DECKey="DEC", method="sigma", sigmas=[5,10], multi_processing=True, mute_plots=False, SNlimit=[5,5],min_c_size=3,zKey=None)

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
```

  @author: Nikolaj Bjerregaard Sillassen
