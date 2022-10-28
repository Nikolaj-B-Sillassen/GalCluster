# -*- coding: utf-8 -*-
"""
Example of how to run galcluster

@author: Nikolaj Bjerregaard Sillassen
"""
import matplotlib.pyplot as plt
import GalCluster as gc

plt.close('all') # Close all plt figures
data_path = "PATH TO FITS FILE" # Define path of data
dataset = "mydata.fits" # Define dataset name
if __name__ == "__main__":
    # Running file as a function
    output = gc.GalCluster(data_path+dataset,RAKey="RA",DECKey="DEC",SNlimit=[3,3],min_c_size=3,zKey="z500")
    # Running file as an interactive program
    #output = gc.galcluster()