# -*- coding: utf-8 -*-
"""
Example of how to run galcluster

@author: Nikolaj Bjerregaard Sillassen
"""
import matplotlib.pyplot as plt
import GalCluster as gc

plt.close('all') # Close all plt figures
data_path = "../../../Data/" # Define path of data
dataset = "COSMOS_Redshift2Selected.fits" # Define dataset name
if __name__ == "__main__":
    # Running file as a function
    output = gc.GalCluster(data_path+dataset,RAKey="ALPHA_J2000_1",DECKey="DELTA_J2000_1",method="dbscan",sigmas=[5],SNlimit=[3.5],min_c_size=3,zKey="lp_zPDF_1")
    # Running file as an interactive program
    #output = gc.galcluster()