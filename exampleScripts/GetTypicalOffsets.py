# GetTypicalOffsets.py
# 
# A code to calculate the position offset of each extracted source
# relative to its average runcat position.
# A histogram is produced using all of the offsets and then fit with
# a Gaussian distribution.

# Import all the dependencies and generic setup
import scipy as sp
import numpy as np
import pandas as pd
import sqlalchemy
from sqlalchemy import *
from sqlalchemy.orm import relationship
import tkp.db
import logging
logging.basicConfig(level=logging.INFO)
query_loglevel = logging.WARNING  # Set to INFO to see queries, otherwise WARNING
import sys
sys.path.append('../')
from dblogin import * # This file contains all the variables required to connect to the database
from databaseTools import dbtools
from tools import tools
from plotting import plot_varib_params as pltvp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties
from tkp.db.model import Runningcatalog
from tkp.db.model import Extractedsource
from tkp.db.model import Assocxtrsource
import astropy.units as u
from astropy.coordinates import SkyCoord
import os
from scipy.optimize import curve_fit

def gaussian_func(x, a, x0, sigma):
    if x0 < 0: # force the mean to be >0
        return 10000.
    else:
        return a * np.exp(-(x-x0)**2/(2*sigma**2))
    
# The input database, dataset and thresholds
dataset_id = 28
database = 'AR_testing6'

outfile='ds'+str(dataset_id)+'_offsets.csv'

if os.path.exists(outfile):
    separations = np.genfromtxt(outfile, delimiter=',')
else:
    # Connect to the database and run the queries
    session = dbtools.access(engine,host,port,user,password,database)

    runcats= session.query(Runningcatalog).filter(Runningcatalog.dataset_id == dataset_id).order_by(Runningcatalog.id).all()

    separations = [] # in arcsec
    for runcat in runcats:
        print(runcat.id)
        flxVals = session.query(Assocxtrsource,Extractedsource).select_from(join(Assocxtrsource,Extractedsource)).filter(Assocxtrsource.runcat_id == runcat.id).all()
        c1 = SkyCoord(runcat.wm_ra, runcat.wm_decl, unit="deg")
        for src in flxVals:
            if src.Extractedsource.extract_type==0: # only blindly fit sources
                if src.Assocxtrsource.type==3: # only 1-to-1 associations
                    c2 = SkyCoord(src.Extractedsource.ra, src.Extractedsource.decl, unit="deg")
                    sep = c1.separation(c2)
                    separations.append(sep.arcsecond)
    np.savetxt(outfile, separations, delimiter=',')

nbins=100
array = plt.hist(separations,bins=nbins,histtype='stepfilled',density=False)#,log=True)

ydata = array[0]
x=array[1]
xdata = [(((x[n]-x[n-1])/2.)+x[n]) for n in range(len(x)-1)]

initial_guess = [max(ydata),np.median(separations),np.std(separations)]
popt, pcov = curve_fit(gaussian_func, xdata, ydata,p0=initial_guess)

xplot = np.linspace(0,60,1000)
plt.plot(xplot,gaussian_func(xplot,*popt))
print('mean = '+str(popt[1])+', sigma = '+str(popt[2]))


plt.xlim(0,30)
plt.xlabel('Offset (arcsec)')
plt.ylabel('Number of sources')
plt.savefig('ds'+str(dataset_id)+'_offsets.png')

exit()

