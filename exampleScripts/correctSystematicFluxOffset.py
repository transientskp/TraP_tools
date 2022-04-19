#
# A code to correct any systemmatic flux density offsets in radio images.
# First run TraP on the images with a high detection threshold to get bright sources.
# Then run this script to obtain the multiplicative flux density correction for each individual image.
# It compares the extracted flux density to the average flux density for each source.
# by fitting a straight line between the datapoints, the gradient gives the multiplicative flux density factor.
# The key assumption is that the majority of sources are stable.
#
# Advice for use
# - select bright, point-like sources within the FWHM of the beam.
# - Use the TraP monitoring list capability (with a very high detection threshold) to only get the sources chosen

# note - only process one frequency at once

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
import pylab
pylab.rcParams['legend.loc'] = 'best'
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties
from tkp.db.model import Runningcatalog
from tkp.db.model import RunningcatalogFlux
from tkp.db.model import Image
from tkp.db.model import Extractedsource
from tkp.db.model import Assocxtrsource
from astropy.io import fits


# The input database, dataset and thresholds
dataset_id = 20
database = 'AR_testing4'
plots=True
saveFigUrl='/scratch/antoniar/JelleVariable/1820/AR_images/tmp/'

# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)
VarParams = dbtools.GetVarParams(session,dataset_id)


# Obtain the average flux density measurement for each source
# RuncatID, Avg Int Flux Density
Runcats = session.query(Runningcatalog,RunningcatalogFlux).select_from(join(Runningcatalog,RunningcatalogFlux)).filter(Runningcatalog.dataset_id == dataset_id).all()
RuncatsData=[[Runcats[i].Runningcatalog.id,Runcats[i].RunningcatalogFlux.avg_f_int,Runcats[i].Runningcatalog.datapoints] for i in range(len(Runcats))]
RuncatsData = pd.DataFrame(data=RuncatsData, columns=['RuncatID','avg_f_int','datapoints'])


# Get all the extracted sources with their runcat ids
listRuncats = RuncatsData.RuncatID
listRuncats = list(listRuncats)
listRuncats.sort()
AssocSrcs = session.query(Assocxtrsource).filter(Assocxtrsource.runcat_id.in_(listRuncats)).all()
AssocSrcs = [[AssocSrcs[i].xtrsrc.id, AssocSrcs[i].runcat.id] for i in range(len(AssocSrcs))]
AssocSrcs = pd.DataFrame(data=AssocSrcs, columns=['xtrsrc','RuncatID'])

# Loop through the images
Images = session.query(Image).filter(Image.dataset_id == dataset_id).all()
numImgs=len(Images)
RuncatsData = RuncatsData.drop(RuncatsData[RuncatsData.datapoints < numImgs].index)

for a in range(len(Images)):
    # For each image obtain the extracted source flux density for each RuncatID
    extractedSrcs = session.query(Extractedsource).filter(Extractedsource.image_id == Images[a].id).all()
    extractedSrcsData = [[extractedSrcs[b].id,extractedSrcs[b].f_int,extractedSrcs[b].f_int_err] for b in range(len(extractedSrcs))]
    extractedSrcsData = pd.DataFrame(data=extractedSrcsData, columns=['xtrsrc','f_int','f_int_err'])
    extractedSrcsData = pd.merge(extractedSrcsData,AssocSrcs, on='xtrsrc')
    plotData = pd.merge(extractedSrcsData,RuncatsData,on='RuncatID')

    # Fit straight line through the origin and output the gradient
    x=plotData.avg_f_int
    y=plotData.f_int
    x = x[:,np.newaxis]
    b, _, _, _ = np.linalg.lstsq(x,y) # fit a straight line going through the origin   
    print Images[a].url, b[0]
    xrng = np.linspace(0,max(plotData.avg_f_int), 10)
    xfit = [x*b for x in xrng]
    if plots==True:
        # Plot scatter graph of data points
        plt.errorbar(plotData.avg_f_int, plotData.f_int, yerr=plotData.f_int_err, fmt='o')
        plt.plot(xrng,xfit,'-')
        plt.show()

    #Open each image and multiply the pixel values by the gradient - save each image to new location
    hdu = fits.open(Images[a].url)
    data = hdu[0].data[0][0]
    hdu[0].data[0][0] = data / b[0]
    newFilename = Images[a].url.split('/')[-1]
    newFilename = saveFigUrl+'corrected_'+newFilename
    hdu.writeto(newFilename)
 
