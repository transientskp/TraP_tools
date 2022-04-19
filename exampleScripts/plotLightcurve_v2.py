# plotLightcurve.py
#
# a code to plot a multi-frequency lightcurve for a source

# Import all the dependencies and generic setup
import scipy as sp
import numpy as np
import pandas as pd
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
from astropy import units as u
from astropy.coordinates import SkyCoord
from datetime import datetime

from tkp.db.model import Extractedsource
from tkp.db.model import Assocxtrsource
from tkp.db.model import Image

# The input database, dataset and thresholds
dataset = 25
sourceID = 59613

freqs = pd.DataFrame(data=[[0.963,'#e41a1c','o'],
                                [1.177,'#a65628','s'],
                                [1.391,'#377eb8','v'],
                                [1.605,'#984ea3','D']],
                          columns = ['freq','colour','marker'])
database = 'AR_testing4'
outname = str(sourceID)+'_lightcurve.png' # name of the plotdata datafile

# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)


# setting up the plot
nullfmt   = NullFormatter()         # no labels
fontP = FontProperties()
fontP.set_size('xx-large')
plt.figure(1,figsize=(12,6))
plt.tight_layout()

flxVals = session.query(Assocxtrsource,Extractedsource).select_from(join(Assocxtrsource,Extractedsource)).filter(Assocxtrsource.runcat_id == sourceID).all()

lightcurve = pd.DataFrame(data=[[flxVals[x].Extractedsource.image.id, flxVals[x].Extractedsource.f_int, flxVals[x].Extractedsource.f_int_err, flxVals[x].Extractedsource.det_sigma] for x in range(len(flxVals))], columns = ['Image','Flux','FluxErr','SNR'])

images= session.query(Image).all()
images = pd.DataFrame(data=[[images[x].id,images[x].taustart_ts,np.around(images[x].band.freq_central/1e9, decimals=3)] for x in range(len(images))], columns=['Image','Time','Freq'])

lightcurve = pd.merge(lightcurve, images, on="Image")


for index, row in freqs.iterrows():
    print row.freq
    lightcurveTMP = lightcurve.loc[lightcurve['Freq'] == row.freq]
    plt.errorbar(lightcurveTMP.Time,lightcurveTMP.Flux*1e3, yerr=lightcurveTMP.FluxErr*1e3, c=row.colour, fmt=row.marker, markersize=7, linestyle='-')

#xmin, xmax, ymin, ymax = plt.axis()
#ymin_new = 10.**(np.log10(ymin)*.9)
#plt.ylim(ymin_new,ymax)
    
plt.legend(['0.96 GHz','1.18 GHz','1.39 GHz','1.61 GHz'], loc=4)

plt.xlabel('Date', fontsize=16)
plt.ylabel('Flux Density (mJy)', fontsize=16)
plt.savefig(outname)

exit()
    
    
