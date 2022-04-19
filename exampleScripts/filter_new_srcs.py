# filter_new_srcs.py
# 
# A code to conduct simple filtering steps on the new source list output by TraP.
# Sources are removed that are associated with other sources in the initial image
# and then filtered on their reduced weighted chi^2 to remove sources close to the
# background noise.
#

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
import pylab
pylab.rcParams['legend.loc'] = 'best'
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.stats import norm

from tkp.db.model import Extractedsource
from sqlalchemy import *
from sqlalchemy.orm import relationship
from tkp.db.model import Image
from tkp.db.model import Assocxtrsource
from tkp.db.model import Varmetric
from tkp.db.model import Skyregion
from tkp.db.model import Frequencyband


def gaussian_fit(data,param,xmin,xmax,nbins):
    range_data=np.linspace(xmin,xmax,nbins)
    fit=norm.pdf(range_data,loc=param[0],scale=param[1])
    return range_data,fit


# The input database, dataset and thresholds
dataset_id = 38
database = 'AR_testing4'
websiteURL = 'http://banana.transientskp.org/r3/vlo_'+database+'/runningcatalog/'
BMaj = 0.002185 # in degrees
beamwidths=5.
SrcAssocRadius = BMaj * beamwidths # in degrees, 10" is ~2.8e-3
#print SrcAssocRadius
sigma1 = 1.5 # Threshold on the reduced weighted chi^2
sigma2 = 1.5 # Threshold on the variability parameter
sigmaExcessLimit = 5
FilterRadius = 0.051 # in degrees 0.017~1'
maxSigma = 5.56 # Threshold for the transient search
minSigma = 4 # Threshold for associated sources in different frequencies
DeepCatalog = pd.DataFrame(data=[[0.963,'/scratch/antoniar/JelleVariable/1820/AR_images/freq0_combined.csv'],[1.177,'/scratch/antoniar/JelleVariable/1820/AR_images/freq1_combined.csv'],[1.391,'/scratch/antoniar/JelleVariable/1820/AR_images/freq2_combined.csv'],[1.605,'/scratch/antoniar/JelleVariable/1820/AR_images/freq3_combined.csv']], columns = ['freq','path']) # links to the extracted sources in a deep image at each observing frequency considered


# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)
NewSrcs = dbtools.GetNewSrcs(session,dataset_id)       # Get the new source table
VarParams = dbtools.GetVarParams(session,dataset_id)   # Get the running catalogue and varmetric catalogues and combine
Img1Srcs = dbtools.GetImg1Srcs(session,dataset_id)     # Get all the sources identified in the first image of the dataset

freqBands = session.query(Frequencyband).all()
freqBands = [[freqBands[x].id,np.around(freqBands[x].freq_central/1e9, decimals=3)] for x in range(len(freqBands))]
freqBands = pd.DataFrame(data=freqBands, columns=['id','freq'])

# Get co-ordinates, fluxes and other information for the new sources
NewSrcData=[[NewSrcs[i].Runningcatalog.id,NewSrcs[i].Runningcatalog.wm_ra,NewSrcs[i].Runningcatalog.wm_decl,NewSrcs[i].Newsource.trigger_xtrsrc.id,NewSrcs[i].Newsource.trigger_xtrsrc.image.taustart_ts,NewSrcs[i].Newsource.trigger_xtrsrc.image.band.id,NewSrcs[i].Newsource.trigger_xtrsrc.det_sigma] for i in range(len(NewSrcs))]
NewSrcDataFrame = pd.DataFrame(data=NewSrcData, columns=['RuncatID','ra','decl','xtrsrc','TimeDetect','Band','detSigma'])
NewSrcDataFrame = NewSrcDataFrame.sort_values(by=['RuncatID'])

print "Number of new sources: "+str(len(NewSrcDataFrame))


# Filter out all new sources that were below the detection threshold for transient sources
NewSrcDataFrame = NewSrcDataFrame.loc[NewSrcDataFrame['detSigma'] >= maxSigma]

# Make a list of all the positions of the new sources 
NewSrcs_coord = SkyCoord(ra=(NewSrcDataFrame.ra*u.degree).values,dec=(NewSrcDataFrame.decl*u.degree).values)

VarData = [[VarParams[i].Runningcatalog.id, VarParams[i].Varmetric.eta_int, VarParams[i].Varmetric.v_int, VarParams[i].Varmetric.lightcurve_max, VarParams[i].Varmetric.lightcurve_median, (VarParams[i].Varmetric.band.freq_central/1e6), VarParams[i].Runningcatalog.datapoints, VarParams[i].Varmetric.newsource] for i in range(len(VarParams))]
VarData = pd.DataFrame(data=VarData,columns=['RuncatID','eta','V','maxFlx','avgFlx','freq','dpts','newSrc'])
VarData = VarData.sort_values(by=['RuncatID'])
VarData = VarData.fillna('N')
NewSrcDataFrame = pd.merge(NewSrcDataFrame, VarData, on="RuncatID")

print "Number of new sources: "+str(len(NewSrcs_coord))


# Get co-ordinates of all sources in the first image
Img1VarParamsData=[[VarParams[i].Runningcatalog.id,VarParams[i].Runningcatalog.wm_ra,VarParams[i].Runningcatalog.wm_decl] for i in range(len(VarParams)) if VarParams[i].Runningcatalog.id in Img1Srcs]
Img1VarParamsDataFrame = pd.DataFrame(data=Img1VarParamsData, columns=['RuncatID','ra','decl'])
Img1VarParamsDataFrame = Img1VarParamsDataFrame.sort_values(by=['RuncatID'])
Img1VarParams_coord = SkyCoord(ra=(Img1VarParamsDataFrame.ra*u.degree).values,dec=(Img1VarParamsDataFrame.decl*u.degree).values)

# From CatalogueMatching.ipynb
# idx, d2d, d3d = aart_coord.match_to_catalog_sky(green_coord)
#idx is an array of indices for the green Dataframe, such that green.iloc[idx[i]] is the nearest source to aart.iloc[i].  len(idx) = len(aart)
#d2d is the pointwise 2D angular distances. len(d2d) = len(aart)
#d3d is the 3D distance, usefull if distance to each source is known, otherwise they're assumed to be on a unit sphere. To add physical distances to the catalogued soueces see http://docs.astropy.org/en/stable/coordinates/matchsep.html. len(d3d) = len(aart)


###### Filter 1 - remove any poorly associated sources
###### i.e. those that were detected in the first image but TraP incorrectly labeled as new sources
idx, d2d, d3d = NewSrcs_coord.match_to_catalog_sky(Img1VarParams_coord)
UnassocNewSrcs = NewSrcDataFrame[d2d.deg > SrcAssocRadius] # Unassociated sources - those not found in Img1
UnassocNewSrcs_coord = NewSrcs_coord[d2d.deg > SrcAssocRadius] # Unassociated sources - those not found in Img1

print "Number of new sources after Filter 1: "+str(len(UnassocNewSrcs_coord))

uniqueIDlist=list(UnassocNewSrcs.RuncatID)


###### Filter 2 - Reject sources too close to the source extraction radius.
# Not general enough at the moment...

skyreg = session.query(Skyregion).select_from(Skyregion).filter(Skyregion.dataset_id == dataset_id).one()
centre = SkyCoord(ra=(skyreg.centre_ra*u.degree),dec=(skyreg.centre_decl*u.degree))
xtrRadius = skyreg.xtr_radius

uniqueIDlistTMP=[]
for a in range(len(UnassocNewSrcs)):
    sep = UnassocNewSrcs_coord[a].separation(centre)
    if sep.degree < xtrRadius - FilterRadius:
        uniqueIDlistTMP.append(UnassocNewSrcs.iloc[a].RuncatID)

uniqueIDlist=uniqueIDlistTMP

print "Number of new sources after Filter 2: "+str(len(uniqueIDlist))


###### Filter 3 - Needs to be detected in at least 2 frequencies at the same time
# get all images from this time - i.e. matching taustart_ts
# get list of extracted sources from those images and associated runcat ids
# find number of associated extracted sources - i.e. where runcat matches that of new source
# 


uniqueIDlist_multi = []


for src in range(len(UnassocNewSrcs)):
    time= UnassocNewSrcs.iloc[src].TimeDetect
    runcat = UnassocNewSrcs.iloc[src].RuncatID
    band = UnassocNewSrcs.iloc[src].Band
    if runcat in uniqueIDlist:
        numDetect1stTime=1
    
        images= session.query(Image).select_from(Image).filter(Image.taustart_ts == time).filter(Image.dataset_id == dataset_id).all()
        for img in images:
            if band != img.band.id:
               #imgSrcs= session.query(Assocxtrsource,Extractedsource).select_from(join(Assocxtrsource,Extractedsource)).filter(Assocxtrsource.runcat_id == runcat).filter(Extractedsource.image_id == img.id).all()
               imgSrcs= session.query(Assocxtrsource,Extractedsource).select_from(join(Assocxtrsource,Extractedsource)).filter(Assocxtrsource.runcat_id == runcat).filter(Extractedsource.image_id == img.id).filter(Extractedsource.det_sigma >= minSigma).all()
               numDetect1stTime = numDetect1stTime + len(imgSrcs)
        
        if numDetect1stTime > 1:
            uniqueIDlist_multi.append(runcat)
print uniqueIDlist_multi
 
#uniqueIDlist_multi=uniqueIDlist

print "Number of new sources after Filter 3: "+str(len(uniqueIDlist_multi))


###### Filter 4 - is either not detected in deep image or is significantly brigter


# conduct a source association between the transient sources and the catalogues
# of the deep images at each observing frequency. If not associated then it is a
# transient candidate and if significantly brighter it is a variable candidate


IDlist_deep = []
bandTMP = []
uniqueIDlistTMP=[]
srclist=[]
counter=0
#for src in uniqueIDlist_multi:
for src in range(len(UnassocNewSrcs)):
    if UnassocNewSrcs.iloc[src].RuncatID in uniqueIDlist_multi:
        tmpDF=0
        runcat = UnassocNewSrcs.iloc[src].RuncatID
        time= UnassocNewSrcs.iloc[src].TimeDetect
        band = UnassocNewSrcs.iloc[src].Band
        maxFlx = UnassocNewSrcs.iloc[src].maxFlx
        freq = freqBands.loc[freqBands.id==band]['freq']
        combImg = DeepCatalog.loc[np.around(DeepCatalog.freq,3)==freq.values[0]]
        DeepSrcs = pd.read_csv(combImg.path.values[0], header=0,skip_blank_lines=True, sep=', ')
    
    
        tmp_coord1=SkyCoord(ra=(UnassocNewSrcs.iloc[src].ra*u.degree),dec=(UnassocNewSrcs.iloc[src].decl*u.degree))

        separations = [[row[0],row[1],tmp_coord1.separation(SkyCoord(ra=(row[0]*u.degree),dec=(row[1]*u.degree))).degree, row[2], row[3]] for row in DeepSrcs[['ra','dec','int_flux','int_flux_err']].to_numpy()]
        associations = [[row[2],row[3],row[4]] for row in separations if row[2] < SrcAssocRadius]
    
        if not associations:
            srclist.append(runcat)
            print 'not associated: '+str(runcat)
            #tmpDF = UnassocNewSrcs.iloc[src]
            #tmpDF['DeepFlx']=0
            #tmpDF['DeepFlxErr']=0
            #tmpDF['sigmaExcess']=0
            #uniqueIDlistTMP.append(tmpDF)
        else:
            sigmaExcess = (maxFlx - associations[0][1])/associations[0][2]
            if sigmaExcess > sigmaExcessLimit:
                srclist.append(runcat)
                #tmpDF = UnassocNewSrcs.iloc[src]
                #tmpDF['DeepFlx']=associations[0][1]
                #tmpDF['DeepFlxErr']=associations[0][2]
                #tmpDF['sigmaExcess']=sigmaExcess
                #uniqueIDlistTMP.append(tmpDF)
        counter=counter+1
uniqueIDlist_multi=srclist

print "Number of new sources after Filter 4: "+str(len(uniqueIDlist_multi))

#IdTrans=uniqueIDlist_multi.RuncatID.unique()
if len(uniqueIDlist_multi)>0:
    for a in uniqueIDlist_multi:
        print websiteURL+str(a)
else:
    print "No Transients"

exit()

###### Filter 5 - Check variability parameters of remaining new sources

      
VarParams = dbtools.GetVarParams(session,dataset_id)   # Get the running catalogue and varmetric catalogues and combine

VarData = [[VarParams[i].Runningcatalog.id, VarParams[i].Varmetric.eta_int, VarParams[i].Varmetric.v_int, VarParams[i].Varmetric.lightcurve_max, VarParams[i].Varmetric.lightcurve_median, (VarParams[i].Varmetric.band.freq_central/1e6), VarParams[i].Runningcatalog.datapoints, VarParams[i].Varmetric.newsource] for i in range(len(VarParams))  if VarParams[i].Runningcatalog.id in uniqueIDlist_multi]
VarData = pd.DataFrame(data=VarData,columns=['RuncatID','eta','V','maxFlx','avgFlx','freq','dpts','newSrc'])
VarData = VarData.sort_values(by=['RuncatID'])
VarData = VarData.fillna('N')
VarData = VarData[(VarData['eta']>0) & (VarData['V'] > 0)] 
plotdata=VarData

paramx, paramx2 = tools.SigmaFit(np.log10(plotdata['eta']))
paramy, paramy2 = tools.SigmaFit(np.log10(plotdata['V']))

if sigma1 == 0:
    sigcutx = 0
else:
    sigcutx = paramx[1]*sigma1+paramx[0]

if sigma2 == 0:
    sigcuty = 0
else:
    sigcuty = paramy[1]*sigma2+paramy[0]

print('Gaussian Fit eta: '+str(round(10.**paramx[0],2))+'(+'+str(round((10.**(paramx[0]+paramx[1])-10.**paramx[0]),2))+' '+str(round((10.**(paramx[0]-paramx[1])-10.**paramx[0]),2))+')')
print('Gaussian Fit V: '+str(round(10.**paramy[0],2))+'(+'+str(round((10.**(paramy[0]+paramy[1])-10.**paramy[0]),2))+' '+str(round((10.**(paramy[0]-paramy[1])-10.**paramy[0]),2))+')')
print 'Eta_nu threshold='+str(10.**sigcutx)+', V_nu threshold='+str(10.**sigcuty)

frequencies = plotdata.freq.unique()
col = pltvp.make_cmap(frequencies)


nullfmt   = NullFormatter()         # no labels
fontP = FontProperties()
fontP.set_size('large')
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]
fig = plt.figure(1,figsize=(12,12))
axScatter = fig.add_subplot(223, position=rect_scatter)
plt.xlabel(r'$\eta_{\nu}$', fontsize=28)
plt.ylabel(r'$V_{\nu}$', fontsize=28)
axHistx=fig.add_subplot(221, position=rect_histx)
axHisty=fig.add_subplot(224, position=rect_histy)
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)
axHistx.axes.yaxis.set_ticklabels([])
axHisty.axes.xaxis.set_ticklabels([])

for i in range(len(frequencies)):
    plotdataTMP=plotdata.loc[(plotdata['freq']==frequencies[i])]
    xdata_var=np.log10(plotdataTMP['eta'])
    ydata_var=np.log10(plotdataTMP['V'])
    axScatter.scatter(xdata_var, ydata_var,color=col[i], s=10., zorder=5)

    
x = np.log10(plotdata['eta'])
y = np.log10(plotdata['V'])

axHistx.hist(x, bins=pltvp.make_bins(x), normed=1, histtype='stepfilled', color='b')
axHisty.hist(y, bins=pltvp.make_bins(y), normed=1, histtype='stepfilled', orientation='horizontal', color='b')

freq_labels=[int(f) for f in frequencies]
axScatter.legend(freq_labels,loc=4, prop=fontP)
xmin=int(min(x)-1.1)
xmax=int(max(x)+1.1)
ymin=int(min(y)-1.1)
ymax=int(max(y)+1.1)
xvals=range(xmin,xmax)
xtxts=[r'$10^{'+str(a)+'}$' for a in xvals]
yvals=range(ymin,ymax)
ytxts=[r'$10^{'+str(a)+'}$' for a in yvals]
axScatter.set_xlim([xmin,xmax])
axScatter.set_ylim([ymin,ymax])
axScatter.set_xticks(xvals)
axScatter.set_xticklabels(xtxts, fontsize=20)
axScatter.set_yticks(yvals)
axScatter.set_yticklabels(ytxts, fontsize=20)
axHistx.set_xlim( axScatter.get_xlim())
axHisty.set_ylim( axScatter.get_ylim())

if sigcutx != 0 or sigcuty != 0:
    axHistx.axvline(x=sigcutx, linewidth=2, color='k', linestyle='--')
    axHisty.axhline(y=sigcuty, linewidth=2, color='k', linestyle='--')
    axScatter.axhline(y=sigcuty, linewidth=2, color='k', linestyle='--')
    axScatter.axvline(x=sigcutx, linewidth=2, color='k', linestyle='--')

range_x,fitx = pltvp.gaussian_fit(x,paramx)
axHistx.plot(range_x,fitx, 'k:', linewidth=2)
range_y,fity = pltvp.gaussian_fit(y,paramy)
axHisty.plot(fity,range_y, 'k:', linewidth=2)

#plt.show()
plt.savefig('ds'+str(dataset_id)+'_newsrcs.png')

#tmp=plotdata.loc[(plotdata['eta']>=1.)]# & (plotdata['V']>=10.**sigcuty)]

#tmp=plotdata.loc[(plotdata['eta']>=10.**sigcutx) & (plotdata['V']>=10.**sigcuty)]
tmp=plotdata

IdTrans=tmp.RuncatID.unique()
if len(tmp)>0:
    for a in IdTrans:
        print websiteURL+str(a)
else:
    print "No Transients"

print "Number of new sources after Filter 5: "+str(len(tmp))

    
exit()
