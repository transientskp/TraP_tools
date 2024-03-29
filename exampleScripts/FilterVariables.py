# FilterVariables.py
# See also FilterVariables.ipynb
#
# A code to plot variability parameters for all sources in a given dataset
# and to select the candidate variables using input sigma thresholds.
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
from database_tools import dbtools
from tools import tools
from plotting import plot_varib_params as pltvp
import matplotlib
import matplotlib.pyplot as plt
import pylab
pylab.rcParams['legend.loc'] = 'best'
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties

# The input database, dataset and thresholds
dataset_id = 9
database = 'AR_R4'
sigma1 = 2 # Threshold on the reduced weighted chi^2
sigma2 = 2 # Threshold on the variability parameter
websiteURL = 'http://banana.transientskp.org/r3/vlo_'+database+'/runningcatalog/'
outname = 'freq3_data.csv' # name of the plotdata datafile

# New inputs to search for duplicates. The BMaj is the restoring beam major axis as given in the fits header of the images processed
BMaj = 0.00194102683426239  # in degrees
beamwidths=5.
                                                        


# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)
VarParams = dbtools.GetVarParams(session,dataset_id)

# Example queries:
#print VarParams[1].Varmetric.runcat.dataset_id
#print VarParams[1].Runningcatalog.id
#print [VarParams[i].Runningcatalog.id for i in range(len(VarParams))]

# Remove duplicates --- note this method is SLOW!!!                                                                                                                                                                                     
Srcs=[]
matchSrcs=[]
SrcAssocRadius = BMaj * beamwidths # in degrees, 10" is ~2.8e-3
for i in range(len(VarParams)):
    VarParams_coord = SkyCoord(ra=(VarParams[i].Runningcatalog.wm_ra*u.degree),dec=(VarParams[i].Runningcatalog.wm_decl*u.degree))
    Srcs.append([VarParams[i].Runningcatalog.id,VarParams_coord])
for a in range(len(Srcs)):
    for b in range(len(Srcs)):
        if b>a:
            if Srcs[b][0] not in matchSrcs:
            sep = Srcs[a][1].separation(Srcs[b][1])
                if sep < SrcAssocRadius*u.degree:
                    matchSrcs.append(Srcs[b][0])


# Set up data for plotting
plotdata = [[VarParams[i].Runningcatalog.id, VarParams[i].Varmetric.eta_int, VarParams[i].Varmetric.v_int, VarParams[i].Varmetric.lightcurve_max, VarParams[i].Varmetric.lightcurve_median, (VarParams[i].Varmetric.band.freq_central/1e6), VarParams[i].Runningcatalog.datapoints, VarParams[i].Varmetric.newsource] for i in range(len(VarParams)) if VarParams[i].Runningcatalog.id not in matchSrcs]
plotdata = pd.DataFrame(data=plotdata,columns=['runcat','eta','V','maxFlx','avgFlx','freq','dpts','newSrc'])
plotdata = plotdata.fillna('N')

# Gaussian fitting to population variability parameters
plotdata = plotdata.loc[(plotdata['eta'] > 0) & (plotdata['V'] > 0) & (plotdata['dpts']>1) & (plotdata['newSrc']=='N')]
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

# Save the data for plotting
plotdata.to_csv(outname, index=False)

# Find the unique frequencies for plotting
frequencies = plotdata.freq.unique()

# Create a colourmap for each of the parameters
col = pltvp.make_cmap(frequencies)

# Creating eta V plot
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
plt.savefig('ds'+str(dataset_id)+'_scatter_plots.png')

plt.clf()

# Create diagnostic plots
fig = plt.figure(1,figsize=(12,12))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
fontP = FontProperties()
fontP.set_size('large')
fig.subplots_adjust(hspace = .001, wspace = 0.001)
ax1.set_ylabel(r'$\eta_\nu$', fontsize=28)
ax3.set_ylabel(r'$V_\nu$', fontsize=28)
ax3.set_xlabel('Max Flux (Jy)', fontsize=24)
ax4.set_xlabel('Max Flux / Median Flux', fontsize=24)

for i in range(len(frequencies)):
    plotdataTMP=plotdata.loc[(plotdata['freq']==frequencies[i])]
    xdata_ax3=plotdataTMP['maxFlx']
    xdata_ax4=plotdataTMP['maxFlx']/plotdataTMP['avgFlx']
    ydata_ax1=plotdataTMP['eta']
    ydata_ax3=plotdataTMP['V']
    ax1.scatter(xdata_ax3, ydata_ax1,color=col[i], s=10., zorder=5)
    ax2.scatter(xdata_ax4, ydata_ax1,color=col[i], s=10., zorder=6)
    ax3.scatter(xdata_ax3, ydata_ax3,color=col[i], s=10., zorder=7)
    ax4.scatter(xdata_ax4, ydata_ax3,color=col[i], s=10., zorder=8)
    ax4.legend(freq_labels, loc=4, prop=fontP)

Xax3=plotdata['maxFlx']
Xax4=plotdata['maxFlx']/plotdataTMP['avgFlx']
Yax1=plotdata['eta']
Yax3=plotdata['V']

    
if sigcutx != 0 or sigcuty != 0:
    ax1.axhline(y=10.**sigcutx, linewidth=2, color='k', linestyle='--')
    ax2.axhline(y=10.**sigcutx, linewidth=2, color='k', linestyle='--')
    ax3.axhline(y=10.**sigcuty, linewidth=2, color='k', linestyle='--')
    ax4.axhline(y=10.**sigcuty, linewidth=2, color='k', linestyle='--')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax4.set_yscale('log')
xmin_ax3=10.**(int(np.log10(min(Xax3))-1.1))
xmax_ax3=10.**(int(np.log10(max(Xax3))+1.2))
xmin_ax4=0.8
xmax_ax4=int(max(xdata_ax4)+0.5)
ymin_ax1=10.**(int(np.log10(min(Yax1))-1.1))
ymax_ax1=10.**(int(np.log10(max(Yax1))+1.2))
ymin_ax3=10.**(int(np.log10(min(Yax3))-1.1))
ymax_ax3=10.**(int(np.log10(max(Yax3))+1.2))
ax1.set_ylim(ymin_ax1,ymax_ax1)
ax3.set_ylim(ymin_ax3,ymax_ax3)
ax3.set_xlim(xmin_ax3,xmax_ax3)
ax4.set_xlim(xmin_ax4,xmax_ax4)
ax1.set_xlim( ax3.get_xlim() )
ax4.set_ylim( ax3.get_ylim() )
ax2.set_xlim( ax4.get_xlim() )
ax2.set_ylim( ax1.get_ylim() )
ax1.xaxis.set_major_formatter(nullfmt)
ax4.yaxis.set_major_formatter(nullfmt)
ax2.xaxis.set_major_formatter(nullfmt)
ax2.yaxis.set_major_formatter(nullfmt)
plt.savefig('ds'+str(dataset_id)+'_diagnostic_plots.png')

# Print out the URL to the Banana webpage for each candidate
tmp=plotdata.loc[(plotdata['eta']>=10.**sigcutx) & (plotdata['V']>=10.**sigcuty)]
tmp.to_csv('candidate_variables.csv',index=False)

IdTrans=tmp.runcat.unique()
if len(tmp)>0:
    for a in IdTrans:
        print websiteURL+str(a)
else:
    print "No Variables"

