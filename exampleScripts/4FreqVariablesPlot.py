# This script outputs 4 plots of the variability parameters for 4 different observing frequencies.
# Run FilterVariables.py separately for each unique observing frequency first

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

matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['ytick.major.size'] = 10
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['ytick.major.width'] = 2

# Add your observing frequencies here
Freq0='0.96 GHz'
Freq1='1.18 GHz'
Freq2='1.39 GHz'
Freq3='1.61 GHz'

# Add the thresholds obtained from FilterVariables.py
Freq0_eta_cut=8.22
Freq1_eta_cut=16.5
Freq2_eta_cut=13.5
Freq3_eta_cut=50.2
Freq0_V_cut=0.203
Freq1_V_cut=0.196
Freq2_V_cut=0.204
Freq3_V_cut=0.179

# Path to the data table output by FilterVariables.py
Freq0_data=pd.read_csv('freq0_data.csv')
Freq1_data=pd.read_csv('freq1_data.csv')
Freq2_data=pd.read_csv('freq2_data.csv')
Freq3_data=pd.read_csv('freq3_data.csv')

# Creating eta V plot
nullfmt   = NullFormatter()         # no labels
fontP = FontProperties()
fontP.set_size('large')
fig = plt.figure(1,figsize=(12,12))
plt.tight_layout()
fig.subplots_adjust(hspace=0.05,wspace=0.05)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

ax1.xaxis.set_major_formatter(nullfmt)
ax1.axes.xaxis.set_ticklabels([])

ax2.xaxis.set_major_formatter(nullfmt)
ax2.yaxis.set_major_formatter(nullfmt)
ax2.axes.xaxis.set_ticklabels([])
ax2.axes.yaxis.set_ticklabels([])

ax4.yaxis.set_major_formatter(nullfmt)
ax4.axes.yaxis.set_ticklabels([])

ax1.scatter(Freq0_data['eta'],Freq0_data['V'], s=20.)
ax1.axhline(Freq0_V_cut, linewidth=2, color='k', linestyle='--')
ax1.axvline(Freq0_eta_cut, linewidth=2, color='k', linestyle='--')

ax1.text(45,1.56,'1', fontsize=14)
ax1.text(24,0.203,'2', fontsize=14)
ax1.text(490,0.185,'3', fontsize=14)
ax1.text(14,0.0992,'4', fontsize=14)


ax2.scatter(Freq1_data['eta'],Freq1_data['V'], s=20.)
ax2.axhline(Freq1_V_cut, linewidth=2, color='k', linestyle='--')
ax2.axvline(Freq1_eta_cut, linewidth=2, color='k', linestyle='--')

ax2.text(85,1.43,'1', fontsize=14)
ax2.text(31,0.265,'2', fontsize=14)
ax2.text(314,0.201,'3', fontsize=14)
ax2.text(39,0.157,'4', fontsize=14)

ax3.scatter(Freq2_data['eta'],Freq2_data['V'], s=20.)
ax3.axhline(Freq2_V_cut, linewidth=2, color='k', linestyle='--')
ax3.axvline(Freq2_eta_cut, linewidth=2, color='k', linestyle='--')

ax3.text(211,1.66,'1', fontsize=14)
ax3.text(79,0.304,'2', fontsize=14)
ax3.text(249,0.167,'3', fontsize=14)
ax3.text(146,0.222,'4', fontsize=14)

ax4.scatter(Freq3_data['eta'],Freq3_data['V'], s=20.)
ax4.axhline(Freq3_V_cut, linewidth=2, color='k', linestyle='--')
ax4.axvline(Freq3_eta_cut, linewidth=2, color='k', linestyle='--')

ax4.text(109,1.45,'1', fontsize=14)
ax4.text(19,0.410,'2', fontsize=14)
ax4.text(152,0.206,'3', fontsize=14)
ax4.text(41,0.205,'4', fontsize=14)

ax1.text(0.3,1.5,Freq0, fontsize=20)
ax2.text(0.3,1.5,Freq1, fontsize=20)
ax3.text(0.3,1.5,Freq2, fontsize=20)
ax4.text(0.3,1.5,Freq3, fontsize=20)



all_etas=list(Freq0_data['eta'])+list(Freq1_data['eta'])+list(Freq2_data['eta'])+list(Freq3_data['eta'])
all_Vs=list(Freq0_data['V'])+list(Freq1_data['V'])+list(Freq2_data['V'])+list(Freq3_data['V'])

xmin=min(all_etas)*0.8
xmax=max(all_etas)*1.4
ymin=min(all_Vs)*0.8
ymax=max(all_Vs)*1.4

ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin,ymax])

ax1.set_ylabel(r'$V$', fontsize=20)
ax3.set_ylabel(r'$V$', fontsize=20)
ax3.set_xlabel(r'$\eta$', fontsize=20)
ax4.set_xlabel(r'$\eta$', fontsize=20)

ax1.set_yscale('log')
ax1.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xscale('log')

ax1.set_xlim( ax3.get_xlim() )
ax1.set_ylim( ax3.get_ylim() )
ax2.set_xlim( ax3.get_xlim() )
ax2.set_ylim( ax3.get_ylim() )
ax4.set_ylim( ax3.get_ylim() )
ax4.set_xlim( ax3.get_xlim() )
ax1.xaxis.set_major_formatter(nullfmt)
ax4.yaxis.set_major_formatter(nullfmt)
ax2.xaxis.set_major_formatter(nullfmt)
ax2.yaxis.set_major_formatter(nullfmt)

plt.savefig('eta_V_plots.png')
