import matplotlib 
matplotlib.use('Agg') # There was originally a reason why these two lines were first, but I don't remember why anymore, may not matter now
import numpy as np
import psycopg2
import datetime
import os 
import sys
import matplotlib.pyplot as plt
import glob
import argparse
from multiprocessing import Pool, TimeoutError
from scipy.stats import norm, rankdata
from scipy.special import erf
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import Quadrangle
from astropy import units as u
import imageio
from tqdm import tqdm 
from psycopg2.extensions import register_adapter, AsIs
# psycopg2 complains about numpy datatypes, this avoids that error. Got this from stackoverflow but don't recall where
def addapt_numpy_float64(numpy_float64):
    return AsIs(numpy_float64)
def addapt_numpy_int64(numpy_int64):
    return AsIs(numpy_int64)
register_adapter(np.float64, addapt_numpy_float64)
register_adapter(np.int64, addapt_numpy_int64)

# Do arguments here because I want the default argument values to be stored even if this is called as a library and not main
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dataset', nargs="*",  type=int,help='plot --dataset ??? only')
parser.add_argument('--srcid', nargs="*",  type=int, help='plot --srcid ??? only')
parser.add_argument('--rmsqcmax', nargs=1, type=float, help='Max rmsqc to consider for images')
parser.add_argument('--rmsqcmin', nargs=1, type=float, help='Min rmsqc to consider for images')
parser.add_argument('--angrestrict',  type=float, default=1.0, help='Sets max angular sep from image center for candidate sources')
parser.add_argument('--topV', nargs=1, type=int, default=[20], help='Plot top V sources, if no srcid or dataset specified')
parser.add_argument('--topeta', nargs=1, type=int, default=[20], help='Plot top eta sources, if no srcid or dataset specified')
parser.add_argument('--bottomV', nargs=1, type=int, default=[0], help='Plot bottom V sources, if no srcid or dataset specified')
parser.add_argument('--bottometa', nargs=1, type=int, default=[0], help='Plot bottom eta sources, if no srcid or dataset specified')
parser.add_argument('--skysearch', nargs=3, type=float, help='Search sky for srcs expects: RA (decimal degrees) DEC (decimal degrees) radius (arcsec)')
parser.add_argument('--exclusionsfile', nargs=1, type=str, help='File containing regions on sky to exclude in format: RA (decimal degrees), DEC (decimal degrees), radius (arcseconds)')
parser.add_argument('--makemovies', action='store_true', help='Make movies of all sources (slow, default)')
parser.add_argument('--makeplots', action='store_true', help='Make v-eta plots (fast)')
parser.add_argument('--makeimages', action='store_true', help='Make multi-plot of a time series of images around detection')
parser.add_argument('--usefpk', action='store_true', help='Plot Fpeak rather than Fint')
parser.add_argument('--deduplicate', action='store_true', help='Remove sources within  --beamwidths of a bright source')
parser.add_argument('--beamwidths', type=float, default=5., help='Beamwidths to deduplicate sources')
parser.add_argument('--FDlimits', nargs=1, type=str, help='File containing dataset ids and false detection limits in the format: dataset, FD (Jy) ')
parser.add_argument('--noimage', action='store_true', help='Skip all plots of images')
parser.add_argument('--localimage', action='store_true', help='Correct for images on macbook instead of cluster')
parser.add_argument('--parallel', action='store_true', help='Run in parallel (run on full node only)')
parser.add_argument('--parallelanimatebysrc', action='store_true', help='Run in parallel (run on full node only)')
parser.add_argument('--paralleldataset', action='store_true', help='Run in parallel making movies over datasets (run on full node only)')
parser.add_argument('--bigpicture', action='store_true', help='Make movie of srcs in field')
parser.add_argument('--skipsrcmovies', action='store_true', help='Skip movies of sources')
parser.add_argument('--dumpvarmetric', action='store_true', help='Dump varmetric data into CSV file of name "ds???_varmetric.csv"')
parser.add_argument('--associatesources', action='store_true', help='Associate a provided list of sources (via --srcid) with other sources nearby in other images')
parser.add_argument('--resume', action='store_true', help='Resume a previously failed or aborted run based on filenames')
parser.add_argument('--animateinteresting', action='store_true', help='Determine interesting sources using mean and sigma thresholds and animate them')
parser.add_argument('--vmin',  default=-2e-5, type=float, help='Set vmin for color scaling of images. If negative, include equals sign: --vmin=-2e-5')
parser.add_argument('--vmax', default=1e-3, type=float, help='Set vmax for color scaling of images')
parser.add_argument('--nsigma', default=1, type=float, help='Set nsigma threshold to add to mean')
parser.add_argument('--sigthresh', default=5, type=str, help='Limit to sources detected at --sigthresh at least once')
args = parser.parse_args()
# In order to use MJD
start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)

# PDF of the log10 normal distribution
# https://stats.stackexchange.com/questions/22607/base-10-lognormal-pdf-integrated-over-log10x
def lognormal(x, mu, sigma):
    """https://stats.stackexchange.com/questions/22607/base-10-lognormal-pdf-integrated-over-log10x"""
    return (np.log10(np.e) / (x * sigma * np.sqrt(2*np.pi)) ) * np.exp(- ((np.log10(x) - mu)**2 / (2 * sigma**2)))

# Used mathematica to get this, surprisingly nasty
def lognormalmean(mu, sigma):
    return 10**(mu + 0.5 * sigma**2 * np.log(10))

# Also used mathematica here 
def lognormalstdev(mu, sigma):
    return np.sqrt(-10**(2*mu) * (2**(1 + 0.5 * sigma**2 * np.log(500)) * 5**(0.5 * sigma**2 * np.log(50)) - 10**(sigma**2 * np.log(10)) - np.exp(2 * sigma**2 * np.log(10)**2)))

# Get data from varmetric table 
def dovarmetricquery(sqlstring, sqltuple):
    """take in datasetid and output an array with corresponding varmetric data from TraP"""
    cur.execute("""SELECT DISTINCT varmetric.eta_int, varmetric.v_int, image.dataset, runningcatalog.id, MAX(extractedsource.det_sigma) FROM assocxtrsource 
    JOIN runningcatalog ON runningcatalog.id=assocxtrsource.runcat 
    JOIN extractedsource ON extractedsource.id=assocxtrsource.xtrsrc 
    JOIN image ON image.id=extractedsource.image JOIN varmetric ON 
    varmetric.runcat=runningcatalog.id JOIN runningcatalog_flux ON 
    runningcatalog.id=runningcatalog_flux.id JOIN skyregion ON skyregion.id = image.skyrgn WHERE """+sqlstring,sqltuple )
    varlist = cur.fetchall()
    vardat = np.array(varlist, dtype=[('eta','f8'),('v','f8'),('dataset','i8'),('id','i8'),('det_sigma','f8')])
    return vardat


def queryfordatasets(refinesqlstring,inserttuple):
    """Takes in a string to append to the query and the data to insert into the query, outputs datasets to plot"""
    
    cur.execute("""select distinct image.dataset from assocxtrsource 
    join runningcatalog on runningcatalog.id=assocxtrsource.runcat 
    join extractedsource on extractedsource.id=assocxtrsource.xtrsrc 
    join image on image.id=extractedsource.image join varmetric on 
    varmetric.runcat=runningcatalog.id join runningcatalog_flux on 
    runningcatalog.id=runningcatalog_flux.id join skyregion on skyregion.id = image.skyrgn where  """+refinesqlstring,inserttuple)
    datasetstoanalyze = cur.fetchall()
    # creates an array of tuples, need to make into an array of ints
    return tuple([d[0] for d in datasetstoanalyze])


def getetavtoplot(varmetricdata):
    """Input list of numpy arrays of variability data, outputs eta v values we are interested in plotting"""
    global topVindex
    global topetaindex
    global bottomVindex
    global bottometaindex

    relevantetaV = []
    if topVindex > 0 and not bool(args.skysearch) and not bool(args.srcid):
        for d in np.unique(varmetricdata['dataset']):
            var = varmetricdata[varmetricdata['dataset']==d]
            vsortedvar = np.sort(var, axis=0, order='v')
            for i in range(np.maximum(-topVindex,-len(vsortedvar)),0,1):
                 relevantetaV.append(vsortedvar[i])
    if topetaindex > 0 and not bool(args.skysearch) and not bool(args.srcid):
        for d in np.unique(varmetricdata['dataset']):
            var = varmetricdata[varmetricdata['dataset']==d]
            etasortedvar = np.sort(var, axis=0, order='eta')
            for i in range(np.maximum(-topetaindex,-len(etasortedvar)),0,1):
                 relevantetaV.append(etasortedvar[i])
    if bottomVindex > 0 and not bool(args.skysearch) and not bool(args.srcid):
        for d in np.unique(varmetricdata['dataset']):
            var = varmetricdata[varmetricdata['dataset']==d]
            vsortedvar = np.sort(var[var['v'] > 0], axis=0, order='v')
            for i in range(bottomVindex):
                 relevantetaV.append(vsortedvar[i])
    if bottometaindex > 0 and not bool(args.skysearch) and not bool(args.srcid):
        for d in np.unique(varmetricdata['dataset']):
            var = varmetricdata[varmetricdata['dataset']==d]
            etasortedvar = np.sort(var[var['eta'] > 0], axis=0, order='eta')
            for i in range(bottometaindex):
                 relevantetaV.append(etasortedvar[i])
    if bool(args.skysearch) or bool(args.srcid):
        for var in varmetricdata:
            relevantetaV.append(var)
    return np.unique(np.array(relevantetaV, dtype=[('eta', 'f8'),('v','f8'),('dataset','i8'),('id','i8'),('det_sigma','f8')]))


def querystablesrcs(refinesqlstring, inserttuple):
    """Take in a string to add to query and a tuple to insert, return a dataset of stable sources.
    Where stabarr is made"""
    
    cur.execute("""SELECT image.taustart_ts , extractedsource.f_int, extractedsource.ra, extractedsource.decl, image.url, runningcatalog.id, extractedsource.f_peak_err, extractedsource.f_int_err, extractedsource.f_peak, image.dataset, varmetric.v_int, varmetric.eta_int  FROM assocxtrsource 
    JOIN runningcatalog ON runningcatalog.id=assocxtrsource.runcat 
    JOIN extractedsource ON extractedsource.id=assocxtrsource.xtrsrc 
    JOIN image ON image.id=extractedsource.image JOIN varmetric ON 
    varmetric.runcat=runningcatalog.id JOIN runningcatalog_flux ON 
    runningcatalog.id=runningcatalog_flux.id JOIN skyregion ON skyregion.id = image.skyrgn WHERE (extractedsource.f_int > 0.01) AND (sqlseparation(extractedsource.ra, extractedsource.decl, skyregion.centre_ra, skyregion.centre_decl) < (0.5))  """+refinesqlstring, inserttuple )
    stab = cur.fetchall()
    return np.array(stab, dtype=[('date',object),('f_int','f8'),('ra','f8'),('dec','f8'),('image','U128'),('id','i8'),('f_pk_err','f8'),('f_int_err','f8'),('f_pk','f8'),('dataset','i8'),('v','f8'),('eta','f8')])

def queryallsrcs(refinesqlstring, inserttuple):
    """Take in a string to addd to query and a tuple to insert, and return a dataset of all detected sources
    Where bigarr is made """

    cur.execute("""SELECT image.taustart_ts , extractedsource.f_int, extractedsource.ra, extractedsource.decl, image.url, runningcatalog.id, extractedsource.f_peak_err, extractedsource.f_int_err, extractedsource.f_peak, image.dataset, varmetric.v_int, varmetric.eta_int FROM assocxtrsource 
    JOIN runningcatalog ON runningcatalog.id=assocxtrsource.runcat 
    JOIN extractedsource ON extractedsource.id=assocxtrsource.xtrsrc 
    JOIN image ON image.id=extractedsource.image JOIN varmetric ON 
    varmetric.runcat=runningcatalog.id JOIN runningcatalog_flux ON 
    runningcatalog.id=runningcatalog_flux.id JOIN skyregion ON skyregion.id = image.skyrgn WHERE   """+refinesqlstring, inserttuple)
    big = cur.fetchall()
    return np.array(big, dtype=[('date',object),('f_int','f8'),('ra','f8'),('dec','f8'),('image','U128'),('id','i8'),('f_pk_err','f8'),('f_int_err','f8'),('f_pk','f8'),('dataset','i8'),('v','f8'),('eta','f8')])

def getinterestingsrcs(d):
    v = varmetricdata[varmetricdata['dataset']==d]
    if v.size > 0:
        etavalsnoNaN = np.nan_to_num(v['eta'][v['eta']>0])
        # If we fit the log10 of the data, we can use the parameters to calculate the
        # mean and stdev of the non-logged data, which gives all necessary parameters for
        # the log10 normal distribution
        (etamu, etasigma) = norm.fit(np.log10(etavalsnoNaN))
        histetabins = np.geomspace(np.amin(etavalsnoNaN), np.amax(etavalsnoNaN), num =20)
        etap = lognormal(histetabins,etamu,etasigma)
        mean = lognormalmean(etamu,etasigma)
        stdev = lognormalstdev(etamu, etasigma)
        etaidthresh = v['id'][v['eta'] > (mean+args.nsigma*stdev)] 
        vvalsnoNaN= np.nan_to_num(v['v'][v['v'] > 0])
        (vmu, vsigma) = norm.fit(np.log10(vvalsnoNaN))
        histvbins = np.geomspace(np.amin(vvalsnoNaN), np.amax(vvalsnoNaN), num =20)
        vp = lognormal(histvbins,vmu,vsigma)
        mean = lognormalmean(vmu,vsigma)
        stdev = lognormalstdev(vmu,vsigma)
        vidthresh = v['id'][v['v'] > (mean+args.nsigma*stdev)]
        bothidthresh = np.intersect1d(vidthresh,etaidthresh) 
        return bothidthresh
    else:
        return np.array([]) 
def makeplots(d):   
    """Input a dataset, output a plot"""
    # We aren't animating, so any date within the dataset is fine for now, we'll get more exact later
    if dates[bigarr['dataset']==d].size > 0:
        mydate = dates[bigarr['dataset']==d][0]
        # Hopefully this provides a useful name, could break later, then we need to just remove this
        imagefamily = bigarr['image'][bigarr['dataset']==d][0].rsplit('/',1)[0]
        v = varmetricdata[varmetricdata['dataset']==d]
        # Our sources of interest need to appear in both arrays to be something for us to plot
        allsrcid = np.intersect1d(bigarr['id'],v['id'])
        sortbyid = np.argsort(allsrcid)
        vplot = []
        if args.usefpk:
            for id in allsrcid[sortbyid]:
                vplot.append((np.nan_to_num(v['eta'][v['id']==id][0]), np.nan_to_num(v['v'][id==v['id']][0]), np.amax(bigarr['f_pk'][bigarr['id']==id])))
            vplotarr = np.array(vplot, dtype=[('eta','f8'),('v','f8'),('fpeak','f8')])
        else:
            # Make an array with V, eta, and max flux values, we use max flux here just as is done in https://ui.adsabs.harvard.edu/abs/2019A%26C....27..111R/abstract Rowlinson et al. 2019
            for id in allsrcid[sortbyid]:
                vplot.append((v['eta'][v['id']==id][0], v['v'][id==v['id']][0], np.amax(bigarr['f_int'][bigarr['id']==id])))
            vplotarr = np.array(vplot, dtype=[('eta','f8'),('v','f8'),('fint','f8')])
        fig, axs = plt.subplots(2,3, figsize=(22.5,15), sharex=False, sharey=False)
        fig.suptitle("From images in "+imagefamily)
        axs[0][0].scatter(v['eta'], v['v'], s=15, marker='x', color='navy')
        myline = np.geomspace(np.minimum(vmin,etamin), np.maximum(vmax,etamax),num=100)
       #  mynextline = np.sqrt(myline)/np.average(vplotarr['fint'])
        axs[0][0].plot(myline,myline)
        allerr = np.average(bigarr['f_int_err'][bigarr['dataset']==d])
        avgint = np.average(bigarr['f_int'][bigarr['dataset']==d])
        mynextline = np.sqrt(myline)/(avgint/allerr)
        axs[0][0].plot(myline,mynextline, c='red')
        axs[0][0].set_xscale('log')
        axs[0][0].set_yscale('log')
        axs[0][0].set_xlim(etamin,etamax)
        axs[0][0].set_ylim(vmin,vmax)
        axs[0][0].set_title('dataset '+str(d))
        axs[0][0].set_xlabel('eta')
        axs[0][0].set_ylabel('V')
        
        # Plot Fint stability 
        
        reldates = dates[bigarr['dataset']==v['dataset'][0]]
        mindate = np.amin(reldates)
        maxdate = np.amax(reldates)
        # determine the earliest possible date that exists in the dataset that is within 1/2 day of mydate
        # determine the latest possible date that exists in the dataset that is within 1/2 day of mydate
        # Together they will make up the plot range for x axis 
        beginepoch = np.amin(stabdates[stabdates > (mydate-0.5)])
        endepoch = np.amax(stabdates[stabdates < (mydate+0.5)])
        axs[0][1].scatter(stabdates,stabydat1, marker='.', s=1.5, color='black')
        axs[0][1].set_xlabel('OBS Date (MJD)')
        if args.usefpk:
            axs[0][1].set_ylabel('Peak Flux (Jy)')
            axs[0][1].title.set_text("dataset "+str(d)+' Fpk')
        else:
            axs[0][1].set_ylabel('Integrated Flux (Jy)')
            axs[0][1].title.set_text("dataset "+str(d)+' Fint')
        axs[0][1].set_yscale('log')
        axs[0][1].set_xlim(beginepoch, endepoch)
        axs[0][1].axvspan(mindate, maxdate, alpha=0.5, color='red')
        
        # eta histogram
        # Do a log binning
        etavalsnoNaN = np.nan_to_num(v['eta'][v['eta']>0])
        histetabins = np.geomspace(np.amin(etavalsnoNaN), np.amax(etavalsnoNaN), num =20)
        etahistvals, _, _ = axs[0][2].hist(etavalsnoNaN, bins=histetabins, density=False, alpha=0.6, color='g')
        binwidth = histetabins[1::2] - histetabins[::2]
        # If we fit the log10 of the data, we can use the parameters to calculate the
        # mean and stdev of the non-logged data, which gives all necessary parameters for
        # the log10 normal distribution
        (etamu, etasigma) = norm.fit(np.log10(etavalsnoNaN))
        etap = lognormal(histetabins,etamu,etasigma)
        # need to adjust the height of the pdf to match the histogram since we 
        # do not want a density hisogram. We correct for the bin width by taking 
        # a log spaced correction from the smallest to the largest binwidth, 
        # not sure how correct this is
        etapcorr = [len(etavalsnoNaN)*bin*etapval for bin,etapval in zip(np.geomspace(np.amin(binwidth),np.amax(binwidth),num=len(etap)),etap)]
        axs[0][2].plot(histetabins,etapcorr, 'k', linewidth=2)
        axs[0][2].set_xlabel('eta')
        axs[0][2].set_ylabel('Number of sources')
        title = "Fit results: mu = %.2f,  sigma = %.2f" % (etamu, etasigma)
        axs[0][2].set_title(title)
        axs[0][2].set_xscale('log')
        mean = lognormalmean(etamu,etasigma)
        stdev = lognormalstdev(etamu, etasigma)
        axs[0][2].axvline(mean+args.nsigma*stdev )
        # also plot on eta-v
        axs[0][0].axvline(mean +args.nsigma*stdev)
        etaidthresh = v['id'][v['eta'] > (mean+args.nsigma*stdev)] 
        print("Ids greater than eta threshold: ",etaidthresh)

        # Eta - flux
        if args.usefpk:
            axs[1][0].scatter(vplotarr['fpeak'], vplotarr['eta'], s=15, marker='x', color='navy')
        else:
            axs[1][0].scatter(vplotarr['fint'], vplotarr['eta'], s=15, marker='x', color='navy')
        axs[1][0].set_xscale('log')
        axs[1][0].set_yscale('log')
        #  axs[1][0].set_xlim(xmin,xmax)
        axs[1][0].set_ylim(etamin,etamax)
        axs[1][0].set_title('dataset'+str(d))
        if args.usefpk:
            axs[1][0].set_xlabel('F_peak')
        else:
            axs[1][0].set_xlabel('F_int')
        axs[1][0].set_ylabel('eta')
        axs[1][0].axhline(mean +args.nsigma*stdev)

        # V - flux
        if args.usefpk:
            axs[1][1].scatter(vplotarr['fpeak'], vplotarr['v'], s=15, marker='x', color='navy')
            axs[1][1].set_xlabel('F_peak')
        else:
            axs[1][1].scatter(vplotarr['fint'], vplotarr['v'], s=15, marker='x', color='navy')
            axs[1][1].set_xlabel('F_int')
        axs[1][1].set_xscale('log')
        axs[1][1].set_yscale('log')
        #  axs[1][1].set_xlim(xmin,xmax)
        axs[1][1].set_ylim(vmin,vmax)
        axs[1][1].set_title('dataset'+str(d))
        axs[1][1].set_ylabel('V')

        # v histogram
        # Do a log binning
        vvalsnoNaN= np.nan_to_num(v['v'][v['v'] > 0])
        histvbins = np.geomspace(np.amin(vvalsnoNaN), np.amax(vvalsnoNaN), num =20)
        vhistvals, _, _ = axs[1][2].hist(vvalsnoNaN, bins=histvbins, density=False, alpha=0.6, color='g')
        binwidth = histvbins[1::2] - histvbins[::2]
        # If we fit the log10 of the data, we can use the parameters to calculate the
        # mean and stdev of the non-logged data, which gives all necessary parameters for
        # the log10 normal distribution
        (vmu, vsigma) = norm.fit(np.log10(vvalsnoNaN))
        vp = lognormal(histvbins,vmu,vsigma)
        # need to adjust the height of the pdf to match the histogram since we 
        # do not want a density hisogram. We correct for the bin width by taking 
        # a log spaced correction from the smallest to the largest binwidth, 
        # not sure how correct this is
        vpcorr = [len(vvalsnoNaN)*bin*vpval for bin,vpval in zip(np.geomspace(np.amin(binwidth),np.amax(binwidth),num=len(vp)),vp)]
        axs[1][2].plot(histvbins, vpcorr, 'k', linewidth=2)
        axs[1][2].set_xlabel('V')
        axs[1][2].set_ylabel('Number of sources')
        axs[1][2].set_xscale('log')
        title = "Fit results: mu = %.2f,  sigma = %.2f" % (vmu, vsigma)
        axs[1][2].set_title(title)
        mean = lognormalmean(vmu,vsigma)
        stdev = lognormalstdev(vmu,vsigma)
        axs[1][2].axvline(mean+args.nsigma*stdev)
        # also plot on eta-v
        axs[0][0].axhline(mean+args.nsigma*stdev)
        axs[1][1].axhline(mean+args.nsigma*stdev)
        vidthresh = v['id'][v['v'] > (mean+args.nsigma*stdev)]
        print("Ids greater than v threshold: ",vidthresh)
        bothidthresh = np.intersect1d(vidthresh,etaidthresh) 
        print("Ids greater than both thresholds: ",bothidthresh)
        np.savetxt("ds"+str(d)+"interestingsrcs.csv",bothidthresh,delimiter=',',fmt='%i')
        print('Saved interesting sources to ',"ds"+str(d)+"interestingsrcs.csv")
        plt.savefig('dataset'+str(d)+'scatter.png')
        plt.close()
        print("plotted "+'dataset'+str(d)+'scatter.png')
            # Some integration images aren't finished running though TraP yet. These will show an
            # IndexError. We skip these and delete the axis. A simple way of dealing with a problem
            # that I don't want to worry about


def makeimages(src):
    """Plots the specified source and animates it, may be in parallel"""

    print('MY SRC: ',src)
    rel = reletaVarr[reletaVarr['id']==src]
    datasetkey = (bigarr['dataset'] == v['dataset'][0])
    stabdatasetkey = (stabarr['dataset'] == v['dataset'][0])
    srcidscondition = (bigarr['id'] == src)
    files = np.unique(stabarr['image'][stabarr['dataset']==d])
    if args.usefpk:
        stabydat1=stabarr['f_pk']
    else:
        stabydat1 = stabarr['f_int']
    srcid = bigarr['id'][srcidscondition][0]
    bigdates = np.array([(b['date'] - start_epoch).total_seconds()/3600/24 for b in bigarr])
    stabdates = np.array([(s['date'] - start_epoch).total_seconds()/3600/24 for s in stabarr])
    sortedstabdatetimearr = np.sort(np.unique(stabarr[['image','date']]),order='date')
    plotfilenames = []
    detloc = SkyCoord(bigarr['ra'][srcidscondition][0], bigarr['dec'][srcidscondition][0], unit='deg')
    dates = []
    dummyint = 0
    detectiondate = np.amin(bigarr[srcidscondition]['date'])
    detectionMJD = (detectiondate - start_epoch).total_seconds()/3600/24
    detectionfile = bigarr['image'][bigarr['date']==detectiondate]
    mindateind = max(0,np.where(detectiondate==sortedstabdatetimearr['date'])[0][0] - 2 )
    maxdateind = min(len(sortedstabdatetimearr)-1, np.where(detectiondate==sortedstabdatetimearr['date'])[0][0] + 2)
    fig = plt.figure(figsize=(25,15))
    axlist = []
    numims = maxdateind - mindateind + 1
    for i, datefile in zip(range(numims), sortedstabdatetimearr[mindateind:(maxdateind+1)]):
        myfile = datefile['image']
        if args.localimage:
            myfile = '/Users/schastain' + myfile
        hdul = fits.open(myfile)
        hdu = hdul[0]
        wcs = WCS(hdu.header, naxis=2)
        axlist.append(fig.add_subplot(2,3,i+1,projection= wcs))
        wcslocation = wcs.world_to_pixel(detloc)
        axlist[i].set_title(str(i+1))
        axlist[i].imshow(hdu.data[0:][0:][0][0], vmin=args.vmin, vmax=args.vmax, origin='lower')
        axlist[i].set_ylim(wcslocation[1] - 50, wcslocation[1] + 50 )
        axlist[i].set_xlim(wcslocation[0] - 50 , wcslocation[0] + 50)
        if datefile['date']==detectiondate:
            axlist[i].scatter(wcslocation[0], wcslocation[1], color='red', marker='x')
        else:
            axlist[i].scatter(wcslocation[0], wcslocation[1], color='black', marker='x')
        hdul.close()
    
    axlist.append(fig.add_subplot(2,3,6,projection='rectilinear'))
    ax2xmin = np.amin(stabdates[stabdatasetkey])
    ax2xmax = np.amax(stabdates[stabdatasetkey])
    ax2ymin = 1e-7
    if args.usefpk:
        ax2ymax = np.maximum(1e-2,np.amax(bigarr['f_pk'][srcidscondition] + bigarr['f_pk_err'][srcidscondition]))
        axlist[numims].scatter(bigdates[srcidscondition],bigarr['f_pk'][srcidscondition], marker='.', s=10)
        axlist[numims].errorbar(bigdates[srcidscondition],bigarr['f_pk'][srcidscondition], yerr=bigarr['f_pk_err'][srcidscondition], fmt='none')
        axlist[numims].set_ylabel('Peak Flux (Jy)')
    else:
        ax2ymax = np.maximum(1e-2,np.amax(bigarr['f_int'][srcidscondition] + bigarr['f_int_err'][srcidscondition]))
        axlist[numims].scatter(bigdates[srcidscondition],bigarr['f_int'][srcidscondition], marker='.', s=10)
        axlist[numims].errorbar(bigdates[srcidscondition],bigarr['f_int'][srcidscondition], yerr=bigarr['f_int_err'][srcidscondition], fmt='none')
        axlist[numims].set_ylabel('Integrated Flux (Jy)')
    axlist[numims].set_xlabel('OBS Date (MJD)')
    axlist[numims].set_ylim(ax2ymin, ax2ymax)
    axlist[numims].set_xlim(ax2xmin, ax2xmax)
    axlist[numims].set_yscale('log')
    axlist[numims].axvline(detectionMJD, linestyle=':')
    name = 'src'+str(srcid)+'_multi.png'
    plt.savefig(name)
    print("plotted "+name)
    plt.close()

def dsparallelanimatesrc(d):
    """Plots the specified source and animates it, may be in parallel"""
    fdl = np.unique(fdlimits['fdlimit'][fdlimits['dataset']==d])
    v = varmetricdata[varmetricdata['dataset']==d]
    alreadyplotted = []
    sortedv = np.sort(v['v'])
    sortedeta = np.sort(v['eta'])
    intersectbigrel = np.intersect1d(bigarr['id'],reletaVarr['id'])
    srcstoconsider = np.intersect1d(intersectbigrel,v['id'])
    if srcstoconsider.size==0:
        return
    for src in srcstoconsider:
        rel = reletaVarr[reletaVarr['id']==src]
        datasetkey = (bigarr['dataset'] == v['dataset'][0])
        stabdatasetkey = (stabarr['dataset'] == v['dataset'][0])
        srcidscondition = (bigarr['id'] == src)
        files = np.unique(stabarr['image'][stabarr['dataset']==d])
        if args.usefpk:
            stabydat1=stabarr['f_pk']
        else:
            stabydat1 = stabarr['f_int']
        srcid = bigarr['id'][srcidscondition][0]
        bigdates = np.array([(b['date'] - start_epoch).total_seconds()/3600/24 for b in bigarr])
        stabdates = np.array([(s['date'] - start_epoch).total_seconds()/3600/24 for s in stabarr])
        plotfilenames = []
        detloc = SkyCoord(bigarr['ra'][srcidscondition][0], bigarr['dec'][srcidscondition][0], unit='deg')
        dates = []
        dummyint = 0
        if srcid not in alreadyplotted:
            for myfile in files:
                if args.localimage:
                    myfile = '/Users/schastain' + myfile
                fig = plt.figure(figsize=(25,15))
                if not args.noimage:
                    hdul = fits.open(myfile)
                    hdu = hdul[0]
                    wcs = WCS(hdu.header, naxis=2)
                    wcslocation = wcs.world_to_pixel(detloc)
                    mydate = (datetime.datetime.strptime(hdu.header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") - start_epoch).total_seconds()/3600/24 
                else:
                    mydate = stabdates[dummyint] 
                    dummyint+=1
                dates.append(mydate) 
                # Show movie of source location
                if not args.noimage:
                    ax1 = fig.add_subplot(2,3,1,projection=wcs)
                    ax1.imshow(hdu.data[0:][0:][0][0], vmin=args.vmin, vmax=args.vmax, origin='lower')
                    ax1.set_ylim(wcslocation[1] - 50, wcslocation[1] + 50 )
                    ax1.set_xlim(wcslocation[0] - 50 , wcslocation[0] + 50)
                    if bigarr['ra'][srcidscondition][bigdates[srcidscondition]==mydate].size > 0:
                        skysrc = SkyCoord(bigarr['ra'][srcidscondition][bigdates[srcidscondition]==mydate], bigarr['dec'][srcidscondition][bigdates[srcidscondition]==mydate], unit='deg')
                        wcscurloc = wcs.world_to_pixel(skysrc)
                        ax1.scatter(wcscurloc[0], wcscurloc[1], color='red', marker='x')
                else:
                    ax1 = fig.add_subplot(2,3,1)
                fig.suptitle('Dataset '+str(d)+', Date: '+str(mydate))
               
               # Plot lightcurve  
                ax2 = fig.add_subplot(2,3,2, projection='rectilinear')
                ax2xmin = np.amin(stabdates[stabdatasetkey])
                ax2xmax = np.amax(stabdates[stabdatasetkey])
                ax2ymin = 1e-7
                if args.usefpk:
                    ax2ymax = np.maximum(2,np.amax(bigarr['f_pk'][srcidscondition] + bigarr['f_pk_err'][srcidscondition]))
                    ax2.scatter(bigdates[srcidscondition],bigarr['f_pk'][srcidscondition], marker='.', s=10)
                    ax2.errorbar(bigdates[srcidscondition],bigarr['f_pk'][srcidscondition], yerr=bigarr['f_pk_err'][srcidscondition], fmt='none')
                    ax2.set_ylabel('Peak Flux (Jy)')
                else:
                    ax2ymax = np.maximum(2,np.amax(bigarr['f_int'][srcidscondition] + bigarr['f_int_err'][srcidscondition]))
                    ax2.scatter(bigdates[srcidscondition],bigarr['f_int'][srcidscondition], marker='.', s=10)
                    ax2.errorbar(bigdates[srcidscondition],bigarr['f_int'][srcidscondition], yerr=bigarr['f_int_err'][srcidscondition], fmt='none')
                    ax2.set_ylabel('Integrated Flux (Jy)')
                ax2.axvline(mydate, linestyle=':')
                ax2.set_xlabel('OBS Date (MJD)')
                ax2.set_ylim(ax2ymin, ax2ymax)
                ax2.set_xlim(ax2xmin, ax2xmax)
                ax2.set_yscale('log')
                if fdl.size > 0:
                    ax2.axhline(fdl, linestyle=':')

                # fit variability metrics to lognormal dist
                (etamu, etasigma) = norm.fit(np.log10(v['eta'][v['eta']>0]))
                etamean = lognormalmean(etamu,etasigma)
                etastdev = lognormalstdev(etamu, etasigma)
                (vmu, vsigma) = norm.fit(np.log10(v['v'][v['v']>0]))
                vmean = lognormalmean(vmu,vsigma)
                vstdev = lognormalstdev(vmu, vsigma)

                #  Plot eta v 
                ax3 = fig.add_subplot(2,3,3, projection='rectilinear')
                ax3.scatter(v['eta'], v['v'], s=10, color='navy', marker='x')
                ax3.scatter(rel['eta'], rel['v'], s=25, color='red', marker='X')
                ax3.set_xscale('log')
                ax3.set_yscale('log')
                ax3.set_xlim(xmin,xmax)
                ax3.set_ylim(ymin,ymax)
                ax3.set_title('V: '+str(rel['v'])+', eta: '+str(rel['eta']))
                ax3.set_xlabel('eta')
                ax3.set_ylabel('V')
                ax3.axvline(etamean + args.nsigma*etastdev)
                ax3.axhline(vmean + args.nsigma*vstdev)

                # Plot Fint stability 
                # get dates of the chunk being plotted so that we can highlight the chunk
                reldates = stabdates[stabarr['dataset']==v['dataset'][0]]
                mindate = np.amin(reldates)
                maxdate = np.amax(reldates)
                # determine the earliest possible date that exists in the dataset that is within 1/2 day of mydate
                # determine the latest possible date that exists in the dataset that is within 1/2 day of mydate
                # Together they will make up the plot range for x axis 
                beginepoch = np.amin(stabdates[stabdates > (mydate-0.5)])
                endepoch = np.amax(stabdates[stabdates < (mydate+0.5)])
                ax4 = fig.add_subplot(2,3,4, projection='rectilinear')
                ax4.scatter(stabdates,stabydat1, marker='.', s=1.5, color='black')
                if not args.noimage:
                    ax4.title.set_text(hdu.header['OBJECT']+' Fint')
                ax4.set_xlabel('OBS Date (MJD)')
                if args.usefpk:
                    ax4.set_ylabel('Peak Flux (Jy)')
                else:
                    ax4.set_ylabel('Integrated Flux (Jy)')
                ax4.set_yscale('log')
                ax4.set_xlim(beginepoch, endepoch)
                ax4.axvspan(mindate, maxdate, alpha=0.5, color='red')

                # Show big image
                if not args.noimage:
                    ax5 = fig.add_subplot(2,3,5, projection=wcs)
                    ax5.imshow(hdu.data[0:][0:][0][0], vmin=args.vmin, vmax=args.vmax, origin='lower')
                    ax5.add_patch(Quadrangle((detloc.ra, detloc.dec), 100*u.arcsec, 100*u.arcsec,
                        edgecolor='white', facecolor='none',transform=ax5.get_transform('fk5')))
                else:
                    ax5 = fig.add_subplot(2,3,5)
                
                # Put text in empty plot
                ax6 = fig.add_subplot(2,3,6, projection='rectilinear')
                ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
                rankedV = rankdata(-v['v'], method='ordinal')
                rankedeta = rankdata(-v['eta'], method='ordinal')
                ax6.text(0.1,0.5,ordinal(rankedeta[rel['eta'] == v['eta']][0])+' highest '+'eta', fontsize='large')
                ax6.text(0.1,0.3,ordinal(rankedV[rel['v'] == v['v']][0])+' highest '+'V', fontsize='large')
                ax6.set_xticks([])
                ax6.set_yticks([])
                # Name and save plots
                name = 'variablesANDtransients/src'+str(srcid)+"_"+myfile.rsplit('/')[-1].replace('.fits','')+'chunk'+str(d)+'scatter.png'
                plotfilenames.append(name)
                plt.savefig(name)
                print("plotted "+name)
                if not args.noimage:
                    hdul.close()
                plt.close()
            plotfilenamearr = np.array(plotfilenames)
            datesort = np.argsort(dates)
            with imageio.get_writer('src'+str(srcid)+'.mp4', mode='I',fps=3) as writer:
                for myfile in plotfilenamearr[datesort]:
                    image = imageio.imread(myfile)
                    writer.append_data(image)
            print("wrote video: ",'src'+str(srcid)+'.mp4')
            alreadyplotted.append(srcid)
            for myfile in plotfilenames:
                print('removing: ',myfile)
                os.remove(myfile)


def animatesrc(src):
    """Plots the specified source and animates it, may be in parallel"""
    
    def makeimagesforanimation(myfile):
                return plotfilename
            
    rel = reletaVarr[reletaVarr['id']==src]
    datasetkey = (bigarr['dataset'] == v['dataset'][0])
    stabdatasetkey = (stabarr['dataset'] == v['dataset'][0])
    srcidscondition = (bigarr['id'] == src)
    files = np.unique(stabarr['image'][stabarr['dataset']==d])
    if args.usefpk:
        stabydat1=stabarr['f_pk']
    else:
        stabydat1 = stabarr['f_int']
    srcid = bigarr['id'][srcidscondition][0]
    bigdates = np.array([(b['date'] - start_epoch).total_seconds()/3600/24 for b in bigarr])
    stabdates = np.array([(s['date'] - start_epoch).total_seconds()/3600/24 for s in stabarr])
    plotfilenames = []
    detloc = SkyCoord(bigarr['ra'][srcidscondition][0], bigarr['dec'][srcidscondition][0], unit='deg')
    dates = []
    dummyint = 0
    if srcid not in alreadyplotted:
        for myfile in files:
            if args.localimage:
                myfile = '/Users/schastain' + myfile
            fig = plt.figure(figsize=(25,15))
            if not args.noimage:
                hdul = fits.open(myfile)
                hdu = hdul[0]
                wcs = WCS(hdu.header, naxis=2)
                wcslocation = wcs.world_to_pixel(detloc)
                mydate = (datetime.datetime.strptime(hdu.header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") - start_epoch).total_seconds()/3600/24 
            else:
                mydate = stabdates[dummyint] 
                dummyint+=1
            dates.append(mydate) 
            # Show movie of source location
            if not args.noimage:
                ax1 = fig.add_subplot(2,3,1,projection=wcs)
                ax1.imshow(hdu.data[0:][0:][0][0], vmin=args.vmin, vmax=args.vmax, origin='lower')
                ax1.set_ylim(wcslocation[1] - 50, wcslocation[1] + 50 )
                ax1.set_xlim(wcslocation[0] - 50 , wcslocation[0] + 50)
                if bigarr['ra'][srcidscondition][bigdates[srcidscondition]==mydate].size > 0:
                    skysrc = SkyCoord(bigarr['ra'][srcidscondition][bigdates[srcidscondition]==mydate], bigarr['dec'][srcidscondition][bigdates[srcidscondition]==mydate], unit='deg')
                    wcscurloc = wcs.world_to_pixel(skysrc)
                    ax1.scatter(wcscurloc[0], wcscurloc[1], color='red', marker='x')
            else:
                ax1 = fig.add_subplot(2,3,1)
            fig.suptitle('Dataset '+str(d)+', Date: '+str(mydate))

            #  Plot lightcurve  
            # get dates of the chunk being plotted so that we can highlight the chunk
            reldates = stabdates[stabarr['dataset']==d]
            mindate = np.amin(reldates)
            maxdate = np.amax(reldates)
            # determine the earliest possible date that exists in the dataset that is within 1/2 day of mydate
            # determine the latest possible date that exists in the dataset that is within 1/2 day of mydate
            # Together they will make up the plot range for x axis 
            beginepoch = np.amin(stabdates[stabdates > (mydate-0.5)])
            endepoch = np.amax(stabdates[stabdates < (mydate+0.5)])

            ax2 = fig.add_subplot(2,3,2, projection='rectilinear')
            ax2xmin = np.amin(stabdates[stabdatasetkey])
            ax2xmax = np.amax(stabdates[stabdatasetkey])
            ax2ymin = 1e-7
            if args.usefpk:
                ax2ymax = np.maximum(2,np.amax(bigarr['f_pk'][srcidscondition] + bigarr['f_pk_err'][srcidscondition]))
                ax2.scatter(bigdates[srcidscondition],bigarr['f_pk'][srcidscondition], marker='.', s=10)
                ax2.errorbar(bigdates[srcidscondition],bigarr['f_pk'][srcidscondition], yerr=bigarr['f_pk_err'][srcidscondition], fmt='none')
                ax2.set_ylabel('Peak Flux (Jy)')
            else:
                ax2ymax = np.maximum(2,np.amax(bigarr['f_int'][srcidscondition] + bigarr['f_int_err'][srcidscondition]))
                ax2.scatter(bigdates[srcidscondition],bigarr['f_int'][srcidscondition], marker='.', s=10)
                ax2.errorbar(bigdates[srcidscondition],bigarr['f_int'][srcidscondition], yerr=bigarr['f_int_err'][srcidscondition], fmt='none')
                ax2.set_ylabel('Integrated Flux (Jy)')
            ax2.axvline(mydate, linestyle=':')
            ax2.set_xlabel('OBS Date (MJD)')
            ax2.set_ylim(ax2ymin, ax2ymax)
            ax2.set_xlim(ax2xmin, ax2xmax)
            ax2.set_yscale('log')
            if fdl.size > 0:
                ax2.axhline(fdl, linestyle=':')

            # fit variability metrics to lognormal dist
            (etamu, etasigma) = norm.fit(np.log10(v['eta'][v['eta']>0]))
            etamean = lognormalmean(etamu,etasigma)
            etastdev = lognormalstdev(etamu, etasigma)
            (vmu, vsigma) = norm.fit(np.log10(v['v'][v['v']>0]))
            vmean = lognormalmean(vmu,vsigma)
            vstdev = lognormalstdev(vmu, vsigma)

            #  Plot eta v 
            ax3 = fig.add_subplot(2,3,3, projection='rectilinear')
            ax3.scatter(v['eta'], v['v'], s=10, color='navy', marker='x')
            ax3.scatter(rel['eta'], rel['v'], s=25, color='red', marker='X')
            ax3.set_xscale('log')
            ax3.set_yscale('log')
            ax3.set_xlim(xmin,xmax)
            ax3.set_ylim(ymin,ymax)
            ax3.set_title('V: '+str(rel['v'])+', eta: '+str(rel['eta']))
            ax3.set_xlabel('eta')
            ax3.set_ylabel('V')
            ax3.axvline(etamean + args.nsigma*etastdev)
            ax3.axhline(vmean + args.nsigma*vstdev)

            # Plot Fint stability 
            ax4 = fig.add_subplot(2,3,4, projection='rectilinear')
            ax4.scatter(stabdates,stabydat1, marker='.', s=1.5, color='black')
            if not args.noimage:
                ax4.title.set_text(hdu.header['OBJECT']+' Fint')
            ax4.set_xlabel('OBS Date (MJD)')
            if args.usefpk:
                ax4.set_ylabel('Peak Flux (Jy)')
            else:
                ax4.set_ylabel('Integrated Flux (Jy)')
            ax4.set_yscale('log')
            ax4.set_xlim(beginepoch, endepoch)
            ax4.axvspan(mindate, maxdate, alpha=0.5, color='red')

            # Show big image
            if not args.noimage:
                ax5 = fig.add_subplot(2,3,5, projection=wcs)
                ax5.imshow(hdu.data[0:][0:][0][0], vmin=args.vmin, vmax=args.vmax, origin='lower')
                ax5.add_patch(Quadrangle((detloc.ra, detloc.dec), 100*u.arcsec, 100*u.arcsec,
                    edgecolor='white', facecolor='none',transform=ax5.get_transform('fk5')))
            else:
                ax5 = fig.add_subplot(2,3,5)

            # Put text in empty plot
            ax6 = fig.add_subplot(2,3,6, projection='rectilinear')
            ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
            rankedV = rankdata(-v['v'], method='ordinal')
            rankedeta = rankdata(-v['eta'], method='ordinal')
            ax6.text(0.1,0.5,ordinal(rankedeta[rel['eta'] == v['eta']][0])+' highest '+'eta', fontsize='large')
            ax6.text(0.1,0.3,ordinal(rankedV[rel['v'] == v['v']][0])+' highest '+'V', fontsize='large')
            ax6.set_xticks([])
            ax6.set_yticks([])
            # Name and save plots
            name = 'variablesANDtransients/src'+str(srcid)+"_"+myfile.rsplit('/')[-1].replace('.fits','')+'chunk'+str(d)+'scatter.png'
            plt.savefig(name)
            print("plotted "+name)
            if not args.noimage:
                hdul.close()
            plt.close()
            plotfilenames.append(name)
        plotfilenamearr = np.array(plotfilenames)
        datesort = np.argsort(dates)
        with imageio.get_writer('src'+str(srcid)+'.mp4', mode='I',fps=2) as writer:
            for myfile in plotfilenamearr[datesort]:
                image = imageio.imread(myfile)
                writer.append_data(image)
        print("wrote video: ",'src'+str(srcid)+'.mp4')
        alreadyplotted.append(srcid)
        for myfile in plotfilenames:
            print('removing: ',myfile)
            os.remove(myfile)
                    
def animatefield(srcstoconsider):
    """Plots the specified source and animates it, may be in parallel"""

    bigdates = np.array([(b['date'] - start_epoch).total_seconds()/3600/24 for b in bigarr])
    files = np.unique(stabarr['image'][stabarr['dataset']==d])
    dates = [] 
    plotfilenames = []
    for myfile in files:
        if args.localimage:
            myfile = '/Users/schastain'+myfile
        fig = plt.figure(figsize=(15,15))
        hdul = fits.open(myfile)
        hdu = hdul[0]
        wcs = WCS(hdu.header, naxis=2)
        mydate = (datetime.datetime.strptime(hdu.header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") - start_epoch).total_seconds()/3600/24 
        dates.append(mydate)
        ax1 = fig.add_subplot(1,1,1,projection=wcs)
        ax1.imshow(hdu.data[0:][0:][0][0], vmin=args.vmin, vmax=args.vmax, origin='lower')
        for src in srcstoconsider:
            srcidscondition = (bigarr['id'] == src)
            detloc = SkyCoord(bigarr['ra'][srcidscondition][0], bigarr['dec'][srcidscondition][0], unit='deg')
            if bigarr['ra'][srcidscondition][bigdates[srcidscondition]==mydate].size > 0:
                skysrc = SkyCoord(bigarr['ra'][srcidscondition][bigdates[srcidscondition]==mydate], bigarr['dec'][srcidscondition][bigdates[srcidscondition]==mydate], unit='deg')
                wcscurloc = wcs.world_to_pixel(skysrc)
                ax1.add_patch(Quadrangle((detloc.ra, detloc.dec), 100*u.arcsec, 100*u.arcsec,
                    edgecolor='white', facecolor='none',transform=ax1.get_transform('fk5')))
                    # Name and save plots
        name = 'variablesANDtransients/ds'+str(d)+"_"+myfile.rsplit('/')[-1].replace('.fits','')+'bigfield.png'
        plotfilenames.append(name)
        plt.savefig(name)
        print("plotted "+name)
        hdul.close()
        plt.close()
    plotfilenamearr = np.array(plotfilenames)
    datesort = np.argsort(dates)
    with imageio.get_writer('ds'+str(d)+'_bigfield.mp4', mode='I') as writer:
        for myfile in plotfilenamearr[datesort]:
            image = imageio.imread(myfile)
            writer.append_data(image)
    print("wrote video: ",'ds'+str(d)+'_bigfield.mp4')
    for myfile in plotfilenames:
        print('removing: ',myfile)
        os.remove(myfile)


if __name__ == '__main__':



    pathprefix = 'variablesANDtransients' 
    # Check whether the specified path exists or not
    isExist = os.path.exists(pathprefix)
    if not isExist:
    # Create a new directory because it does not exist 
        os.makedirs(pathprefix)

    host = os.getenv('TKP_DBHOST')
    port = os.getenv('TKP_DBPORT')
    user = os.getenv('TKP_DBUSER')
    password = os.getenv('TKP_DBPASSWORD')
    database = os.getenv('TKP_DBNAME')
    if password is None:
        conn = psycopg2.connect("dbname="+database+" user="+user+" host="+host)
    else:
        conn = psycopg2.connect("dbname="+database+" user="+user+" password="+password+" host="+host)
        
    ################################################################################
    # Typical psycopg2 setup 
    conn = psycopg2.connect("dbname="+database+" user="+user+" password="+password+" host="+host)
    cur = conn.cursor()
    # All the integration runs of interest in this database are greater than this dataset id 
    mindatasetid = "171"

    # Create a function in sql to calculate accurate angular separation given ra and dec of two locations 
    # on the sky 
    try: 
        cur.execute("""CREATE FUNCTION sqlseparation(raa float, decla float, rab float, declb float) RETURNS float
            AS 'select 180.*acos((sin(decla*PI()/180.)*sin(declb*PI()/180.)) +  (cos(decla*PI()/180.)*cos(declb*PI()/180.)*cos((raa*PI()/180.)-(rab*PI()/180.))))/PI() ;'
            LANGUAGE SQL 
            IMMUTABLE
            RETURNS NULL ON NULL INPUT;""",)
    except psycopg2.errors.DuplicateFunction:
        cur.execute("ROLLBACK")
        print("sqlseparation already exists. Skipping.")





    # Get pointings from database in order to establish targets

    
    if bool(args.FDlimits):
        fdlimits = np.loadtxt(args.FDlimits[0], delimiter=',', dtype=[('dataset',int),('fdlimit','f8')], converters = {1: lambda s: float(s or 0)})
    else:
        fdlimits = np.zeros(2,dtype=[('dataset',int),('fdlimit','f8')])
    searchlocation = np.array(args.skysearch)
    topVindex = args.topV[0]
    topetaindex = args.topeta[0]
    bottomVindex = args.bottomV[0]
    bottometaindex = args.bottometa[0]
    if bool(args.exclusionsfile):
        exclusiondata = np.loadtxt(args.exclusionsfile[0], delimiter=',', dtype=[('ra','f8'),('dec','f8'),('radius','f8')])
        excludedid = []
        for e in exclusiondata:
            cur.execute("""SELECT runningcatalog.id FROM runningcatalog WHERE sqlseparation(runningcatalog.wm_ra, runningcatalog.wm_decl, %s, %s) < %s;""", (e['ra'], e['dec'], e['radius']/3600,))
            tmpexclude = cur.fetchall()
            excludedid.extend([element for tup in tmpexclude for element in tup])
        
    if bool(args.rmsqcmax):
        rmsqcmax = args.rmsqcmax[0]
    else:
        rmsqcmax = args.rmsqcmax
    if bool(args.rmsqcmin):
        rmsqcmin = args.rmsqcmin[0]
    else:
        rmsqcmin = args.rmsqcmin
    # Begin a bunch of sql calls to get candidates of interest
    if bool(args.dataset):
        datasetTOexamine = tuple(args.dataset,)
        cur.execute(""" SELECT skyregion.centre_ra, skyregion.centre_decl  FROM skyregion WHERE skyregion.dataset IN %s;""", (datasetTOexamine,))
    elif bool(args.srcid):
        srcidlist = []
        assocsrcs = []
        if args.associatesources:
            print("Associating Sources")
            for src in args.srcid:
                cur.execute("""SELECT DISTINCT runningcatalog.id FROM runningcatalog
                    JOIN image ON image.dataset=runningcatalog.dataset WHERE (runningcatalog.wm_ra < (SELECT runningcatalog.wm_ra+2*image.rb_smaj 
                    FROM image JOIN runningcatalog ON runningcatalog.dataset=image.dataset WHERE runningcatalog.id=%s LIMIT 1)) AND 
                    (runningcatalog.wm_ra > (SELECT runningcatalog.wm_ra-2*image.rb_smaj FROM image JOIN 
                    runningcatalog ON runningcatalog.dataset=image.dataset WHERE runningcatalog.id=%s LIMIT 1)) AND  
                    (runningcatalog.wm_decl < (SELECT runningcatalog.wm_decl+2*image.rb_smaj FROM image JOIN runningcatalog ON 
                    runningcatalog.dataset=image.dataset WHERE runningcatalog.id=%s LIMIT 1)) AND (runningcatalog.wm_decl >
                    (SELECT runningcatalog.wm_decl-2*image.rb_smaj FROM image JOIN runningcatalog ON runningcatalog.dataset=image.dataset
                    WHERE runningcatalog.id=%s LIMIT 1));""",(src,src,src,src,))
                srcappend = cur.fetchall()
                srcidlist.extend([element for tup in srcappend for element in tup])
                assocsrcs.append([element for tup in srcappend for element in tup])
            for src,srcgrp in zip(args.srcid, assocsrcs):
                cur.execute("""SELECT assocxtrsource.runcat, image.url, extractedsource.ra, extractedsource.decl, extractedsource.f_int, 
                    extractedsource.f_int_err, image.rms_qc FROM extractedsource JOIN assocxtrsource ON 
                    assocxtrsource.xtrsrc=extractedsource.id JOIN image ON image.id=extractedsource.image JOIN skyregion ON 
                    skyregion.id=image.skyrgn JOIN runningcatalog ON runningcatalog.id=assocxtrsource.runcat WHERE assocxtrsource.runcat 
                    IN %s;""",(tuple(srcgrp),))
                imposfetch = cur.fetchall()
                imposarr = np.array(imposfetch, dtype=[('id','i8'),('image','U128'),('ra','f8'),('dec','f8'),('fint','f8'),('finterr','f8'),('rmsqc','f8')])
                for subsrc in srcgrp:
                    if subsrc != src:
                       if np.array_equal(imposarr[['image','ra','dec']][imposarr['id']==subsrc] , imposarr[['image','ra','dec']][imposarr['id']==src]):
                           srcgrp.remove(subsrc)
                if len(srcgrp) > 1:
                    mainsc = SkyCoord(np.average(imposarr['ra'][imposarr['id']==src]),np.average(imposarr['dec'][imposarr['id']==src]), unit='deg')
                    filedates = []
                    plotfilenames = []
                    imnames = np.unique(imposarr['image'])
                    print(srcgrp)
                    for myim in imnames:
                        imagesrcs = imposarr[imposarr['image']==myim]
                        fig = plt.figure(figsize=(25,15))
                        try:
                            hdul = fits.open(myim)
                            hdu = hdul[0]
                            wcs = WCS(hdu.header, naxis=2)
                            wcslocation = wcs.world_to_pixel(mainsc)
                            mydate = (datetime.datetime.strptime(hdu.header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") - start_epoch).total_seconds()/3600/24 
                            filedates.append(mydate)
                            ax1 = fig.add_subplot(1,1,1,projection=wcs)
                            myvmax = imagesrcs['rmsqc'][0]*4
                            ax1.imshow(hdu.data[0:][0:][0][0], vmin=args.vmin, vmax=myvmax, origin='lower')
                            ax1.set_ylim(wcslocation[1] - 50, wcslocation[1] + 50 )
                            ax1.set_xlim(wcslocation[0] - 50 , wcslocation[0] + 50)
                            wcscurloc = wcs.world_to_pixel(mainsc)
                            ax1.scatter(wcscurloc[0], wcscurloc[1], color='red', marker='s')
                            for curimagesrc in imagesrcs:
                                wcstmploc = wcs.world_to_pixel(SkyCoord(curimagesrc['ra'],curimagesrc['dec'],unit='deg'))
                                ax1.scatter(wcstmploc[0], wcstmploc[1], color='red', marker='x')
                            ax1.set_title('Image '+myim.split('/')[-1]+', Date: '+str(mydate)+', Fint: ' + str(curimagesrc['fint'])+' +/- '+str(curimagesrc['finterr']))
                            name = 'mytmpimage'+myim.rsplit('/')[-1].replace('.fits','')+'associate.png'
                            print('Plotted ',name)
                            plt.savefig(name)
                            plt.close()
                            plotfilenames.append(name)
                            hdul.close()
                        except Exception as e:
                            print(e)
                            print('skipping file')
                    plotfilenamearr = np.array(plotfilenames)
                    datesort = np.argsort(np.array(filedates))
                    with imageio.get_writer('associate_'+str(src)+'.mp4', mode='I', fps=2) as writer:
                        for myfile in plotfilenamearr[datesort] :
                            image = imageio.imread(myfile)
                            writer.append_data(image)
                    print("wrote video: ",'associate_'+str(src)+'.mp4')
                    for myfile in np.unique(plotfilenames):
                        print('removing: ',myfile)
                        os.remove(myfile)
                for src in args.srcid:
                    try:
                        srcidlist.remove(src)
                    except Exception as e:
                        print(e)
        else:
            srcidlist = args.srcid

        srcidTOexamine = tuple(srcidlist,)
        datasetTOexamine = queryfordatasets("""runningcatalog.id IN %s;""",(srcidTOexamine,))
        cur.execute(""" SELECT DISTINCT skyregion.centre_ra, skyregion.centre_decl FROM extractedsource JOIN 
            assocxtrsource ON assocxtrsource.xtrsrc=extractedsource.id JOIN image ON image.id=extractedsource.image
            JOIN skyregion ON skyregion.id=image.skyrgn WHERE assocxtrsource.runcat IN
            %s;""",(srcidTOexamine,)) 
    elif bool(args.srcid) and bool(args.dataset):
        sys.exit("Specifying a src id necessarily restricts the dataset id. I\'m not sure what you want by specifying both at the same time. Pick either one or the other")
    else:
        cur.execute(""" SELECT skyregion.centre_ra, skyregion.centre_decl  FROM skyregion;"""
                        )
    print("Gathering pointings")
    pointing = cur.fetchall()
    pointing = np.array(pointing)
    pointingfloat = np.array(pointing[:,0:2], dtype=float)
    # I f'ed up and didn't keep my leap second database up to date, therefore there are small offsets 
    # in the reported pointing. We'll do a little rounding to work around this problem 
    uniquepointingfloat = np.unique(pointingfloat.round(decimals=2),axis=0)
    # Iterate over the pointings getting the image names from the pointing in question
    print(uniquepointingfloat)
    for point in uniquepointingfloat:
        print('Iterating over pointings')
        if bool(rmsqcmax) and bool(rmsqcmin):
            refinesqlstring = """ (image.rms_qc < %s) AND (image.rms_qc > %s)"""
            inserttuple = (rmsqcmax, rmsqcmin,)
        elif bool(rmsqcmax) and not bool(rmsqcmin):
            refinesqlstring= """ (image.rms_qc < %s)"""
            inserttuple = (rmsqcmax,)
        elif not bool(rmsqcmax) and bool(rmsqcmin):
            refinesqlstring = """  (image.rms_qc > %s)"""
            inserttuple = (rmsqcmin,)
        appendsqlstring = """;"""

        if bool(args.exclusionsfile):
            refinesqlstring = """ (runningcatalog.id NOT IN %s) AND """ + refinesqlstring
            insertlist = list(inserttuple)
            insertlist.insert(0, tuple(excludedid))
            inserttuple = tuple(insertlist)
        # if datasets not specified, get them now to include in sql query for varmetric array
        if not bool(args.dataset) and not (args.srcid):
            print("Getting datasets")
            datasetTOexamine = queryfordatasets(refinesqlstring+appendsqlstring,inserttuple)
            datasetlist = list(datasetTOexamine)
            # iterate over the datasets, some are identical, we assume that if the images inputted and config are identical, they are duplicates
            for dat1index in tqdm(range(len(datasetTOexamine))):
                for dataset2 in datasetlist:
                    if (dataset2 != datasetTOexamine[dat1index]) and (datasetTOexamine[dat1index] in datasetlist):
                        cur.execute("""SELECT config.section, config.key, config.value,  config.type, image.url FROM image JOIN config ON config.dataset=image.dataset WHERE image.dataset=%s""",(datasetTOexamine[dat1index],))
                        dat1info = cur.fetchall()
                        cur.execute("""SELECT config.section, config.key, config.value,  config.type, image.url FROM image JOIN config ON config.dataset=image.dataset WHERE image.dataset=%s""",(dataset2,))
                        dat2info = cur.fetchall()
                        if dat1info==dat2info:
                            cur.execute("""SELECT COUNT(runningcatalog.id) FROM runningcatalog WHERE runningcatalog.dataset=%s""", (datasetTOexamine[dat1index],))
                            dat1runcatlength = cur.fetchall()
                            cur.execute("""SELECT COUNT(runningcatalog.id) FROM runningcatalog WHERE runningcatalog.dataset=%s""", (dataset2,))
                            dat2runcatlength = cur.fetchall()
                            # Some runs were incomplete, so we throw away the dataset with the smaller runningcatalog
                            if dat1runcatlength[0][0] >= dat2runcatlength[0][0]:
                                # We update the list that the inner for loop is pulling from here, otherwise the loop will try to check from already deleted dataset numbers
                                datasetlist.remove(dataset2)
                            else:
                                datasetlist.remove(datasetTOexamine[dat1index])
                                # We break here to stop checking against dataset1 since it is now thrown away. It goes to the next item in the outer for loop
                                break
            if args.resume:
                existingfiles = glob.glob('dataset*scatter.png')
                existingnums = [int(f.split('dataset')[1].split('scatter')[0]) for f in existingfiles]
                print("before: ",datasetlist)
                beforedatasets = tuple(datasetlist)
                for dsnum in beforedatasets:
                    if dsnum in existingnums:
                        datasetlist.remove(dsnum)
                print('after: ',datasetlist)
            datasetTOexamine = tuple(datasetlist)
        if len(datasetTOexamine)==0:
            continue

        # set sql addendum for angular restriction of sources
        try:
            inserttuple
            varinsertlist = list(inserttuple)
        except NameError:
            varinsertlist = []
        varinsertlist.insert(0,str(args.angrestrict))
        varinsertlist.insert(0, datasetTOexamine)
        varinserttuple = tuple(varinsertlist)
        if bool(args.skysearch):
            furtherrefinedSQLstring = """ (image.dataset IN %s) AND (sqlseparation(extractedsource.ra, extractedsource.decl, %s, %s) < %s) ;""" 
            insertlist = list(inserttuple)
            insertlist.insert(0,str(searchlocation[2]/3600))
            insertlist.insert(0,str(searchlocation[1]))
            insertlist.insert(0,str(searchlocation[0]))
            insertlist.insert(0,datasetTOexamine)
            if bool(args.srcid):
                furtherrefinedSQLstring = """(runningcatalog.id IN %s) AND """+furtherrefinedSQLstring
                insertlist.insert(0,srcidTOexamine)
            BIGinserttuple = tuple(insertlist)
        else:
            furtherrefinedSQLstring = """ (image.dataset IN %s) AND (sqlseparation(extractedsource.ra, extractedsource.decl, skyregion.centre_ra, skyregion.centre_decl) < %s) """ 
            insertlist = list(varinserttuple)
            BIGinserttuple = tuple(insertlist) 
        # get datasets for pointing and get varmetric data for each dataset 
        varinsertlist = list(varinserttuple)
        if args.sigthresh is not None:
            varinsertlist.append(args.sigthresh)
            appendsqlstring = """ GROUP BY runningcatalog.id,varmetric.eta_int, varmetric.v_int, image.dataset HAVING MAX(extractedsource.det_sigma) > %s""" + appendsqlstring
        varinserttuple = tuple(varinsertlist)
        print('Acquiring varmetric data')
        try:
            refinesqlstring
            doand = """ AND""" + """refinesqlstring"""
        except NameError:
            doand = """ """
        varmetricdata = dovarmetricquery(""" (image.dataset in %s) AND (sqlseparation(extractedsource.ra, extractedsource.decl, skyregion.centre_ra, skyregion.centre_decl) < %s) """+doand+appendsqlstring, varinserttuple)
        print('Acquiring bright source data')
        try:
            stabarr = querystablesrcs(""" AND """+refinesqlstring+""";""", inserttuple)   
        except NameError:
            stabarr = querystablesrcs(""";""", 0)   
        print('Acquiring bigarr')
        try:
            bigarr = queryallsrcs(furtherrefinedSQLstring+""" AND """+refinesqlstring+""";""", BIGinserttuple)
        except NameError:
            bigarr = queryallsrcs(furtherrefinedSQLstring+""";""", BIGinserttuple)

        duplicates = []
        if args.deduplicate:
            print('Deduplicating')
            dummyi = 1
            for d in datasetTOexamine:
                if d in varmetricdata['dataset']:
                    for id in tqdm(np.unique(varmetricdata['id'][varmetricdata['dataset']==d])):
                        mysource = np.unique(bigarr[bigarr['id']==id])
                        if mysource.size > 0:
                            # Test query, make sure only ids are grabbed, add statement for f_peak option
                            cur.execute("""SELECT DISTINCT runningcatalog.id, runningcatalog_flux.avg_f_int FROM runningcatalog JOIN
                            image ON image.dataset=runningcatalog.dataset JOIN runningcatalog_flux ON runningcatalog.id=runningcatalog_flux.runcat 
                            WHERE image.dataset=%s AND SQLSEPARATION(runningcatalog.wm_ra, runningcatalog.wm_decl, %s, %s) < %s*image.rb_smaj 
                            ORDER BY runningcatalog_flux.avg_f_int DESC;""",
                            (d, np.around(np.average(mysource['ra']), decimals=4), np.around(np.average(mysource['dec']), decimals=4),args.beamwidths,))
                            fetchedsrcs = cur.fetchall()
                            if len(fetchedsrcs) > 1:
                                duplicates.extend([f for fet in fetchedsrcs[1:] for f in fet])                       
                        else:
                            print("no sources to deduplicate")
            duplicatesarr = np.array(duplicates[::2])
            beforelen = len(np.unique(bigarr['id']))
            bigarr = bigarr[(~np.isin(bigarr['id'], duplicatesarr))]
            print('Deduplicated ',beforelen-len(np.unique(bigarr['id'])),' sources')
            varmetricdata = varmetricdata[(~np.isin(varmetricdata['id'], duplicatesarr))]
        if bool(args.FDlimits):
            commondatasets = np.intersect1d(fdlimits['dataset'], bigarr['dataset'])
            excludedid = []
            for d in commondatasets:
                fdl = np.unique(fdlimits['fdlimit'][fdlimits['dataset'] == d])[0]
                if args.usefpk:
                    cur.execute("""SELECT runningcatalog.id FROM runningcatalog JOIN extractedsource ON extractedsource.ff_runcat=runningcatalog.id WHERE (runningcatalog.dataset = %s) GROUP BY runningcatalog.id HAVING MAX(extractedsource.f_peak) < %s;""", (d, fdl,))
                else:
                    cur.execute("""SELECT runningcatalog.id FROM runningcatalog JOIN extractedsource ON extractedsource.ff_runcat=runningcatalog.id  WHERE (runningcatalog.dataset = %s) GROUP BY runningcatalog.id HAVING MAX(extractedsource.f_int) < %s;""", (d, fdl,))
                fetchedids = cur.fetchall()
                excludedid.extend([f for fet in fetchedids for f in fet])
            excludedidarr = np.array(excludedid)
            bigarr = bigarr[~np.isin(bigarr['id'],excludedidarr)]
            varmetricdata = varmetricdata[(~np.isin(varmetricdata['id'],excludedidarr))]
        reletaVarr = getetavtoplot(varmetricdata)
        if bool(args.srcid):
            bigarr = bigarr[np.isin(bigarr['id'], srcidTOexamine)]
        print('dumpvarmetric',point)
        if args.dumpvarmetric:
            for d in datasetTOexamine:
                print('using dataset ',d)
                print(varmetricdata['dataset'])
                if d in varmetricdata['dataset']:
                    print('Now saving varmetric binary')
                    np.save(f'ds{d}pointra{point[0]}pointdec{point[1]}varmetric.npy',varmetricdata[(varmetricdata['dataset']==d) & (varmetricdata['eta'] > 0) & (varmetricdata['v'] > 0)])
  #                   np.savetxt(f'ds{d}_point_{point}_varmetric.csv', varmetricdata[(varmetricdata['dataset']==d) & (varmetricdata['eta'] > 0) & (varmetricdata['v'] > 0)],delimiter=',')
        if args.makeplots:
                # If this is zero, there is nothing to plot
            if bigarr.shape[0] != 0 : 
                # Get an array of dates, preserve order, make MJD
                dates = np.array([(b - start_epoch).total_seconds()/3600/24 for b in bigarr['date']])
                # Do the same, but with stabarr, since stabarr has a more complete dataset
                stabdates = np.array([(s - start_epoch).total_seconds()/3600/24 for s in stabarr['date']])
                if args.usefpk:
                    stabydat1 = stabarr['f_peak']
                else:
                    stabydat1 = stabarr['f_int']
                # Get a common x/y range for eta-v plot for all plots in an observation
                etamin = np.amin(varmetricdata['eta'][varmetricdata['eta']>0])
                etamax = np.amax(varmetricdata['eta'][varmetricdata['eta']>0])
                vmin = np.amin(varmetricdata['v'][varmetricdata['v']>0])
                vmax = np.amax(varmetricdata['v'][varmetricdata['v']>0])
                if args.parallel and len(datasetTOexamine) > 2:
                    with Pool() as pool:
                        pool.map(makeplots,datasetTOexamine)
                else:
                    for d in datasetTOexamine:
                        makeplots(d)
            else: 
                print("no sources")
        if args.makemovies:
            moviedir = 'point'+str(point).replace('[','').replace(']','').replace(' ','_').replace('.','p')
            if not os.path.exists(moviedir):
                os.makedirs(moviedir)

            if bigarr.shape[0] != 0 : 
                # Get a common x/y range for eta-v plot for all plots in an observation
                xmin = np.amin(varmetricdata['eta'][varmetricdata['eta']>0])
                xmax = np.amax(varmetricdata['eta'][varmetricdata['eta']>0])
                ymin = np.amin(varmetricdata['v'][varmetricdata['v']>0])
                ymax = np.amax(varmetricdata['v'][varmetricdata['v']>0])
                if args.paralleldataset and (len(datasetTOexamine) > 2):
                    with Pool() as pool:
                        pool.map(dsparallelanimatesrc, datasetTOexamine)
                else:
                    for d in datasetTOexamine:
                        fdl = np.unique(fdlimits['fdlimit'][fdlimits['dataset']==d])
                        v = varmetricdata[varmetricdata['dataset']==d]
                        alreadyplotted = []
                        sortedv = np.sort(v['v'])
                        sortedeta = np.sort(v['eta'])
                        intersectbigrel = np.intersect1d(bigarr['id'],reletaVarr['id'])
                        if args.animateinteresting:
                            srcstoconsider = getinterestingsrcs(d)
                            args.srcid = list(srcstoconsider)
                            reletaVarr = getetavtoplot(varmetricdata)
                        else:
                            srcstoconsider = np.intersect1d(intersectbigrel,v['id'])
                        if not args.skipsrcmovies:
                            print("making ", len(srcstoconsider)," movies")
                            if args.parallel and len(srcstoconsider) > 2:
                                with Pool() as pool:
                                    pool.map(animatesrc, srcstoconsider)
                            else:
                                for src in srcstoconsider:
                                    animatesrc(src)
                        if not args.noimage and args.bigpicture:
                            animatefield(srcstoconsider)
                        for f in glob.glob('src*mp4'):
                            os.rename(f,moviedir+'/'+f)
                        if args.animateinteresting:
                            args.srcid = []

                                 # Some integration images aren't finished running though TraP yet. These will show an
                                 # IndexError. We skip these and delete the axis. A simple way of dealing with a problem
                                 # that I don't want to worry about
            else: 
                print("no sources")
        if args.makeimages:
            print("making multi-images")
            
            if bigarr.shape[0] != 0 : 
                # Get a common x/y range for eta-v plot for all plots in an observation
                xmin = np.amin(varmetricdata['eta'][varmetricdata['eta']>0])
                xmax = np.amax(varmetricdata['eta'][varmetricdata['eta']>0])
                ymin = np.amin(varmetricdata['v'][varmetricdata['v']>0])
                ymax = np.amax(varmetricdata['v'][varmetricdata['v']>0])
                for d in datasetTOexamine:
                    fdl = np.unique(fdlimits['fdlimit'][fdlimits['dataset']==d])
                    v = varmetricdata[varmetricdata['dataset']==d]
                    alreadyplotted = []
                    sortedv = np.sort(v['v'])
                    sortedeta = np.sort(v['eta'])
                    intersectbigrel = np.intersect1d(bigarr['id'],reletaVarr['id'])
                    srcstoconsider = np.intersect1d(intersectbigrel,v['id'])
                    if not args.skipsrcmovies:
                        if args.parallel and len(srcstoconsider) > 2:
                            with Pool() as pool:
                                pool.map(makeimages, srcstoconsider)
                        else:
                            for src in srcstoconsider:
                                makeimages(src)

                             # Some integration images aren't finished running though TraP yet. These will show an
                             # IndexError. We skip these and delete the axis. A simple way of dealing with a problem
                             # that I don't want to worry about
            else: 
                print("no sources")

            
if args.associatesources: 
    for srcgrp,mainsrc in zip(assocsrcs,args.srcid):
        if not os.path.exists('src'+str(mainsrc)):
            os.makedirs('src'+str(mainsrc))
        for src in srcgrp:
            myfile = glob.glob('*src*' + str(src) + '*mp4')
            os.replace(myfile[0], 'src'+str(mainsrc)+'/'+myfile[0])
    

