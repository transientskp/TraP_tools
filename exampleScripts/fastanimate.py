import datetime
import imageio
import matplotlib.pyplot as plt
import os
import numpy as np 
import argparse
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import Quadrangle
from astropy import units as u
import psycopg2
from psycopg2.extensions import register_adapter, AsIs
from multiprocessing import Pool

# psycopg2 complains about numpy datatypes, this avoids that error. Got this from stackoverflow but don't recall where
def addapt_numpy_float64(numpy_float64):
    return AsIs(numpy_float64)
def addapt_numpy_int64(numpy_int64):
    return AsIs(numpy_int64)
register_adapter(np.float64, addapt_numpy_float64)
register_adapter(np.int64, addapt_numpy_int64)

def makepreimages(imnum,frame):
    fig = plt.figure(figsize=(25,10))
    axlist = []
    with fits.open(frame['image']) as hdul:
        src = srcdat[0]
        hdu = hdul[0]
        imdat = hdu.data[0:][0:][0][0]
        wcs = WCS(hdu.header, naxis=2)
        sc = SkyCoord(src['ra'],src['dec'],unit='deg')
        wcslocation = wcs.world_to_pixel(sc)
        axlist.append(fig.add_subplot(1,3,1,projection= wcs))
        axlist[0].imshow(imdat, vmin=-2e-5, vmax=np.amax(srcdat['fint']), origin='lower')
        axlist[0].set_xlim(wcslocation[0] - 50 , wcslocation[0] + 50)
        axlist[0].set_ylim(wcslocation[1] - 50, wcslocation[1] + 50 )
        axlist.append(fig.add_subplot(1,3,2))
        curdates = srcdat['date'][np.abs(srcdat['date'] - src['date']) < 24]
        curfints = srcdat['fint'][np.abs(srcdat['date'] - src['date']) < 24]
        curfinterrs = srcdat['finterr'][np.abs(srcdat['date'] - src['date']) < 24]
        axlist[1].scatter(curdates,curfints)
        axlist[1].set_box_aspect(imdat.shape[0]/imdat.shape[1])
        axlist[1].errorbar(curdates,curfints,yerr=curfinterrs, fmt='none')
        axlist[1].set_ylim(1e-6,curfints.max())
        axlist[1].set_xlim(curdates.min(),curdates.max())
        axlist[1].set_xlabel('MJD')
        axlist[1].set_ylabel('$F_{int}$')
        axlist[1].axvline(src['date'])
        axlist.append(fig.add_subplot(1,3,3,projection= wcs))
        axlist[2].imshow(imdat, vmin=-2e-5, vmax=np.amax(srcdat['fint']), origin='lower')
        axlist[2].add_patch(Quadrangle((sc.ra, sc.dec), 100*u.arcsec, 100*u.arcsec,
                    edgecolor='white', facecolor='none',transform=axlist[2].get_transform('fk5')))

        plt.savefig(f'tmpim{imnum}.png')
        imnum += 1 
    plt.close()
def plotsrcs(imnum,src,srcdat):
    fig = plt.figure(figsize=(25,10))
    axlist = []
    with fits.open(src['image']) as hdul:
        hdu = hdul[0]
        imdat = hdu.data[0:][0:][0][0]
        wcs = WCS(hdu.header, naxis=2)
        sc = SkyCoord(src['ra'],src['dec'],unit='deg')
        sc2 = SkyCoord(srcdat['ra'][0],srcdat['dec'][0],unit='deg')
        wcslocation = wcs.world_to_pixel(sc)
        wcslocation2 = wcs.world_to_pixel(sc2)
        axlist.append(fig.add_subplot(1,3,1,projection= wcs))
        axlist[0].imshow(imdat, vmin=-2e-5, vmax=np.amax(srcdat['fint']), origin='lower')
        axlist[0].scatter(wcslocation[0],wcslocation[1],marker='x',color='red')
        axlist[0].set_xlim(wcslocation2[0] - 50 , wcslocation2[0] + 50)
        axlist[0].set_ylim(wcslocation2[1] - 50, wcslocation2[1] + 50 )
        axlist.append(fig.add_subplot(1,3,2))
        curdates = srcdat['date'][np.abs(srcdat['date'] - src['date']) < 24]
        curfints = srcdat['fint'][np.abs(srcdat['date'] - src['date']) < 24]
        curfinterrs = srcdat['finterr'][np.abs(srcdat['date'] - src['date']) < 24]
        axlist[1].scatter(curdates,curfints)
        axlist[1].errorbar(curdates,curfints,yerr=curfinterrs, fmt='none')
        axlist[1].set_ylim(1e-6,curfints.max())
        axlist[1].set_xlim(curdates.min(),curdates.max())
        axlist[1].set_xlabel('MJD')
        axlist[1].set_ylabel('$F_{int}$')
        axlist[1].set_box_aspect(imdat.shape[0]/imdat.shape[1])
        axlist[1].axvline(src['date'])
        axlist.append(fig.add_subplot(1,3,3,projection= wcs))
        axlist[2].imshow(imdat, vmin=-2e-5, vmax=np.amax(srcdat['fint']), origin='lower')
        axlist[2].add_patch(Quadrangle((sc.ra, sc.dec), 100*u.arcsec, 100*u.arcsec,
                    edgecolor='white', facecolor='none',transform=axlist[2].get_transform('fk5')))
        plt.savefig(f'tmpim{imnum}.png')
    plt.close()
# In order to use MJD
start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)

parser = argparse.ArgumentParser()
parser.add_argument('srcid',nargs='*',help='sources to animate')
args = parser.parse_args()
if __name__ == '__main__':

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
    for myid in args.srcid:
    
    
        query = """SELECT extractedsource.f_int, extractedsource.f_int_err,extractedsource.ra, extractedsource.decl, image.url, image.taustart_ts FROM extractedsource JOIN assocxtrsource ON extractedsource.id=assocxtrsource.xtrsrc JOIN image ON extractedsource.image=image.id WHERE assocxtrsource.runcat=%s;"""
        cur.execute(query,(myid,))
    
    
        convdate = [(fint,finterr,ra,dec,image,(d-start_epoch).total_seconds()/3600/24) for fint,finterr,ra,dec,image,d in cur.fetchall()]
        srcdat =np.sort(np.array(convdate, dtype=[('fint','f8'),('finterr','f8'),('ra','f8'),('dec','f8'),('image','<U256'),('date','f8')]),order='date')
    
        query = """SELECT image.taustart_ts, image.url FROM image JOIN runningcatalog ON image.dataset=runningcatalog.dataset WHERE runningcatalog.id=%s;"""
    
        cur.execute(query,(myid,))
    
        convalldate = [((d-start_epoch).total_seconds()/3600/24,image) for d,image in cur.fetchall()]
    
        alldat = np.sort(np.array(convalldate, dtype=[('date','f8'),('image','<U256')]),order='date')
    
        startindex = max(int(np.unique((alldat['date']==srcdat['date'].min()).nonzero()).astype(int) - 10),0)
        endindex = min(int(np.unique((alldat['date']==srcdat['date'].min()).nonzero()).astype(int)),startindex+10)
        print(myid)
        if endindex > startindex:
            with Pool() as pool:
                pool.starmap(makepreimages, enumerate(alldat[startindex:endindex]))


        preimnum = endindex - startindex
        endimnum = len(srcdat) + preimnum
        with Pool() as pool:
            pool.starmap(plotsrcs,[(preimnum+ind,src,srcdat) for ind,src in enumerate(srcdat)])
        imnum = endimnum 
        with imageio.get_writer(f'src{myid}.mp4', mode='I',fps=3) as writer:
            for ind in range(imnum):
                image = imageio.imread(f'tmpim{ind}.png')
                writer.append_data(image)
        print(f"wrote video: src{myid}.mp4")
        for myfile in [f'tmpim{ind}.png' for ind in range(imnum)]:
            print('removing: ',myfile)
            os.remove(myfile)









