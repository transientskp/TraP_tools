import numpy as np 
import pandas as pd
import argparse
from multiprocessing import Pool
from astropy.coordinates import SkyCoord
parser = argparse.ArgumentParser(
                    prog = 'deduplicate',
                    description = 'Given a csv file containing columns with at least  wm_ra, wm_decl, datapoints, rb_smaj, runcat reduce source density by eliminating sources within 5 beamwidths, keeping the source with the highest number of datapoints')
parser.add_argument('--csvfile',type=str)
parser.add_argument('--pullsql',action='store_true', help='option to pull data from the database, assuming trapenvvar.bash is already sourced.')
parser.add_argument('--dataset',type=int,help='dataset to pull from sql')
parser.add_argument('--getoutput',action='store_true',help='Also dump csv file for input into slidingetawindow.py')
parser.add_argument('--minrms',type=float, default=1e-16, help='Exclude measurements below this rmsqc value')
parser.add_argument('--maxrms',type=float, default=1e16,help='Exclude measurements above this rmsqc value')
args = parser.parse_args()
def labelrgn(myid):
#     print(myid)
    myra = float(dat['wm_ra'][dat['runcat']==myid].unique())
    mydec = float(dat['wm_decl'][dat['runcat']==myid].unique())
    mypts = int(dat['datapoints'][dat['runcat']==myid].unique())
    mymaj = np.amax(dat['rb_smaj'][dat['runcat']==myid])
    mysc = SkyCoord(myra,mydec,unit='deg')
    otherid = dat['runcat'][dat['runcat']!=myid].unique()
    otherra = dat['wm_ra'][dat['runcat']!=myid].unique()
    otherdec = dat['wm_decl'][dat['runcat']!=myid].unique()
    otheridpts = np.unique(dat[['runcat', 'datapoints']].to_numpy(),axis=0)
    alldist = mysc.separation(SkyCoord(otherra,otherdec,unit='deg')).degree
    samergn = alldist < (5*mymaj)
    if np.sum(samergn) > 0:
        outarr = np.zeros(np.sum(samergn)+1,dtype=[('id','i8'),('rgnid','i8')])
  #       print(np.sum(samergn),samergn.shape,otherid.shape,)
   #      print(len(outarr))
        outarr['id'][0:-1] = otherid[samergn]
        outarr['id'][-1] = myid
        outarr['rgnid'] = np.amin(outarr['id'])
    else:
        outarr = np.zeros(1,dtype=[('id','i8'),('rgnid','i8'),('pts','i8')])
        outarr['id'] = myid
        outarr['rgnid']=myid
    return [(o['id'],o['rgnid']) for o in outarr]

def initdb():
    import psycopg2
    import os
    from psycopg2.extensions import register_adapter, AsIs
    # psycopg2 complains about numpy datatypes, this avoids that error. Got this from stackoverflow but don't recall where
    def addapt_numpy_float64(numpy_float64):
        return AsIs(numpy_float64)
    def addapt_numpy_int64(numpy_int64):
        return AsIs(numpy_int64)
    register_adapter(np.float64, addapt_numpy_float64)
    register_adapter(np.int64, addapt_numpy_int64)
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
    return conn.cursor()

if __name__=='__main__':
    if args.pullsql:
        cur = initdb()

        
        query = """SELECT image.rms_qc, image.rb_smaj,runningcatalog.wm_ra, runningcatalog.wm_decl,runningcatalog.datapoints, image.taustart_ts, image.band, image.freq_eff, extractedsource.f_peak, extractedsource.f_peak_err, extractedsource.f_int, extractedsource.f_int_err, assocxtrsource.runcat, extractedsource.det_sigma, extractedsource.ra, extractedsource.decl, sqlseparation(extractedsource.ra,extractedsource.decl, skyregion.centre_ra, skyregion.centre_decl), image.dataset FROM extractedsource JOIN assocxtrsource ON assocxtrsource.xtrsrc=extractedsource.id JOIN image ON extractedsource.image=image.id JOIN skyregion ON image.skyrgn=skyregion.id  JOIN runningcatalog ON assocxtrsource.runcat=runningcatalog.id WHERE image.dataset=%s AND image.rms_qc<%s AND image.rms_qc>%s GROUP BY image.rb_smaj,image.taustart_ts, image.band, image.freq_eff, extractedsource.f_peak, extractedsource.f_peak_err, extractedsource.f_int, extractedsource.f_int_err, assocxtrsource.runcat, extractedsource.det_sigma, extractedsource.ra, extractedsource.decl, skyregion.centre_ra, skyregion.centre_decl, image.dataset, runningcatalog.wm_ra, runningcatalog.wm_decl, runningcatalog.datapoints, image.rms_qc HAVING sqlseparation(extractedsource.ra, extractedsource.decl, skyregion.centre_ra, skyregion.centre_decl) < 0.8"""
        print('Querying the database')
        cur.execute(query,(args.dataset,args.maxrms, args.minrms))
        
        dat = pd.DataFrame(cur.fetchall(), columns=['rms_qc','rb_smaj','wm_ra','wm_decl','datapoints','taustart_ts',
                                                  'band','freq_eff','f_peak','f_peak_err','f_int','f_int_err',
                                                  'runcat','det_sigma','ra','decl','dist','dataset'])
        print('Put data into dataframe!')
        print(dat['rms_qc'].min(), dat['rms_qc'].max())
    elif args.csvfile:
        dat = pd.read_csv(args.csvfile)
    
    
    else:
        print("Please specify data source")
        import sys
        sys.exit(1)
        
    uniqueids = dat['runcat'].unique()
    with Pool() as pool:
        results = pool.map(labelrgn, uniqueids)

    flatresults = []
    for r in results:
        flatresults.extend(r)
    locarray = np.array(flatresults)
    keepids = []
    for uniqueloc in np.unique(locarray):
        idinlocation =  uniqueids[np.isin(uniqueids,locarray[locarray[:,1]==uniqueloc,0])]
        subdat = dat[np.isin(dat['runcat'],idinlocation)]
        keepdat = subdat[subdat['datapoints']==subdat['datapoints'].max()]
        for myid in keepdat['runcat'].unique():
            keepids.append(myid)
    print(keepids)
    with open('dumpdedupedids.txt','w') as f:
        f.write(str(tuple(keepids)))
        
    if args.getoutput:
        import os
        database = os.getenv('TKP_DBNAME')
        if not args.pullsql:
            cur = initdb()
        query = """SELECT image.taustart_ts, image.band, image.freq_eff, extractedsource.f_peak, extractedsource.f_peak_err, extractedsource.f_int, extractedsource.f_int_err, assocxtrsource.runcat, extractedsource.det_sigma, extractedsource.ra, extractedsource.decl, sqlseparation(extractedsource.ra,extractedsource.decl, skyregion.centre_ra, skyregion.centre_decl), image.dataset FROM extractedsource JOIN assocxtrsource ON assocxtrsource.xtrsrc=extractedsource.id JOIN image ON extractedsource.image=image.id JOIN skyregion ON image.skyrgn=skyregion.id  WHERE  assocxtrsource.runcat IN %s AND image.rms_qc>%s AND image.rms_qc<%s ORDER BY assocxtrsource.runcat ASC ;"""
        cur.execute(query, (tuple(keepids),args.minrms, args.maxrms,))
        outdf = pd.DataFrame(cur.fetchall(),columns=['taustart_ts','band','freq_eff','f_peak','f_peak_err',
                                                   'f_int','f_int_err','runcat','det_sigma','ra','decl','dist','dataset'])
        with open(f'{database}deduped.csv','w') as f:
            f.write(outdf.to_csv(index=False))
            
   


    
