# INSTRUCTIONS: 
# 1. import all the environment variables by executing "source trapenvvar.bash" or whatever your file is called
# 2. get the dataset number from banana (should read dataset #X) where X is the number you want
# 3. execute by running python fetchallsources.py X
# 4. output file will be named datasetXdump.csv
# NOTE, file output may be very large in some cases and could take some time. 

import numpy as np
import psycopg2
import datetime
import os 
import sys
import glob
import argparse    

host = os.getenv('TKP_DBHOST')
port = os.getenv('TKP_DBPORT')
user = os.getenv('TKP_DBUSER')
password = os.getenv('TKP_DBPASSWORD')
database = os.getenv('TKP_DBNAME')
if password is None:
    conn = psycopg2.connect("dbname="+database+" user="+user+" host="+host)
else:
    conn = psycopg2.connect("dbname="+database+" user="+user+" password="+password+" host="+host)
cur = conn.cursor()


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('dataset', type=int, help='Dataset to retrieve')
args = parser.parse_args()
# In order to use MJD
start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)

cur.execute("""SELECT DISTINCT assocxtrsource.runcat, image.taustart_ts, extractedsource.f_int, 
    extractedsource.f_int_err, extractedsource.f_peak, extractedsource.f_peak_err FROM extractedsource JOIN 
    image ON image.id=extractedsource.image JOIN assocxtrsource ON assocxtrsource.xtrsrc=extractedsource.id 
    WHERE image.dataset=%s ;""",
    (args.dataset,))
fetched = cur.fetchall()
with open(f'dataset{args.dataset}dump.csv', 'w') as f:
    f.write("id,date,F_int,F_int_err,F_pk,F_pk_err\n")
    for dat in fetched:
        f.write(f'{dat[0]},{(dat[1] - start_epoch).total_seconds()/3600./24.},{dat[2]},{dat[3]},{dat[4]},{dat[5]}\n')

