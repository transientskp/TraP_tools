# runcats_to_ds9_region.py
# 
# A code to calculate the position offset of each extracted source
# relative to its average runcat position.
# A histogram is produced using all of the offsets and then fit with
# a Gaussian distribution.

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
from tkp.db.model import Runningcatalog

# The input database, dataset and thresholds
dataset_id = 31
database = 'AR_testing6'

# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)
runcats= session.query(Runningcatalog).filter(Runningcatalog.dataset_id == dataset_id).order_by(Runningcatalog.id).all()

with open('ds'+str(dataset_id)+'_runcat_regions.reg', 'w') as f:
    f.write('# Region file format: DS9 version 4.1 \n')
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
    f.write('fk5 \n')
    for src in runcats:
        f.write('circle('+str(src.wm_ra)+','+str(src.wm_decl)+',30") \n')

