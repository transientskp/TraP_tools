import tkp.db
from tkp.db.model import Varmetric
from tkp.db.model import Runningcatalog
from tkp.db.model import Newsource
from tkp.db.model import Extractedsource
from tkp.db.model import Image
from sqlalchemy import *
from sqlalchemy.orm import relationship
import pandas as pd



def access(engine,host,port,user,password,database):
    """ Access the database using sqlalchemy"""
    # make db global in order to be used in GetPandaExtracted
    global db
    db = tkp.db.Database(engine=engine, host=host, port=port,
                     user=user, password=password, database=database)
    db.connect()
    session = db.Session()
    print 'connected!'
    return session

def GetVarParams(session,dataset_id):
    # Returns all the variability parameters for sources in a given dataset
    VarParams = session.query(Varmetric,Runningcatalog).select_from(join(Varmetric,Runningcatalog)).filter(Runningcatalog.dataset_id == dataset_id).all()
    return VarParams

def MergeTables(session,dataset_id,value,dx):
    """ Merges Tables for pandas output"""
    # if table name = extractedsource
    if value.lower() == "extractedsource":
        y = session.query(Extractedsource).join(Runningcatalog.xtrsrc).filter(Runningcatalog.dataset_id == dataset_id)
        dy = pd.read_sql_query(y.statement,db.connection)

    # if table name = varmetric also make sure to change id name to varmetric
    if value.lower() == "varmetric":
        y = session.query(Varmetric).join(Runningcatalog).filter(Runningcatalog.dataset_id == dataset_id)
        dy = pd.read_sql_query(y.statement,db.connection)
        dy = dy.rename(index=str,columns={'id':'varmetric'})

    # if table name = newsource also make sure to rename columns id and trigger_xtrsrc
    if value.lower() == "newsource":
        y = session.query(Newsource).join(Runningcatalog).filter(Runningcatalog.dataset_id == dataset_id)
        dy = pd.read_sql_query(y.statement,db.connection)
        dy = dy.rename(index=str,columns={'id':'newsource','trigger_xtrsrc':'xtrsrc'})

    # if tablename = image make sure to rename id
    if value.lower() == "image":
        y = session.query(Image).filter(Image.dataset_id == dataset_id)
        dy = pd.read_sql_query(y.statement,db.connection)
        dy = dy.rename(index=str,columns={'id':'image'})

    try:
        # check for same columns and then merge on those
        same_col = dy.columns.intersection(dx.columns)
        dx = pd.merge(dx,dy,on=list(same_col))
    except:
        "please enter a valid table or add the table you want to dbtools.py"

    return dx

def GetPandaExtracted(session,dataset_id,tablelist):
    """Function for extracting data in pandas format"""

    # First we always get Runcat in order to be able to filter on dataset_id more easily
    x = session.query(Runningcatalog).filter(Runningcatalog.dataset_id == dataset_id)
    dx = pd.read_sql_query(x.statement,db.connection)
    dx = dx.rename(index=str,columns={'id' : 'runcat'})

    if isinstance(tablelist,list):
        # loop through all the tables provided
        for value in tablelist:
            dx = MergeTables(session,dataset_id,value,dx)
    else:
        dx = MergeTables(session,dataset_id,tablelist,dx)


    return dx
