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
    # Access the database using sqlalchemy
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

def GetPandaExtracted(session,dataset_id,**kwargs):
    x = session.query(Runningcatalog).filter(Runningcatalog.dataset_id == dataset_id)
    dx = pd.read_sql_query(x.statement,db.connection)
    # here we rename entries in runcat to their proper names
    dx = dx.rename(index=str,columns={'id' : 'runcat','xtrsrc':'id'})
    # here we drop duplicate columns from runcat
    dx = dx.drop(columns = {'x','y','z','zone'})
    print dx.keys()

    if kwargs is not None:
        for key,value in kwargs.iteritems():
            if value.lower() == "extractedsource":
                y = session.query(Extractedsource).join(Runningcatalog.xtrsrc).filter(Runningcatalog.dataset_id == dataset_id)
            if value.lower() == "varmetric":
                y = session.query(Extractedsource).join(Runningcatalog).filter(Runningcatalog.dataset_id == dataset_id)
            if value.lower() == "extracteource":
                y = session.query(Extractedsource).join(Runningcatalog.xtrsrc).filter(Runningcatalog.dataset_id == dataset_id)
            if value.lower() == "extradsource":
                y = session.query(Extractedsource).join(Runningcatalog.xtrsrc).filter(Runningcatalog.dataset_id == dataset_id)

            dy = pd.read_sql_query(y.statement,db.connection)
            try:
                dx = pd.merge(dx,dy,on=['id','runcat'])
            except:
                print "error no column runcat&id trying merge on id:"
                try:
                    dx = pd.merge(dx,dy,on=['id'])
                except: "error no column id trying merge on runcat:"
                    dx = pd.merge(dx,dy,on=['runcat'])
            dy = []

    # dx = pd.read_sql_query(x.statement,db.connection)
    # dy = pd.read_sql_query(y.statement,db.connection)

    # # here we rename entries in runcat to their proper names
    # dx = dx.rename(index=str,columns={'id' : 'runcat','xtrsrc':'id'})
    # # here we drop duplicate columns from runcat
    # dx = dx.drop(columns = {'x','y','z','zone'})
    # here we merge runcat and extractedsource
    # merged = pd.merge(dx,dy,on=['id'])
    return merged
