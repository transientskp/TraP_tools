# TraP_tools
A python library of tools to process and analyse TraP data, including example scripts and Jupyter notebooks

## Folders

### databaseTools
Tools to access and query the TraP databases
* **dbtools.py** - contains:
  * access(engine,host,port,user,password,database) - opens a connection to the TraP database using SQLalchemy
  * GetVarParams(session,dataset_id) - returns the varmetric and runningcatalogue databases for a specific dataset id

### exampleScripts
Example Jupyter Notebooks and standalone scripts using these tools
* **FilterVariables.py** and **FilterVariables.ipynb** Are equivalent, one is a standalone script and the other is a Jupyter notebook. They plot the variability parameters for a specific dataset.
* **dblogin.py** contains the main login parameters required for the database

### PreTraPimageQC
Tools to prepare images for TraP processing and to conduct initial image quality control
* **script1.py** This is a temporary summary

### tools
Various other useful tools 
* **tools.py** Includes:
  * SigmaFit(data) - fits a gaussian distribution to data and returns the fit parameters
  * extract_data(filename) - extracts data from a csv file
  * write_data(filename,tmp) - writes data into a csv file