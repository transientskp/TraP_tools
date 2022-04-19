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
* **4FreqVariablesPlot.py** uses the output of FilterVariables.py run
  separately on each observing frequency in the dataset for 4
  different observing frequencies. Outputs a 2x2 plot of the
  variability parameters at each observing frequency.
* **CatalogueMatching.py** is a jupyter notebook that can associate
two catalogues, one from TraP, with  each other.
* **correctSystematicFluxOffset.py** Sometimes there can be a
  systematic flux offest between images of the same field. This script
  uses the average flux density of point sources, as calculated using
  TraP, to correct the flux  density scale in all the images.
* **FilterVariables.py** and **FilterVariables.ipynb** Are equivalent,
  one is a standalone script and the other is a Jupyter notebook. They
  plot the variability parameters for a specific
  dataset. FilterVariables.py is a more up to date version of the
  code - as needed for 4FreqVariablesPlot.py.
* **filter_new_sources.py** A filtering strategy using TraP outputs
  and deep catalogues of the field. For more details see
  [Rowlinson et al. 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv220316918R/abstract)
* **SensitivityPlot.ipynb** Code to make an image noise map using the
  TraP outputs
* **dblogin.py** contains the main login parameters required for the database

### PreTraPimageQC
Tools to prepare images for TraP processing and to conduct initial image quality control
* **script1.py** This is a temporary summary

### plotting
Various plotting tools used by the example scripts

### tools
Various other useful tools 
* **tools.py** Includes:
  * SigmaFit(data) - fits a gaussian distribution to data and returns the fit parameters
  * extract_data(filename) - extracts data from a csv file
  * write_data(filename,tmp) - writes data into a csv file
  * animation_zoom.ipynb - IPython notebook which creates a movie from a collection of images, allowing for the user to specify a zoom-in window on a target of interest. Based on the animation prototype created by Mark in his scratchpad repo.
