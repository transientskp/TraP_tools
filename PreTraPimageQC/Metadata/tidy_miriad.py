# Corrects the metadata from miriad images in the working folder so that TraP can read them and recognise them as LOFAR fits images
# Following this you will need to add the metadata to the images using the standard tkp-inject.py script (http://docs.transientskp.org/tkp/r1.1/tools/tkp-inject.html)

import pyfits
import glob
import os
import sys

images= glob.glob(os.path.expanduser("*.fits"))

if len(sys.argv) != 1:
    print 'tidy_miriad.py <8CharTelescope>'
    exit()

tscope = str(sys.argv[1])

if len(tscope) !=7:
    print 'Telescope name must be 8 characters'
    exit()

for img in images:
    hdulist=pyfits.open(img, mode='update')
    prihdr=hdulist[0].header
    del prihdr['BLANK']
    frq = hdulist[0].header['CRVAL3']
    prihdr.update('RESTFRQ',frq)
    prihdr.update('RESTFREQ',frq)
    prihdr.update('TELESCOP', tscope)
    prihdr.update('SPECSYS', 'LSRK    ')
    hdulist.flush()
