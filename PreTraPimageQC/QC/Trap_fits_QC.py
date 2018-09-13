import sys
import os
import pyfits
import numpy as np
import glob
import getRMS
import tools

hdu=0
if len(sys.argv) != 5:
    print 'python TraP_fits_QC.py <sigmaClip> <sigmaRej> <f> <frequencies>'
    exit()
#Usually 3 is setting from TraP (is rms_est_sigma)
sigmaClip = float(sys.argv[1])
#The value you are going to clip the RMS values from (how accurate you gonna be)(here big is bad(assuming most images are good))
sigmaRej = float(sys.argv[2])
#fraction image you use to calculate rms (is in trap file).. (is rms_est_fraction) IS 8 USUALLY
f = float(sys.argv[3])
# F since we are averaging over all frequencies (F=False)
allFreqs=sys.argv[4]

freqs=[]
imageData=[]
for filename in glob.glob('*image-pb.fits'):
    print "read: "  + str(filename)
    data = getRMS.read_data(pyfits.open(filename)[hdu],filename, plane=None)
    hdulist=pyfits.open(filename)
    prihdr=hdulist[0].header
    freq=int((float(prihdr['CRVAL3'])/1e6))
    imageData.append([filename, freq, getRMS.rms_with_clipped_subregion(data, sigmaClip, f)*1000.])
    if freq not in freqs:
        freqs.append(freq)

thresholds={}

if allFreqs=='T':
    for frequency in freqs:
        noise_avg_log, noise_scatter_log, noise_threshold_log = tools.fit_hist([np.log10(n[2]) for n in imageData if n[1]==frequency], sigmaRej, r'Observed RMS (mJy)', 'rms', frequency)
        noise_avg=10.**(noise_avg_log)
        noise_max=10.**(noise_avg_log+noise_scatter_log)-10.**(noise_avg_log)
        noise_min=10.**(noise_avg_log)-10.**(noise_avg_log-noise_scatter_log)
        print 'Average RMS Noise in images (1 sigma range, frequency='+str(frequency)+' MHz): '+str(frequency)+' MHz): '+str(noise_avg)+' (+'+str(noise_max)+',-'+str(noise_min)+') mJy'
        thresholds.append([frequency,noise_max*sigma,noise_min*sigma])


frequency='all'
TMPdata = np.array([np.log10(n[2]) for n in imageData])
TMPdata = TMPdata[np.isfinite(TMPdata)]
noise_avg_log, noise_scatter_log, noise_threshold_log = tools.fit_hist(TMPdata, sigmaRej, r'Observed RMS (mJy)', 'rms', frequency)
print noise_avg_log
noise_avg=10.**(noise_avg_log)
noise_max=10.**(noise_avg_log+noise_scatter_log)-10.**(noise_avg_log)
noise_min=10.**(noise_avg_log)-10.**(noise_avg_log-noise_scatter_log)
print 'Average RMS Noise in images (1 sigma range, frequency='+str(frequency)+' MHz): '+str(noise_avg)+' (+'+str(noise_max)+',-'+str(noise_min)+') mJy'
thresholds[frequency]=[noise_avg,noise_max*sigmaRej,noise_min*sigmaRej*-1.]
print thresholds
goodImg=[]
for image in imageData:
    if image[2]<thresholds['all'][0]+thresholds['all'][1] and image[2]>thresholds['all'][0]+thresholds['all'][2]:
        goodImg.append(image[0])
    else:
        print 'Bad image:',image[0],image[1],image[2]

trapmin = image[2]>thresholds['all'][0]+thresholds['all'][2]
trapmax = image[2]<thresholds['all'][0]+thresholds['all'][1]
print 'Instead of using the images_to_process file,'
print 'Minimum Value for TraP='+str(trapmin)
print 'Maximum Value for TraP='+str(trapmax)

f = open('images_to_process.py', 'w')
f.write('images = ')
f.write('[')
f.writelines(["'%s',\n" % good for good in goodImg])
f.write(']\n')
f.write('''#Just for show:
print "******** IMAGES: ********"
for f in images:
    print f
print "*************************"
''')
f.close()
