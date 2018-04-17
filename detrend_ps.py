from __future__ import division,print_function
import numpy as np
import scipy
#from scipy import stats
#from scipy import interpolate
import scipy.optimize
#import sys
import os
from os.path import expanduser
import pyfits
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import fnmatch
import gatspy.periodic as gp
import scipy.ndimage
import re
#from emcee import autocorr
#import batman
#from pybls import BLS

def rebin(x, r):
    m = len(x) // r
    return x[:m*r].reshape((m,r)).mean(axis=1),r


def polysmooth(time, flux, degree):
   """Takes a time series and fits a polynomial
      of input degree, then divides by polynomial.
   """

   pfit = np.polyfit(time, flux, degree)
   fit = np.poly1d(pfit)
   xfit = np.linspace(min(time), max(time), len(time))

   flux2 = flux / fit(xfit)

   return flux2#, fit, xfit


def med_filt(x, y, dt):
    """
    De-trend a light curve using a windowed median.
    """
    x, y = np.atleast_1d(x), np.atleast_1d(y)
    assert len(x) == len(y)
    r = np.empty(len(y))
    for i, t in enumerate(x):
        inds = (x >= t - 0.5 * dt) * (x <= t + 0.5 * dt)
        r[i] = np.nanmedian(y[inds])
    return r

def ps(time,flux,N=1):
	time=time-time[0] #start time at zero
	if time[1]<1:
		time=time*86400 #change time to seconds (rather than days)

	c=[]
	for i in range(len(time)-1):
		c.append(time[i+1]-time[i])
	c=np.median(c)
	#print(c)
	nyq=1/(2*(time[1]-time[0]))
	nyq=1/(2*c) #calculate nyquist frequency
	#print(nyq*1e6)
	df=1/time[-1] #spacing between frequencies for lomb scargle

	f,p=gp.lomb_scargle_fast.lomb_scargle_fast(time,flux,f0=0,df=df,Nf=N*(nyq/df))
	lhs=(1/len(time))*np.sum(flux**2) 
	rhs= np.sum(p)
	ratio=lhs/rhs
	p*=ratio/(df*1e6)# change power to ppm^2/uHz
	f*=1e6 #change frequency to microhertz
	return f,p

def sigclip(time,flux,eflux=None,sigma=3,repeat=False):
	idx=np.ones(len(time))
	if repeat==True:
		loop=0
		while len(idx)>0:
			idx=np.where(np.abs(flux)>sigma*np.std(flux))[0]
			flux=np.delete(flux,idx,axis=0)
			if eflux != None:
				eflux=np.delete(eflux,idx,axis=0)
			time=np.delete(time,idx,axis=0)
			loop+=1
			print(loop)
	else:
		idx=np.where(np.abs(flux)>sigma*np.std(flux))[0]
		flux=np.delete(flux,idx,axis=0)
		if eflux != None:
			eflux=np.delete(eflux,idx,axis=0)
		time=np.delete(time,idx,axis=0)
	if eflux!= None:
		return time, flux, eflux
	else:
		return time,flux


def detrend(mypath,start=0,end=None,median=True,dt=10.,clip=True,fits=False,ktwo=False,sc=False):
	lc=[]
	if fits==True:
		for file in os.listdir(mypath):
			if ktwo==True:
				if fnmatch.fnmatch(file, '*ktwo**c.fits'):
					lc.append(file)
			else:
				if fnmatch.fnmatch(file, '*kplr**c.fits'):
					lc.append(file)
	else:
		for file in os.listdir(mypath):
			if ktwo==True:
				if fnmatch.fnmatch(file, '*ktwo**c.dat'):
					lc.append(file)
			else:
				if fnmatch.fnmatch(file, '*kplr**c.dat'):
					lc.append(file)
				if fnmatch.fnmatch(file, '*KIC**.txt'):
					lc.append(file)
	#print (lc)

	choice=[]
	for file in lc:
		if sc==True:
			if fnmatch.fnmatch(file, '*slc*'):
				choice.append(file)
		else:
			if fnmatch.fnmatch(file, '*llc*'):
				choice.append(file)
			if fnmatch.fnmatch(file, '*.txt'):
				choice.append(file)
	lc=choice
	lc=np.sort(lc)
	lc=lc[start:end]
	print(lc)

	for i in lc:
		if fits==True:
			datafile=pyfits.open(mypath+'/'+str(i))
			lightdata=datafile[1].data
			t=lightdata.field("TIME")            #Barycenter corrected Julian date #Start with the PDCSAP_FLUX if  you're not sure which one to use.
			f=lightdata.field("PDCSAP_FLUX")
			ef=lightdata.field("PDCSAP_FLUX_ERR")
		else:
			t,f,ef=np.loadtxt(mypath+'/'+str(i),unpack=True)#,usecols=(0,3,4))

        # changes any infinite values to nans
		f[np.isfinite(f)==False]=np.nan
		ef=ef/np.nanmedian(f) #divides by the median not including nans
		f=(f/np.nanmedian(f))

		finite=np.isfinite(f)
		f=f[finite] #take only finite values (since nans aren't finite)
		ef=ef[finite]
		t=t[finite]
		if median==True: #moving median filter	
			med=med_filt(t,f,dt)	
			plt.plot(t,f,'.')	
			plt.plot(t,med)	
			plt.show()	
			f=1e6*((f/med)-1)
		else:
			f=1e6*(f-1)
		ef=1e6*ef
		##############################################################
        
		if clip==True:
			print('Clip on')
			t,f,ef=sigclip(t,f,ef,repeat=True,sigma=3)   
			#f = f/np.max(np.abs(f))    # attempt to make the light curve all follow along one level        

		if i == lc[0]:	
			data=np.c_[t,f,ef]
		else:
			data1=np.c_[t,f,ef]
			data=np.append(data,data1,axis=0) #append new data set to combine to one bit time series

	idx=np.argsort(data[:,0]) #order by time
	data=data[idx]
	time,flux,eflux=data.T
	'''if clip==True:
		print('Clip on')
		time,flux,eflux=sigclip(time,flux,eflux,repeat=True,sigma=3)'''

	return time,flux,eflux

home=expanduser('~')
#data location
mypath='../KIC8006161/short_cadence'
start = 13
end = 16
time,flux,eflux=detrend(mypath,dt=3.,clip=True,median=True,fits=True,ktwo=False,sc=True,start=start,end=end)

f,p=ps(time,flux,N=1)
bw=f[1]-f[0]
bins=int(5/bw)
f1=rebin(f,bins)[0]
p1=rebin(p,bins)[0]

#kasoc=pyfits.open('kplr009145861_kasoc-psd_llc_v1.fits')
#fk,pk=kasoc[1].data['FREQUENCY'],kasoc[1].data['PSD']
#fkr,pkr=rebin(fk,bins)[0],rebin(pk,bins)[0]
time *=86400 #time now in seconds
plt.figure()
plt.plot(time,flux,'k.')
plt.show()

plt.figure()
plt.plot(f,p,'k',alpha=0.9)
#plt.plot(f1,p1,'r')
#plt.plot(fkr,pkr,'b')
plt.yscale('log')
plt.show()

# get rid of very small power values
'''p[np.abs(p) < 10E-1] = 0.75
plt.figure()
plt.plot(f,p,'k',alpha=0.9)
plt.yscale('log')
plt.show()'''

'''def get_file_number(filename):#get numbers for naming file
    regex = re.compile(r'\d+') 
    return regex.findall(filename)
    
allfiles = os.listdir(mypath)
start_file = get_file_number(allfiles[start])[1] + '.' + get_file_number(allfiles[start])[2]
end_file = get_file_number(allfiles[end-1])[1] + '.' + get_file_number(allfiles[end-1])[2]

filename = os.path.join('../eta_parameter','ps_' + start_file + '_' + end_file + '.csv')
print (filename)
#np.savetxt(filename, np.c_[f,p])'''
#np.savetxt('../KIC8006161/individual_quarters/KIC8006161_Q17.txt',np.column_stack((f,p)))
    