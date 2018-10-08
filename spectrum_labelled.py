#!/usr/bin/env python3

"""A short script to make a power spectrum from mode data in AADG3
format.

Mostly copied from
https://github.com/warrickball/AADG3/blob/master/python/AADG3.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def lorentz(nu, nu0, w):
    return 1./(1.+4.*(nu-nu0)**2/w**2)

def gaus_env(nu,A,numax,width):
    return A*np.exp(-(nu-numax)**2/(2*width**2))

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

modes = np.loadtxt('./data/s4tess_00771_mode_table.dat', usecols=range(5))

f = np.linspace(1.0, 4000.0, 40001)
p = np.zeros_like(f) + 0.03
p = np.random.exponential(p)
peaks = []

vis = [1.0, 1.5, 0.5, 0.05]  # typical normalised visibilities
i=0
for row in modes:
    l, n, freq, width, power = row
    height = 2.*power/np.pi/width*vis[int(l)]

    power_add = height*lorentz(f, freq, width)
    peaks = np.append(peaks,np.max(power_add))
    p += power_add
    
#need to find peaks
f_smooth = modes[:25,2]
p_smooth = peaks[:25]
p0=(np.max(p_smooth), f_smooth[np.argmax(p_smooth)], 300.)
popt, pcov = curve_fit(gaus_env, f,p, p0)

plt.figure(figsize=(8,6))
noisy_p = np.random.exponential(p)
plt.plot(f, noisy_p, label="exponential noise", alpha=0.5)
#plt.plot(f, p, label="model")
plt.plot(f, gaus_env(f,4,popt[1]-80, popt[2]-50), '-', linewidth=2)
plt.xlabel(r"Frequency ($\mu$Hz)")
plt.ylabel("Power spectral density (arb. units)")
plt.xlim((1000.,3300.))
ymax = np.max(noisy_p)+1
plt.ylim((0.,ymax))
#plt.legend()
#numax
plt.axvline(x=popt[1]-80, ymin=0, ymax=np.max(gaus_env(f,4,popt[1]-80, popt[2]-50))/ymax,linestyle='--', c='xkcd:orange', linewidth=2)
plt.annotate(r'$\nu_{\rm max}$',xy=(popt[1]-150,np.max(gaus_env(f,4,popt[1]-80, popt[2]-50))+0.5), color='xkcd:orange', fontsize=12)

#deltanu
#plt.arrow(modes[12,2],2,modes[12,2]-modes[13,2],0, color='xkcd:orange', width=0.05, head_length=10);
plt.annotate('', xy=(modes[13,2]+10,2), xytext=(modes[12,2]-10,2), arrowprops=dict(arrowstyle='<->', color='xkcd:orange'), va='center')

highest=7
for row in modes:
    l, n, freq, width, power = row
    if freq>1500. and freq<2800.:
        if l==0:
            plt.annotate('0',xy=(freq,highest), color='g')
        elif l==1:
            plt.annotate('1',xy=(freq,highest-0.5), color='r')
        elif l==2:
            plt.annotate('2',xy=(freq,highest-1), color='b')
        elif l==3 and freq>1800:
            plt.annotate('3',xy=(freq,highest-1.5), color='k')

plt.show()