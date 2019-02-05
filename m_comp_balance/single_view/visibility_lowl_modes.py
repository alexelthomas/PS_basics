#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 16:27:04 2018

@author: axt367
"""
import numpy as np
import matplotlib.pyplot as plt


class visibility():
    def __init__(self, ls):
        self.ls = ls
        self.inclinations = np.arange(0,90,0.5)
        self.numls = len(self.ls)

    def calc_indv_visibility(self,l, m, inclination):
        inclination = np.radians(inclination)
        if l==0:
            if m==0:
                return 1
            else:
                return np.nan
        elif l==1:
            if m==0:
                return (np.cos(inclination))**2
            elif np.abs(m)==1:
                return 0.5*(np.sin(inclination))**2
            else:
                return np.nan
        elif l==2:
            if m==0:
                return (0.5*(3*np.cos(inclination)**2-1))**2
            elif np.abs(m)==1:
                return 1./6.* (3*np.sin(inclination)*np.cos(inclination))**2
            elif np.abs(m)==2:
                return 1./24.* (3*np.sin(inclination)**2)**2

    def calc_visibilities(self, l, inclination):
        m0 = self.calc_indv_visibility(l,0,inclination)
        m1 = self.calc_indv_visibility(l,1,inclination)
        m2 = self.calc_indv_visibility(l,2,inclination)
        return [m0, m1, m2]

    def plot_visibilities(self, fig=[], incl=[]):
        if fig == []:
            fig = plt.figure(figsize=(16,9), dpi=100)
        for idxl, l in enumerate(self.ls):
            ax = fig.add_subplot(self.numls,2,2*l+1)
            ax.plot(self.inclinations, [self.calc_indv_visibility(l,0,i) for i in self.inclinations], '-')
            ax.plot(self.inclinations, [self.calc_indv_visibility(l,1,i) for i in self.inclinations], '--')
            ax.plot(self.inclinations, [self.calc_indv_visibility(l,2,i) for i in self.inclinations], '-.')
            if not incl == []:
                ax.axvline(incl, c='k', ls=':', alpha=0.7)
            ax.set_ylabel('l='+str(l)+' visibilities')
            ax.set_xlim(0,90)
        plt.xlabel('Inclination angle (degrees)')

    def lorentzian(self, params):
        """
        Compute Lorentzian profile given input parameters
        """
        height, width, c_freq = params
        return height / (1.0+ (4.0 / width**2)*(self.freqs - c_freq)**2)

    def basic_model(self,l,visibilities):
        height, width, c_freq, splitting = 1., 2., 15., 5.
        self.freqs = np.linspace(0.,30.,1000)
        model = np.zeros(len(self.freqs))
        mfreqs = []
        mpows = []
        for m in np.linspace(0,l,l+1):
            if m==0:
                peak = self.lorentzian([height*visibilities[int(m)], width, c_freq])
                model += peak
                mfreqs = np.append(mfreqs, c_freq)
                mpows = np.append(mpows, np.max(peak))
            else:
                peak = self.lorentzian([height*visibilities[int(m)], width, c_freq-m*splitting])
                model += peak
                mfreqs = np.append(c_freq-m*splitting, mfreqs)
                mpows = np.append(np.max(peak), mpows)
                peak = self.lorentzian([height*visibilities[int(m)], width, c_freq+m*splitting])
                model += peak
                mfreqs = np.append(mfreqs,c_freq+m*splitting)
                mpows = np.append(mpows, np.max(peak))
        return model, mfreqs, mpows

    def plot_model(self, l, visibilities, mx, fig=[]):
        model, mfreqs, mpows = self.basic_model(l, visibilities)
        if fig == []:
            fig = plt.figure(figsize=(16,9), dpi=100)
        mx.plot(self.freqs, model, 'k-', alpha=0.8)
        ms = np.linspace(-l,l,2*l+1)
        for i, m in enumerate(ms):
            #print(mfreqs[i], np.where(self.freqs==mfreqs[i]))
            mx.annotate(str(int(m)), (mfreqs[i]-0.5, mpows[i]+0.1))#model[np.where(self.freqs==mfreqs[i])]))
        mx.set_ylabel('l=' + str(l) + ' modelled power')
        mx.set_ylim([0, 1.3])
        mx.set_xticks([])




if __name__=='__main__':
    ls = [0,1,2]
    inclination = 45.0
    vis = visibility(ls)
    print("Plotting visibilities")
    fig = plt.figure(figsize=(9,5), dpi=100)
    vis.plot_visibilities(fig, incl=inclination)
    
    ax0 = fig.add_subplot(len(ls), 2,2)
    vis.plot_model(0,vis.calc_visibilities(0,inclination), ax0, fig=fig)
    ax1 = fig.add_subplot(len(ls), 2,4)
    vis.plot_model(1,vis.calc_visibilities(1,inclination), ax1, fig=fig)
    ax2 = fig.add_subplot(len(ls), 2,6)
    vis.plot_model(2,vis.calc_visibilities(2,inclination), ax2, fig=fig)
    ax2.set_xlabel('Frequency')
    
    plt.tight_layout()
