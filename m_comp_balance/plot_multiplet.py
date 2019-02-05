#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 16:27:04 2018

@author: axt367
"""
import numpy as np
import matplotlib.pyplot as plt

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QLabel
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from MyMainWindow import MyMainWindow

class visibility():
    def __init__(self, l, inclination):
        self.l = l
        self.inclination = np.radians(inclination)
        self.height = 1.
        self.width = 1.0
        self.c_freq = 0.
        self.splitting = 0.4
        self.asymmetry = 0.25
        

    def calc_indv_visibility(self, m):
        if self.l==0:
            if m==0:
                return 1
            else:
                return np.nan
        elif self.l==1:
            if m==0:
                return (np.cos(self.inclination))**2
            elif np.abs(m)==1:
                return 0.5*(np.sin(self.inclination))**2
            else:
                return np.nan
        elif self.l==2:
            if m==0:
                return (0.5*(3*np.cos(self.inclination)**2-1))**2
            elif np.abs(m)==1:
                return 1./6.* (3*np.sin(self.inclination)*np.cos(self.inclination))**2
            elif np.abs(m)==2:
                return 1./24.* (3*np.sin(self.inclination)**2)**2

    def calc_visibilities(self):
        m0 = self.calc_indv_visibility(0)
        m1 = self.calc_indv_visibility(1)
        m2 = self.calc_indv_visibility(2)
        return [m0, m1, m2]

    def lorentzian(self, params):
        """
        Compute Lorentzian profile given input parameters
        """
        height, width, c_freq = params
        return height / (1.0+ (4.0 / width**2)*(self.freqs - c_freq)**2)

    def basic_model(self):
        visibilities = self.calc_visibilities()
        self.freqs = np.linspace(-1.5,1.5,1000)
        model = np.zeros(len(self.freqs))
        m_indv = np.zeros((len(self.freqs), 2*l+1))
        for m in np.linspace(0,self.l,self.l+1):
            if m==0:
                m_indv[:,self.l] = self.lorentzian([self.height*visibilities[int(m)], self.width, self.c_freq])
                model += m_indv[:,self.l]
            else:
                m_indv[:,int(self.l-m)] = self.lorentzian([self.height*visibilities[int(m)], self.width, self.c_freq-m*self.splitting+self.asymmetry])
                model += m_indv[:,int(self.l-m)]
                m_indv[:,int(self.l+m)] = self.lorentzian([self.height*visibilities[int(m)], self.width, self.c_freq+m*self.splitting+self.asymmetry])
                model += m_indv[:,int(self.l+m)]
        return model, m_indv

    def plot_model(self,fig=[]):
        model, m_indv = self.basic_model()
        linestyles = ['-.','--', ':','-']
        if fig == []:
            fig,ax = plt.subplots(figsize=(7,5))
        ax.plot(self.freqs, model, 'k',ls=linestyles[-1])
        for m in range(2*self.l+1):
            ax.plot(self.freqs, m_indv[:,m], 'k', ls=linestyles[np.abs(m-self.l)])
            print(m)
            if (m-self.l) == 0:
                ax.axvline(x=self.c_freq-(m-self.l)*self.splitting, color='k',ls='-.', alpha=0.2)
            else:
                ax.axvline(x=self.c_freq-(m-self.l)*self.splitting, color='k',ls='--', alpha=0.2)
        ax.set_ylabel('Relative power', fontsize=13)
        ax.set_xlabel(r'Offset from central frequency ($\mu$Hz)', fontsize=13)
        ax.set_xticks([-1.5,-1,-0.5,0,0.5,1,1.5])
        #ax.set_yticks(np.linspace(0,2,4))
        ax.set_xlim(-1.5,1.5)
        print(np.min(model), np.max(model))




if __name__=='__main__':
    l=1
    inclination = 45.
    vis = visibility(l, inclination)
    print("Plotting visibilities")
    vis.plot_model()
    plt.tight_layout()
    plt.show()
