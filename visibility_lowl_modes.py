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

    def plot_visibilities(self, fig=[]):
        if fig == []:
            fig = plt.figure(figsize=(16,9), dpi=100)
        for idxl, l in enumerate(self.ls):
            ax = fig.add_subplot(self.numls,2,2*l+1)
            ax.plot(self.inclinations, [self.calc_indv_visibility(l,0,i) for i in self.inclinations], '-')
            ax.plot(self.inclinations, [self.calc_indv_visibility(l,1,i) for i in self.inclinations], '--')
            ax.plot(self.inclinations, [self.calc_indv_visibility(l,2,i) for i in self.inclinations], '-.')
            ax.set_ylabel('l='+str(l)+' visibilities')
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
        for m in np.linspace(0,l,l+1):
            if m==0:
                model += self.lorentzian([height*visibilities[int(m)], width, c_freq])
            else:
                model += self.lorentzian([height*visibilities[int(m)], width, c_freq-m*splitting])
                model += self.lorentzian([height*visibilities[int(m)], width, c_freq+m*splitting])
        return model

    def plot_model(self, l, visibilities, mx, fig=[]):
        model = self.basic_model(l, visibilities)
        if fig == []:
            fig = plt.figure(figsize=(16,9), dpi=100)
        mx.plot(self.freqs, model, '-')
        mx.set_ylabel('l=' + str(l) + ' modelled power')
        mx.set_ylim([0, 1])




if __name__=='__main__':
    ls = [0,1,2]
    vis = visibility(ls)
    print("Plotting visibilities")
    vis.plot_visibilities()
    plt.show()
