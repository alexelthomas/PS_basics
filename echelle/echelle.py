#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 12:21:10 2019

@author: alex
"""
import numpy as np
import matplotlib.pyplot as plt

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QLabel
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from MyMainWindow import MyMainWindow

class Echelle():
    '''Produces Echelle diagram for specified l'''
    def __init__(self, ls_list, df_freqs):
        self.ls_list = ls_list
        self.df_freqs = df_freqs.where(df_freqs['l'].isin(ls_list)).dropna(0,'any')
        
    def get_freqs (self, l):
        ''' Picks out rows which correspond to l. '''
        return self.df_freqs.where(self.df_freqs['l']==l).dropna(0,'any')
    
    def echelle_prep(self, chosen_freqs, dnu):
        return chosen_freqs['freq']%dnu
    
    def echelle_indvplot(self, x, y, label, mstyle, ax):
        ax.plot(x,y, 'o', label=label, marker=mstyle)
        
    def plot_echelle(self, dnu, ax, fig=[]):
        style = ["d","s","^","o"]
        if fig == []:
            fig = plt.figure(figsize=(16,9), dpi=100)
            
        for idl, l in enumerate(self.ls_list):
            chosen_freq = self.get_freqs(l)
            self.echelle_indvplot(self.echelle_prep(chosen_freq, dnu), chosen_freq['freq'], label='l='+str(l), mstyle=style[idl], ax=ax)
        
        ax.legend(loc=0)
        ax.set_xlim([0.,180.])
        ax.set_ylabel(r'Frequency ($\mu$Hz)', fontsize=13)
        ax.set_xlabel(r'Frequency mod $\Delta \nu$ ($\mu$Hz)', fontsize=13)
        ax.set_title(r'Echelle calculated using $\Delta \nu = $' +str(dnu) + '($\mu$Hz)', fontsize=13)