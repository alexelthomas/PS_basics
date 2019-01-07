#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd
from astropy.convolution import convolve, Box1DKernel

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QLabel
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

app = None

#from keplerdata import Dataset

from MyMainWindow import MyMainWindow

def get_data(ns_list=[]):
    if not ns_list:
        ns_list = np.arange(12,26,1)
        
    df = pd.read_csv('~/Documents/PhD/Projects/separation_ratios/data/simulated_freqs/pristine', sep='\s+', header=0, names=['l','n','freq','freqerr'])
    df = df.loc[df['l']!=3]
    df = df.where(df['n'].isin(ns_list)).dropna(axis=0, how='any')
    return df

def main(ech, dnu_start, verbose=False):
    '''
    app must be defined already!!!
    '''
    global app
    if verbose:
        print('Setting up app')
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    if verbose:
        print('Setting up MyMainWindow')
    w = MyMainWindow(ech, dnu_start)
    w.show()
    if verbose:
        print('Exiting')
    app.exit(app.exec_())

if __name__ == "__main__":
    from echelle import Echelle

    df_freqs = get_data()
    ls_list = ['0','1','2']
    ech = Echelle(ls_list, df_freqs)
    dnu_start = 137.0

    print('Now running main ...')
    main(ech, dnu_start, verbose=True)
