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

def main(vis, verbose=False):
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
    w = MyMainWindow(vis)
    w.show()
    if verbose:
        print('Exiting')
    app.exit(app.exec_())

if __name__ == "__main__":
    from visibility_lowl_modes import *

    ls = [0,1,2]
    vis = visibility(ls)

    print('Now running main ...')
    main(vis, verbose=True)
