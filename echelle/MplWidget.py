
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QLabel
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

#from alexfit import alexfit
#from stanfit import stanfit

class MyMplWidget(FigureCanvas):
    def __init__(self, ech, parent=None, figsize=(16,9), dpi=100):
        self.ech = ech
        self.fig = plt.figure(figsize=figsize, dpi=dpi)
        self.canvas = FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        #self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def plot_new_echelle(self, dnu, del_axis=True):
        self.ax.clear()
        self.ech.plot_echelle(dnu, self.ax, fig=self.fig)
        FigureCanvas.draw(self)

    def plot_echelle_start(self, dnu_start):
        self.ax = self.fig.add_subplot(1,1,1)
        self.ech.plot_echelle(dnu_start, self.ax, fig=self.fig)