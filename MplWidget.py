
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QLabel
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

#from alexfit import alexfit
#from stanfit import stanfit

class MyMplWidget(FigureCanvas):
    def __init__(self, vis, parent=None, figsize=(16,9), dpi=100):
        self.vis = vis
        self.fig = plt.figure(figsize=figsize, dpi=dpi)
        self.canvas = FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        #self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def plot_data(self, factor=1.7):
        self.vis.plot_visibilities(fig=self.fig)

    def plot_mode(self, inclination, del_axis=True):
        self.mx0.clear()
        self.vis.plot_model(0, self.vis.calc_visibilities(0, inclination), self.mx0, fig=self.fig)
        self.mx1.clear()
        self.vis.plot_model(1, self.vis.calc_visibilities(1, inclination), self.mx1, fig=self.fig)
        self.mx2.clear()
        self.vis.plot_model(2, self.vis.calc_visibilities(2, inclination), self.mx2, fig=self.fig)
        plt.xlabel('Frequency')
        FigureCanvas.draw(self)

    def plot_mode_start(self):
        self.mx0 = self.fig.add_subplot(self.vis.numls,2,2)
        self.vis.plot_model(0, self.vis.calc_visibilities(0, 90.), self.mx0, fig=self.fig)
        self.mx1 = self.fig.add_subplot(self.vis.numls,2,4)
        self.vis.plot_model(1, self.vis.calc_visibilities(1, 90.), self.mx1, fig=self.fig)
        self.mx2 = self.fig.add_subplot(self.vis.numls,2,6)
        self.vis.plot_model(2, self.vis.calc_visibilities(2, 90.), self.mx2, fig=self.fig)
        plt.xlabel('Frequency')



    '''def on_click(self, event):
        self.ax.scatter(event.xdata, event.ydata, c='r')
        self.plot_guess(event.xdata)
        self.draw()

    def fit_click(self):
        self.af.mlefit()
        self.plot_timeseries_model()'''
