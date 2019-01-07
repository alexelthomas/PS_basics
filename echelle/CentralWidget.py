from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QLabel, QSlider
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from MplWidget import MyMplWidget

import pandas as pd

class MyCentralWidget(QWidget):

    def __init__(self, main_window, ech, dnu_start):
        super().__init__()
        self.main_window = main_window
        self.ech = ech
        self.dnu_start = dnu_start
        self.initUI()

    def initUI(self):
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setFocusPolicy(Qt.StrongFocus)
        self.slider.setTickPosition(QSlider.TicksBothSides)
        self.slider.setMinimum(2*(self.dnu_start-15))
        self.slider.setMaximum(2*(self.dnu_start+15))
        self.slider.setValue(2*(self.dnu_start))
        self.slider.setTickInterval(10)
        self.slider.setSingleStep(1)

        self.mpl_widget = MyMplWidget(self.ech)
        # define label
        self.label = QLabel(self)
        # Place the buttons - HZ
        # Place the slider
        hbox = QHBoxLayout()
        hbox.addStretch(0.1)
        hbox.addWidget(self.slider)
        hbox.addStretch(0.01)
        self.slider.valueChanged.connect(self.valuechange)

        # place hbox and label into vbox
        vbox = QVBoxLayout()
        vbox.addWidget(self.mpl_widget)
        vbox.addLayout(hbox)
        self.setLayout(vbox)
        self.mpl_widget.plot_echelle_start(self.dnu_start)

    def valuechange(self):
      self.chosen_dnu = 0.5*self.slider.value()
      self.mpl_widget.plot_new_echelle(self.chosen_dnu)

    def on_finished_button_clicked(self):
        self.main_window.close()
