from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QLabel, QSlider
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from MplWidget import MyMplWidget

import pandas as pd

class MyCentralWidget(QWidget):

    def __init__(self, main_window, vis):
        super().__init__()
        self.main_window = main_window
        self.vis = vis
        self.initUI()

    def initUI(self):
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setFocusPolicy(Qt.StrongFocus)
        self.slider.setTickPosition(QSlider.TicksBothSides)
        self.slider.setMinimum(0)
        self.slider.setMaximum(90)
        self.slider.setValue(90)
        self.slider.setTickInterval(10)
        self.slider.setSingleStep(1)

        self.mpl_widget = MyMplWidget(self.vis)
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
        self.mpl_widget.plot_data()
        self.mpl_widget.plot_mode_start()

    def valuechange(self):
      self.chosen_inclination = self.slider.value()
      self.mpl_widget.plot_mode(self.chosen_inclination)

    def on_finished_button_clicked(self):
        self.main_window.close()
