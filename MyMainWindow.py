
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QLabel
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from CentralWidget import MyCentralWidget

class MyMainWindow(QMainWindow):
    def __init__(self, vis):
        super().__init__()
        self.vis = vis
        self.initUI()

    def initUI(self):
        self.resize(1600,900)
        self.move(50,50)
        central_widget = MyCentralWidget(self, self.vis)
        self.setCentralWidget(central_widget)

        self.setWindowTitle('Background Visualisation')
        self.statusBar().showMessage('Waiting ...')
