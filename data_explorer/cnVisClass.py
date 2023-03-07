from pyqtgraph.widgets import MatplotlibWidget
from pyqtgraph import GraphicsLayoutWidget
import pyqtgraph as pg
import numpy as np
class MatplotView(MatplotlibWidget.MatplotlibWidget):
    """
    Matplotlib figure widget
    """
    def __init__(self, parent=None, size = (5.0, 4.0), title = '', dpi=100, plotGrid = None):
        ''' plotGrid is a tuple with (nrow, ncol)
        '''
        super(MatplotView, self).__init__()
        self.setParent(parent)
        self.figure = self.getFigure()
        self.clf()
        self.grid=None
        self.axes=None
    
    def setSubplots(self, nrow, ncol, widthRatio=None, heightRatio=None):
        ## set up subplots
        self.nRow=nrow
        self.nCol=ncol
        if widthRatio is None:
            widthRatio = [1]*ncol
        if heightRatio is None:
            heightRatio = [1]*nrow
        axes = self.figure.subplots(nrow, ncol,\
             gridspec_kw={'width_ratios': widthRatio,'height_ratios':heightRatio})
        if type(axes) is not list:
            if type(axes) is np.ndarray:
                self.axes = list(axes.flatten())
            else:
                self.axes = [axes]
        # self.axes = [self.figure.add_subplot(nrow, ncol, j) for j in np.arange(1, nrow*ncol+1)]
    
    def setPlotGrid(self, nRow, nCol):
        self.grid = self.figure.add_gridspec(nRow, nCol)

    def plotline(self, i, j,  x, y, color = 'r', title = '', linewidth = 0.5):
        if i >= self.grid.nrows or j >= self.grid.ncols:
            raise ValueError("row index or col index out of bound")
        ax = self.axes()   
        ax.plot(x, y, color = color, linewidth = linewidth)
        ax.set_title(title, fontsize = 12, color = 'k')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        return ax
    
    def adjustSubplot(self, left = 0.08, bottom=0.1, top=0.95, right = 0.95, wspace = 0.35, hspace = 0.35):
        self.figure.subplots_adjust(left = left, bottom = bottom, \
                                    top= top, right = right, wspace = wspace, hspace = hspace)
        
    def clf(self):
        self.figure.clear()

class PlotView(GraphicsLayoutWidget):
    """
    Class for plot view.
    """
    def __init__(self, parent, title = '', background = 'k'):
        super(PlotView, self).__init__()
        self.frame = parent
        self.setWindowTitle(title)
        self.setBackground(background)
        self.ci.setBorder()
        self.ci.setBorder(pg.mkPen(None))
        self.ci.setContentsMargins(10, 10, 10, 20)
        self.setCentralItem(self.ci)
        self.show()
        
    def refresh(self):
        self.ci.setContentsMargins(10, 10, 10, 20)
