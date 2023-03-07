from pathlib import Path
from PyQt5 import QtGui, QtCore, QtSvg
from pyqtgraph.Qt import QtWidgets
from PyQt5.QtWidgets import QMessageBox, QTableWidgetItem
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph.parametertree import Parameter, ParameterTree
from pyqtgraph import GraphicsLayoutWidget
import sys, os
import anndata
import pandas as pd
import numpy as np
from cnVisClass import *
from cnVisHelpers import annotatedColor 
import matplotlib.pyplot as plt
import scanpy as sc
import neurom as morphor_nm 
from neurom.features.utilities import (getSomaStats, extractMorhporFeatures)
from neurom.core.morphology import Morphology
from neurom import viewer as morphor_viewer
from morphio import SomaType
import glob
import warnings
import time
warnings.simplefilter('ignore')
## all data
homeDir = 'data'
tenX_meta_fn = os.path.join(homeDir,"snRNA-seq", "neurons","cn_neurons.meta.csv")
tenXdata_fn = os.path.join(homeDir, 'snRNA-seq', 'neurons', 'cn_neurons.h5ad')
patseq_meta_fn = os.path.join(homeDir,"Patchseq", "patchseq.dataset.final.meta.csv")
patseq_fn = os.path.join(homeDir,"Patchseq","patchseq.dataset.final.h5ad")
morphor_fn = os.path.join(homeDir,"Morphology")
morphor_files = glob.glob(os.path.join(homeDir,"Morphology","*.asc"))
morphor_names = [f.split('\\')[-1].split('.')[0] for f in morphor_files]
# patchseqFile = "data\\patchseq.dataset.final.h5ad"
# tenXUmapFile = "data\\tenxumap.csv"
# tenXUmap = np.array(pd.read_csv(tenXUmapFile,index_col=0))
tenX_metaDF = pd.read_csv(tenX_meta_fn, index_col=0)
ref_x,ref_y = tenX_metaDF['UMAP_embedding1'].values, tenX_metaDF['UMAP_embedding2'].values
qtCreatorFile = "resources\\uis\\cnVis.ui" # QT disigner file
Ui_MainWindow, QtBaseClass = pg.Qt.loadUiType(qtCreatorFile)
adata_pq = anndata.read_h5ad(patseq_fn)
patseq_umapX,patseq_umapY = adata_pq.obsm['X_umap'][:,0],adata_pq.obsm['X_umap'][:,1]
metaDF = adata_pq.obs.copy()
metaDF = annotatedColor(adata_pq, metaDF, 'T-cluster') ## add different color scheme if not there
metaDF = metaDF.rename(columns={"Cell": "Cell_ID"})
## selected metainformation for transcitome 
transcritMeta = ['Cell_ID','n_genes_by_counts', 'total_counts', 'total_counts_ERCC',
       'pct_counts_ERCC', 'T-cluster', 'T-cluster_umapX', 'T-cluster_umapY',
        'Seq', 'Expert_label', 'recording', 'Glyt2_eGFP',
       'Sst_tdTomato', 'Penk_tdTomato', '>10 gene number', '>0 gene number',
       'Seq_Depth', 'STAR_mapping_rate', 'Batch', 'UMAP1',
       'UMAP2']

ephyMeta =['Cell_ID','Tau (ms)',
       'SagRatio', 'Rebound (mV)', 'RM (Mohm)', 'RMP (mV)', 'RheoBase (pA)',
       'Spike Delay (ms)', 'Threshold (mV)', 'Amplitude (mV)', 'AHP (mV)',
       'Depolarization Time (ms)', 'Repolarization Time (ms)',
       'Half Width (ms)', 'Max Spike Number', 'AP Number @ 2xRheabase',
       'Initial Adaptation', 'Last Adaptation', 'AP2_Diff', 'AP3_Diff',
       'AP_End_Diff', 'Rebound_abs (mV)', 'Rebound_spikes',
       'Spike_Delay @ 2xRheobase', 'epsc_rise_time (ms)',
       'epsc_half_width (ms)', 'epsc_decay_tau (ms)', 'epsc_freq (Hz)',
       'epsc_amp (pA)', 'ipsc_rise_time (ms)', 'ipsc_half_width (ms)',
       'ipsc_decay_tau (ms)', 'ipsc_freq (Hz)', 'ipsc_amp (pA)']

ephyTracesDir = "data\\Electrophysiology" ## raw traces
morphorlogyDir = "data\\Morphology" ## morphology ASC files

class cnVis(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self,app, parent = None):
        super(cnVis, self).__init__(parent)
        self.app = app
        Ui_MainWindow.__init__(self)

        self.setupUi(self)
        self.setWindowTitle("Cochlear nucleus patch-seq data viewer")
        self.setWindowIcon(pg.QtGui.QIcon('resources\\icons\\ear.jpg'))
        self.statusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)
        self.metaDF = self.get_filtered_DF()
        self.init_comboBox_groupBy()
        self.init_transcritome_plots()
        self.init_morphorlogy_plots()
        self.init_ephy_plots()
        self.init_listWidget_ct()
        self.init_listWidget_cells()
        self.init_listWidget_sweeps()
        self.init_tableWidgets_metaData()
        self.pushButton_show_gene.clicked.connect(self.geneUmap) ## show gene expression
        self.cb_hasMorph.toggled.connect(self.filter_BlockChanged)
        self.checkBox_goodEphy.toggled.connect(self.filter_BlockChanged)
        self.horizontalSlider_mappingrate.valueChanged.connect(self.update_mappingrate)
        self.showMaximized()
        

    def init_comboBox_groupBy(self):
        ''' group cells according to user's choice.
        Default: annotated cell type
        '''
        self.gbLabels = {'Annotated':'annoated_10x_clusters',
        'Expert labeled':'Expert_label',}
        # self.ct_groupby()
        for t in self.gbLabels.keys():
            self.comboBox_groupBy.addItem(t)
        self.comboBox_groupBy.currentIndexChanged.connect(self.gbCBoxChange)

    @property
    def getCurrentGb(self):
        selGbTxt  = self.comboBox_groupBy.currentText()
        return self.gbLabels[selGbTxt]

    def gbCBoxChange(self):
        self.init_listWidget_ct()
        print(self.getCurrentGb)

    def update_mappingrate(self):
        self.spinBox_mappingrate.setValue(self.horizontalSlider_mappingrate.value())
        self.filterCellType()

    def resetItems(self):
        self.listWidget_sweeps.setCurrentItem(None)
        self.clearFigure(self.FP_plot)
        self.clearFigure(self.ephy_model_plot)
        self.clearFigure(self.morphor_plot)
        self.init_tableWidgets_metaData(clear=True)
        self.metaDF = metaDF.copy() ## A copy of the full dataset 


    def clearFigure(self,fig):
        fig.clf()
        fig.canvas.draw()

    def init_transcritome_plots(self):
        self.transcritome_plot = MatplotView()  ## 
        self.tabWidget_Transcitome.addTab(self.transcritome_plot,'Transcritome UMAP')
        self.transcritome_ref_plots()
        self.mapThisCTto10xUmap('All')

    def transcritome_ref_plots(self):
        self.transcritome_plot.clf()
        self.transcritome_plot.setSubplots(1, 2)
        ax1 = self.transcritome_plot.axes[0]
        sc.pl.umap(
            adata_pq,  s=25,ax=ax1,
            color=['T-cluster'], legend_fontsize=6,legend_fontoutline=0.1,title='',
            frameon=False, legend_loc= 'on data'
        )
        ax2 = self.transcritome_plot.axes[1]
        ax2.scatter(ref_x, ref_y, marker='.', c='gray', alpha=0.2, s=3)
        ax2.axis('off')
        plt.show()

    def geneUmap(self):
        ''' show gene expression in umap
        '''
        pass

    def init_ephy_plots(self):
        self.FP_plot  = MatplotView()
        self.PSPC_plot = MatplotView()
        self.ephy_model_plot = MatplotView()
        self.tabWidget_Ephy.addTab(self.FP_plot,'Firng pattern')
        self.tabWidget_Ephy.addTab(self.PSPC_plot ,'PSP/PSC')
        self.tabWidget_Ephy.addTab(self.ephy_model_plot ,'Ephy model')

    def init_morphorlogy_plots(self):
        self.morphor_plot = MatplotView()
        self.slice_plot = MatplotView()
        self.tabWidget_Morphor.addTab(self.morphor_plot,'Morphorlogy')
        self.tabWidget_Morphor.addTab(self.slice_plot,'Slice')

    def init_tableWidgets_metaData(self, clear=False):
        if clear:
            self.tableWidget_EphyMetaData.setRowCount(0)
            self.tableWidget_morphorMetaData.setRowCount(0) 
            self.tableWidget_TransMetaData.setRowCount(0)

        self._init_tableWidget_metaData(self.tableWidget_EphyMetaData,len(ephyMeta),130)
        self._init_tableWidget_metaData(self.tableWidget_morphorMetaData,20,featureColWidth=160) 
        self._init_tableWidget_metaData(self.tableWidget_TransMetaData,len(transcritMeta), featureColWidth=150)

    def _init_tableWidget_metaData(self, tablewidget, nrow,featureColWidth=120):
        tablewidget.setColumnWidth(0,featureColWidth)
        tablewidget.setColumnWidth(1,100)
        tablewidget.setRowCount(nrow)

    def init_listWidget_sweeps(self):
        self.listWidget_sweeps.itemDoubleClicked.connect(self.draw_FP)
        self.listWidget_sweeps.itemSelectionChanged.connect(self.draw_FP)

    def init_listWidget_ct(self):
        self.listWidget_ct.clear() 
        self.listWidget_ct.addItem('All')
        self.cellTypes = sorted(set(list(metaDF[self.getCurrentGb])))
        for ct in self.cellTypes:
            self.listWidget_ct.addItem(ct)
        self.listWidget_ct.setCurrentItem(self.listWidget_ct.item(0))
        self.listWidget_ct.itemDoubleClicked.connect(self.filterCellType)
        self.listWidget_ct.itemSelectionChanged.connect(self.filterCellType)
        
    def filterCellType(self):
        ''' show cell types specified in the cells list
        '''
        self.metaDF = self.get_filtered_DF()
        curCT = self.listWidget_ct.currentItem().text()
        cellIDs = self.load_cellID_list(curCT)
        self.listWidget_cells.clear()
        for c in cellIDs:
             self.listWidget_cells.addItem(c)
        self.resetItems()
        self.label_cellNum.setText(f'n={len(cellIDs)}')
        if hasattr(self, 'transcritome_plot'):
            self.mapThisCTto10xUmap(curCT)  
        self.statusBar.showMessage(curCT)

    def filter_BlockChanged(self):
        '''refilter dependent on check box state'''
        self.filterCellType()
        self.mapThisCTto10xUmap(self.listWidget_ct.currentItem().text())

    def filter_with_mappingRate(self, df):
        mpRateThr = self.horizontalSlider_mappingrate.value() ## lower threhold
        return df[df['STAR_mapping_rate']>=mpRateThr]

    def filter_with_morphorlogy(self, df):
        hasMorph = self.cb_hasMorph.isChecked()
        df_ = df.copy()
        if hasMorph:
            df_ = df[df['Cell_ID'].isin(morphor_names)] ## filter cells with morph
        return df_

    def filter_with_ephyRecording(self, df):
        goodEphy = self.checkBox_goodEphy.isChecked()
        if goodEphy:
            df = df[df.recording.isin(['ok','good'])]
        return df

    def get_filtered_DF(self):
        df = self.filter_with_morphorlogy(metaDF.copy())
        df = self.filter_with_mappingRate(df)
        df = self.filter_with_ephyRecording(df)
        return df     

    def load_cellID_list(self, curCT):

        if curCT == 'All':
            cellID = list(self.metaDF.index)
        else:
            cellID = list(self.metaDF[self.metaDF.annoated_10x_clusters==curCT].index)
        return cellID

    def init_listWidget_cells(self):
        self.listWidget_cells.itemDoubleClicked.connect(self._load_cell_data)
        self.listWidget_cells.itemSelectionChanged.connect(self._load_cell_data)
        self.filterCellType()

    def _load_cell_data(self):
        ''' load data when double click a Cell ID
        '''
        i = self.listWidget_cells.currentItem()
        if i is not None:
            t0 = time.time()
            self.refresh_cellMeta(i.text()) ## update meta table information
            self.refresh_umaps(i.text()) ## update UMAP      
            self.refresh_FP(i.text()) ## update firing pattern
            self.refresh_morphorlogy(i.text())
            statMsg = f'{self.listWidget_ct.currentItem().text()}: {i.text()}'
            self.statusBar.showMessage(statMsg)

    def refresh_morphorlogy(self, cellID):
        self.morphor_plot.clf()
        self.morphor_plot.setSubplots(1, 1)
        ax =  self.morphor_plot.axes[0]
        cellName = metaDF.loc[cellID, 'Cell_ID']
        ascFile = glob.glob(os.path.join(morphorlogyDir, cellName +".asc"))
        if ascFile==[]:
            ascFile = glob.glob(os.path.join(morphorlogyDir, cellName +".ASC"))
            if ascFile==[]: ## could be either .asc or.ASC
                ax.set_title('No morphorlogy')
                self.tableWidget_morphorMetaData.setRowCount(0) ## erase previous one
        else:
            ascFile=ascFile[0]
            # neurons, custom_data = morphor_nm.load_neuron(ascFile) ## neurm_pv version
            neuron = morphor_nm.load_morphology(ascFile, somaType = SomaType.SOMA_CYLINDERS)
            morphor_viewer.draw(
                    neuron, mode="2d", realistic_diameters=True, ax=ax)
            ax.set_title('')
            self.refresh_morphTable(neuron)

        ax.axis("equal") 
        ax.axis('off')
        self.morphor_plot.getFigure().tight_layout()
        self.morphor_plot.canvas.draw()


    def refresh_cellMeta(self, cellID):
        # print(f'cell {cellID} loaded!')
        self.refresh_MetaTable(cellID, self.tableWidget_EphyMetaData, ephyMeta)
        self.refresh_MetaTable(cellID, self.tableWidget_TransMetaData, transcritMeta)
        ## to be integrated into meta
        # self.refresh_MetaTable(cellID, self.tableWidget_morphorMetaData, morphorMeta)

    def refresh_morphTable(self, neuron:Morphology) -> None :
        df = extractMorhporFeatures(neuron)
        features = list(df.keys())
        widget = self.tableWidget_morphorMetaData
        widget.setRowCount(0)
        widget.setRowCount(len(features))
        ridx = 0
        for feature in features:
            value = df[feature]
            widget.setItem(ridx, 0, QtWidgets.QTableWidgetItem(feature))
            if value is None:
                value = 'nan'
            elif not isinstance(value, str):
                value = str(np.round(value,3))
            widget.setItem(ridx, 1, QtWidgets.QTableWidgetItem(value))
            ridx +=1 

    def refresh_MetaTable(self, cellID, widget, features):
        ridx = 0
        for feature in features:
            value = str(metaDF.loc[cellID, feature])
            widget.setItem(ridx, 0, QtWidgets.QTableWidgetItem(feature))
            widget.setItem(ridx, 1, QtWidgets.QTableWidgetItem(value))
            ridx +=1 

    def getEphyMeta(self, cellID):
        ''' Ephy meta can be read from original txt/csv file. Or
            from intergrated andata object: refresh_ephyMeta()
        '''
        self.ephyMeta = metaDF[metaDF.Cell_ID==cellID]
        with  open(os.path.join(ephyTracesDir,cellID+".asc"), "r") as f:
            for c in ['Holding(pA)','current_step','i0']:
                line = f.readline()
                self.ephyMeta[c] = float(line.strip().split(':')[-1])
        self.getStimPars() ## extrat stim parameters

    def getStimPars(self):
        rheoBase = self.ephyMeta['RheoBase (pA)'].values[0]
        self.i0 = self.ephyMeta['i0'].values[0]  ## initial current
        self.stepI = self.ephyMeta['current_step'].values[0] ## current step
        ## rheobase is esmated thus may exceed actuall sweep number
        try:
            self.rheobaseSweepIdx = min(np.ceil((rheoBase-self.i0)//self.stepI), self.fp.shape[1]-1)
        except:
            print(f'rheobase:{rheoBase}, i0:{self.i0},step:{self.stepI},{np.ceil((rheoBase-self.i0)//self.stepI)}')
        self.sampleTime = np.diff(self.fp.index)[0] ## sampling perirod
        self.timeX = np.arange(len(self.fp))*self.sampleTime*1000 ## stimuli time in ms

    def refresh_FP(self, cellID):
        cellID = metaDF.loc[cellID, 'Cell_ID'] ## ephy cell ID
        try:
            self.fp = pd.read_csv(os.path.join(ephyTracesDir,cellID+".asc"), sep='\t',
                                  header=0,index_col=0, skiprows=3)
            self.fp.columns = range(1,self.fp.shape[1]+1)
        except:
            print(f'file not found for cell {cellID}')
            return 
        sweepList = [str(c) for c in list(self.fp.columns)]
        self.getEphyMeta(cellID)
        self.listWidget_sweeps.clear()
        for c in sweepList:
            self.listWidget_sweeps.addItem(c)
        self.label_sweepNum.setText(f'n={len(sweepList)}')
        self.plotSweeps([1, self.rheobaseSweepIdx , self.fp.shape[1]])
        
    def draw_FP(self):
        ''' Dispatcher to draw firing pattern
        '''
        if hasattr(self, 'fp'):
            i = self.listWidget_sweeps.currentItem()
            if i is None: ## draw tree repr. sweeps if none
                self.plotSweeps([1, self.fp.shape[1]//2, self.fp.shape[1]])
            else:
                self.plotSweeps(int(i.text()))

    def plotSweeps(self, i):
        ''' plot sweeps of response and stimulus
            if provide an plotting axes, make it contain at least two subplots panel.
            By default, stimulus will be in the second panel.
        '''
        if type(i) is int:
            i = [i]
        assert type(i) is list, "index to be integer or list of integers"
        self.FP_plot.clf()
        self.FP_plot.setSubplots(2, 1,heightRatio=[4,1])
        fpAxe = self.FP_plot.axes[0]
        stimAxe = self.FP_plot.axes[1]
        for s in i: ## plot all series
            stimulus, response = self.getSweep(int(s))
            fpAxe.plot(self.timeX, response)
            stimAxe.plot(self.timeX, stimulus)
        fpAxe.set_xlabel('')
        fpAxe.spines['top'].set_visible(False)
        fpAxe.spines['bottom'].set_visible(False)
        fpAxe.spines['right'].set_visible(False)
        fpAxe.set_ylabel('mV')
        fpAxe.grid(visible=True, color='k', linestyle='--',alpha=0.3, zorder=1000)
        fpAxe.get_xaxis().set_ticks([])
        stimAxe.set_ylabel('pA')
        stimAxe.set_xlabel('time (mS)')
        stimAxe.spines['top'].set_visible(False)
        stimAxe.spines['right'].set_visible(False)
        stimAxe.grid(visible=True, color='k', linestyle='--',alpha=0.3, zorder=1000)
        self.FP_plot.getFigure().tight_layout()
        self.FP_plot.canvas.draw()

    def getSweep(self,sweepIdx):
        ''' extract response and stim
        '''
        assert hasattr(self,'fp'), 'Firing pattern data not available!'
        resp = self.fp.loc[:,sweepIdx].values
        N = len(self.fp)
        initIdx = int(0.1/self.sampleTime) ## usually 0.1 s 
        stimIdx = int(0.6/self.sampleTime)
        stim = np.ones((N,))*self.ephyMeta['Holding(pA)'].values[0]
        stim[initIdx: initIdx+stimIdx] = self.i0 + (sweepIdx-1)*self.stepI## holding current
        return stim, resp

    def refresh_umaps(self, cellID, removePrevousDot=True, s=30): 
        ''' add highlight to currently selected cell
        '''
        if hasattr(self, 'pqplot_dot') and removePrevousDot:
            try:
                self.pqplot_dot.remove() ## remove previous dot in patchseq reference
                self.tenxplot_dot.remove() ## remove previous dot in tenx reference
            except ValueError:
                pass# print('Nothing to be removed')
        ax1 = self.transcritome_plot.axes[0]
        coords = np.array(adata_pq[cellID].obsm['X_umap'])[0]
        self.pqplot_dot = ax1.scatter(coords[0],coords[1],marker='o', color='r', s=s)
        ax2 = self.transcritome_plot.axes[1]
        x= metaDF.loc[cellID]['10x_mapping_umapX']
        y= metaDF.loc[cellID]['10x_mapping_umapY']
        self.tenxplot_dot = ax2.scatter(x,y,marker='o', color='r', s=s)
        self.transcritome_plot.canvas.draw()

    def mapThisCTto10xUmap(self, curCT):
        ''' show mapped location for specified cell type
        '''
        if hasattr(self, 'pqplot_dot') :
            try:
                self.pqplot_dot.remove() ## remove previous dot in patchseq reference
            except ValueError:
                pass# print('No dots to be removed in patchseq UMAP')
        if hasattr(self, 'tenxplot_dot'):
            try:
                self.tenxplot_dot.remove() ## remove previous dot in tenx reference
            except ValueError:
                pass# print('No dots to be removed in 10x UMAP')
        
        if curCT=='All':
            df = self.metaDF
        else:
            df = self.metaDF[self.metaDF['annoated_10x_clusters']==curCT]
        cellIDs = list(df.index)
        ax1 = self.transcritome_plot.axes[0]
        coords = np.array(adata_pq[adata_pq.obs_names.isin(cellIDs)].obsm['X_umap'])[0]
        self.pqplot_dot = ax1.scatter(coords[0],coords[1],marker='o', color='r', s=10)
        ax2 = self.transcritome_plot.axes[1]
        x, y = df['10x_mapping_umapX'], df['10x_mapping_umapY']
        c = df['T-cluster_colors']
        self.tenxplot_dot = ax2.scatter(x,y,marker='o', color=c, s=6)
        self.transcritome_plot.canvas.draw()

    def closeApp(self):
        sys.exit()


def main():
    app = QtWidgets.QApplication(sys.argv)
    main = cnVis(app)
    main.show()
    sys.exit(app.exec_())
  
if __name__ == '__main__':
    main()
