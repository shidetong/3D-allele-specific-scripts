import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import pyBigWig
#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')
                
def caxis_S_horizontal(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 23)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0)
    ax.spines['left'].set_linestyle('dotted')


def pca_To_20K(fil):
    """
    """
    pca_type = np.dtype({'names':['chr' , 'PCA'],
                     'formats':['S4' , np.float]})
    PCA_Data = np.loadtxt(fil , dtype=pca_type)
    
    chroms = set(PCA_Data['chr'])
    New_Data = {}
    for c in chroms:
        New_Data[c] = {}
        tmp_data = PCA_Data[PCA_Data['chr'] == c]
        New_Data[c] = []
        for i in tmp_data:
            New_Data[c].extend([i['PCA']] * 10)

    return New_Data

def UpdateDI(DI):
    """
    """
    New_DI = []
    New_index = []

    for index in range(len(DI) - 1):
        if DI[index] * DI[index + 1] < 0:
            New_DI.append(DI[index])
            New_DI.append(0)
            New_index.append(index)
            New_index.append(index + 0.5)
        else:
            New_DI.append(DI[index])
            New_index.append(index)
    
    return np.array(New_index), np.array(New_DI)

def bigwig_10bp_Plot(fil,chro,start,end,fig,location,color,label):
    """
    """
    sig_type = np.dtype({'names':['start' , 'end' , 'score'],
                      'formats':[np.int , np.int , np.float]})
    bw = pyBigWig.open(fil)
    bw = bw.intervals(chro, start, end)
    tmp_data = np.array(list(bw) , dtype = sig_type)
    bin_size = (end - start) // 10 + 1
    sig_data = np.zeros((bin_size,))
    for line in tmp_data:
        s = line['start'] // 10 - start // 10
        e = line['end'] // 10 - start // 10
        for i in range(s,e):
            if i >= 0 and i < bin_size: 
                sig_data[i] += line['score']
            else:
                pass
    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')
    ax.set_xlim((0 , len(sig_data)))
    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S_horizontal(ax,color)


PCA_Data =pca_To_20K('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/H6C7/Cooler/Traditional_PC/Traditional_PC_Compartment_500K.txt')
RNA_Data = '/public/home/shidetong/projects/yf/chip-seq/chip_2/bw/Bxpc_CTCF-1.bw'
selected_interval = [('4' , 55200000 , 57000000 , 'Klf4')]
R = 40000

for i in selected_interval:
    print (i)
    g = i[0]
    startHiC = i[1] // R 
    endHiC = i[2] // R 

    chro = 'chr'+g 

    size = (12, 12)
    Left = 0.25 ; HB = 0.05 ; width = 0.6 ; HH = 0.6
    PCAData = PCA_Data[g]
    PCA = np.array(PCAData[startHiC:endHiC])
    fig = plt.figure(figsize = size)
    ax1 = fig.add_axes([Left  , HB , width , HH])

    PCA_index , PCA = UpdateDI(PCA)
    ax7 = fig.add_axes([Left, HB + width + 0.05 + 0.075 , width , 0.075])
    ax7.fill_between(PCA_index , PCA , where = PCA >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
    ax7.fill_between(PCA_index , PCA , where = PCA <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
    ax7.set_xlim(0 , PCA_index.max())
    ax7.set_ylabel('PC1',fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S_horizontal(ax7, 'black') 

    ##RNA Tracks
    location = [Left , HB + HH + 0.05, width , 0.075]
    ax4 = fig.add_axes(location)
    bigwig_10bp_Plot(RNA_Data , chro , i[1] , i[2] , fig , location , 'fuchsia' , 'RNA')





