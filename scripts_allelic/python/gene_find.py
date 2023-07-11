from __future__ import division
import math
import numpy as np
import csv , copy
# import xlrd
import pandas as pd
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['U64' , 'U64' , 'U8' , 'U4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf


def Get_interval_genes(interval , gtf ):
    '''
    '''
    genes = []
    for i in interval:
        g = i[0]
        start = i[1] 
        end = i[2]
        tmp = gtf[gtf['chr'] == "chr"+ g]
        mask = (tmp['start'] <= end) & (tmp['end'] >= start)
        overlap = tmp[mask]
        if overlap.size != 0:
            for j in overlap:
                genes.append(j)
            
            
    genes = np.array(genes , dtype = gtf.dtype)
    return genes



gtf = Load_gtf('/public/home/shidetong/wedata/test/genenum/gencode.v37lift37.annotation.gtf')
interval =[('22' , 22300000 , 23290555)]

genes = Get_interval_genes(interval , gtf)
