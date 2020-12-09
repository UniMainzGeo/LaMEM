# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 09:56:36 2020

@author: AndreaPiccolo
"""
# -*- coding: utf-8 -*-
# import scipy.io as sio
# import sys,os,fnmatch
# import numpy as np
# import matplotlib.pyplot as plt
# plt.switch_backend('agg')
import matplotlib
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import shutil
import time
import glob
import re
# from matplotlib.colors import BoundaryNorm
# from matplotlib.ticker import MaxNLocator
# from pylab import figure, axes, pie, title, show
# import pylab
#from matplotlib import rcParams
# import time
# from mpl_toolkits.mplot3d import Axes3D
#import json
from AUX_FUN import *
import matplotlib.colors as colors
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.ma import masked_array
import matplotlib.cbook as cbook
from  classes_visualization import * 
###############################################################################
# Main script
###############################################################################
# Folder Path
####################################

   # File list 
"""
LITTLE SUGGESTION. TRY TO NOT USE DIFFERENT NAME FOR THE PVD FILES. DEFINE YOUR 
TEST FOLDER NAME. 
"""
# If you follow the suggestion, change the name here.
fname=("Archean_drips")

####################################


dyn             ='%s.pvtr'      %fname
surf            ='%s_surf.pvtr' %fname
phase           ='%s_phase.pvtr'%fname
passive_tracers ='%s_passive_tracers.pvtu'%fname



DYN   = 1          #print dynamic properties
SURF  = 0          #print free surface properties 
PHASE = 1          #print Phase field
PTT   = 1          #print PTT_path


#Zoom and quiver Here you can define several type of zoom and contour. By default 
#quiver and 
####################################
contour1=contour()
contour1.yes='yes'
contour1.vector= np.linspace(0,1600,num=8,endpoint=True)
contour1.color='black'
contour1.text='yes'
contour1.width=0.5
######################################
quiver1=quiver()
quiver1.yes = 'yes'
quiver1.nodes = 6
quiver1.color='k'
######################################
zoom1 = zoom()
zoom1.yes = 'yes'
zoom1.zoom_coord=(1000, 3000, -700, 20)


streaml = stream_line()
streaml.yes = 'yes'
streaml.scale = ['constant',0.2]
streaml.density = 0.6

streaml2 = stream_line()
streaml2.yes = 'yes'
streaml2.field = "tot_displ [km]"
streaml2.scale = ['scale',0.5]
streaml2.density = 0.6
###########################
#PTT
PTT_show1=PTT_show()
PTT_show1.yes = 'yes'
       
######################################
#####################################################################
# Custom option -check visualization classes for the default option #
#####################################################################

strain_rate_ii.clim = [-19,-12]

Pressure.colormap = 'Seismic'
Pressure.zoom     = zoom1

Velocity.quiver = quiver1
Velocity.zoom   = zoom1
Velocity.component_show = ['mag','x','z']

Temperature.stream_line = streaml

Melt_Fraction.PTT_show = PTT_show1
Melt_Fraction.stream_line = streaml2



Melt_Extracted.PTT_show = PTT_show1
Melt_Extracted.stream_line = streaml2



#################################################################
#Customize classes                                              #
#################################################################


Prop=['Temperature','strain_rate_ii','stress_ii','Pressure','Density','Velocity','Visc_Creep','Visc_Total','Melt_Fraction','Melt_Extracted']  #Vector containing the properties that you want to visualize and print

############################################################################################
Phase=P_phase()
##https://matplotlib.org/3.1.0/gallery/color/named_colors.htmlhttps://matplotlib.org/3.1.0/gallery/color/named_colors.html
cmap=[  'w',
        'lavender',       #1
       'lightsteelblue',   #2
       'royalblue',   #3
       'navy',   #4
       'darkblue',   #5
       'indigo',   #6
       'greenyellow',   #7
       'darkseagreen',   #8
       'forestgreen',   #9
       'darkgreen',   #10
       'black',   #11
       'gainsboro',   #12
       'goldenrod',    #13
       'palevioletred',     #14
       'sienna',   #15
       ]

Phase.numPhase= 16
Phase.colorbar = cmap
Phase.shading= 'flat' 
Phase.aspectratio = 0.25
Phase.stream_line = streaml2
Phase.PTT_show = PTT_show1
Phase.labels = ['air','Astenosphere','Mantle_res1','Mantle_res2','Mantle_res3','Mantle_res4','Mantle_res5','HBasalt','NewHBasalt','Bas_res1','Bas_res2','Bas_res3','Intrusion','Felsic','Plume','New_Intrusion']
#                0      1              2              3            4             5              6             7          8          9          10         11          12         13      14       15
Phase.ticks   = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5]
#Phase.contour =0    1   2   3  4    5   6   7   8   9   10   11   12   13   14   15
buff=eval('Phase')
print('Properties to print phase')      
#########################################################################################
# PTT_Path routines
#########################################################################################

PTT_vis = PTT_Visualization ()
PTT_vis.lagrangian_grid = 'no'
PTT_vis.coloring        = ['field','mf']
PTT_vis.size_marker     =  0.5 
PTT_vis.criteria_ptt    = 'number_melting_event' # 'geometrical box', 'melt_fraction', 'number_melting_event' 
PTT_vis.value_criteria  = 1
PTT_vis.field_ptt       = ["Temperature [C]","Pressure [MPa]",'ID','Active','Phase']







