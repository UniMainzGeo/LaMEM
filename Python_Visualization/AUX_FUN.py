#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scipy.io as sio
import sys,os,fnmatch
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib
from matplotlib import cm, pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import shutil
import time
import glob
import re
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from pylab import figure, axes, pie, title, show
import pylab
#from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
#import json
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.ma import masked_array
import matplotlib.cbook as cbook
from initial_input import aspectratio
from  classes_visualization import * 
import copy

"""
AUX_FUNCTION
Created on Wed Mar 25 14:09:07 2020

@author: piccolo

# This is a simply dictionary that allow python to parse the information of the 
# pvd file. In short, the program read the file, and select automatically the path
# to the folder that you require. Do not touch it.
# Define a line search parser
"""

def _parse_line(line,key):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex

    """
    rx = rx_dict.get(key)
    match = rx.search(line)
    if match:
        return match
    return None

rx_dict = {
    'timestep': re.compile(r'timestep=(\".*\") '),
    'file': re.compile(r'file=(\".*\")'),
    'WholeExtent': re.compile(r'WholeExtent=(\".*\")')
    }

   # File list 
   
def _file_list(fname):
   
   
    Time=[]    # Time of the numerical experiment
    Flist=[]   # File list 
   
    F=open(fname,'r')
    d=F.readlines()
    F.seek(0)
    for line in d:
        key = 'timestep'
        match = _parse_line(line,key)
        if match:
            Time.append(match.group(1))
            match =[] 
        key = 'file'
        match = _parse_line(line,key)
        if match:
            Flist.append(match.group(1))
    F.close()
    time=np.zeros(len(Time))
    n_tim_step=[]
    i=0
    for itime in Time:
        tm=itime[1:-1]
        time[i]=np.float(tm)
        i+=1

    for istep in Flist:
        index = [x.start() for x in re.finditer('_', istep)]
        istp=istep[index[0]+1:index[1]]
        n_tim_step.append(istp)
        
        
        
    return time, Flist, n_tim_step

def _parse_grid(Filename,p):
    F=open(Filename,'r')
    d=F.readlines()
    F.seek(0)
    for line in d:
        key = 'WholeExtent'
        match = _parse_line(line,key)
        if match:
            nodes_num=(match.group(1))
            match=[]
            break
           
    F.close()
    nodes_num=nodes_num[1:-1]
    if p == 1:
        nodes_num = nodes_num[:-15]


    nodes_num=[int(s) for s in nodes_num.split(' ')]
    nx = nodes_num[1]
    ny = nodes_num[3]
    nz = nodes_num[5]
    return nx, ny, nz

def _figure_maker(ptsave,buf,x,z,iprop,t_cur,tstep,Pic,label,component):
     
    time_sim = "{:.3f}".format(t_cur)
    tick='Time = %s Myrs' %time_sim

    if(component =='none'):
        ptsave_b=os.path.join(ptsave,iprop)
    else:
        ptsave_b=os.path.join(ptsave,iprop)

        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)
        ptsave_b=os.path.join(ptsave,iprop,component)
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)

    buf=buf[:,0,:]
    if not os.path.isdir(ptsave_b):
        os.mkdir(ptsave_b)
    
    fna='Fig'+str(tstep)+'.png'
    fn=os.path.join(ptsave_b,fna)

    if(Pic.axis_lim[0]=='none'):
        xmin = np.min(x)
        xmax = np.max(x)
        zmin = np.min(z)
        zmax = np.max(z)
    else:
        xmin = Pic.axis_lim[1]
        xmax = Pic.axis_lim[2]
        zmin = Pic.axis_lim[3]
        zmax = Pic.axis_lim[4]
    fg = figure()
    ax0 = fg.gca()
    cmap = plt.cm.get_cmap(Pic.colorbar)
    if(Pic.masked_value[0]=='yes'):
        if(Pic.masked_value[1]=='cmin'):
            cmap.set_under(color='k', alpha=0.0)       
    elif(Pic.masked_value[1]=='cmax'):
            cmap.set_over(color='k', alpha=0.0)       

    elif(Pic.masked_value[1]=='cmax & cmin'):
            cmap.set_under(color='k', alpha=0.0)
            cmap.set_over(color='k', alpha=0.0)       

 
            
    cf = ax0.pcolormesh(x,z,(buf),cmap=cmap,shading=Pic.shading,vmin=Pic.clim[0],vmax=Pic.clim[1])
    fg.colorbar(cf, ax=ax0,orientation='horizontal', label = Pic.name)

    ####################################################
    if (Pic.contour.yes == 'yes') | (Pic.quiver.yes == 'yes') | (Pic.stream_line.yes == 'yes'): 
        if (Pic.contour.yes == 'yes'):
            field=Pic.contour.Scalar
            level=Pic.contour.vector
            CS = ax0.contour(x, z, field, levels = level, colors=Pic.contour.color,linewidths =Pic.contour.width)
            ax0.clabel(CS, fontsize=5, inline=1)
            
        
        if Pic.quiver.yes == 'yes':
            stp = Pic.quiver.nodes
            ax0.quiver(x[0::stp,0::stp], z[0::stp,0::stp], Pic.quiver.cx[0::stp,0::stp], Pic.quiver.cz[0::stp,0::stp])

        if Pic.stream_line.yes == 'yes':
            if(Pic.stream_line.linewidth[0] == 'scale'):
                speed = np.sqrt(Pic.stream_line.cx**2+Pic.stream_line.cz**2)
                lw = Pic.stream_line.linewidth[1]*speed /np.max(speed)
            else:
                lw = Pic.stream_line[1]
            size_x = np.int((np.size(x[0,0:]))) 
            size_z = np.int((np.size(z[0:,0])))
            XG = np.linspace(np.min(x),np.max(x),size_x)
            ZG = np.linspace(np.min(z),np.max(z),size_z)
            ax0.streamplot(XG, ZG,Pic.stream_line.cx, Pic.stream_line.cz,  color=Pic.stream_line.color, linewidth=lw)

           
        if(Pic.PTT_show.yes == 'yes'):
            ax1 = fg.gca()
            fb= ax1.scatter(Pic.PTT_show.PTT_structure.x,Pic.PTT_show.PTT_structure.z,marker = 'P',c = (Pic.PTT_show.PTT_structure.T/(Pic.PTT_show.PTT_structure.p/1000)),cmap = 'YlOrRd' ,s=0.1,)

            fg.colorbar(fb, ax=ax0,orientation='vertical', label = r'$T/P [C/GPa]$')

    ax0.set_aspect((float(xmax-xmin)/float(zmax-zmin)) *Pic.aspectratio)
    ax0.set_title(tick)
    ax0.set_ylim(zmin,zmax)
    ax0.set_xlim(xmin,xmax)
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    plt.draw()    # necessary to render figure before saving
    fg.savefig(fn,dpi=300,transparent=True)
    if(Pic.zoom.yes == 'yes'):
        ax0.set_ylim(Pic.zoom.zoom_coord[2],Pic.zoom.zoom_coord[3])
        ax0.set_xlim(Pic.zoom.zoom_coord[0],Pic.zoom.zoom_coord[1])
        ptsave_b=os.path.join(ptsave,iprop,'zoom')
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)
        fna = 'Zoom_Fig'+str(tstep)+'.png'
        fn=os.path.join(ptsave_b,fna)
        fg.savefig(fna,dpi=300,transparent=False)


        
    fg.clear()
    plt.close()
    
    
    
def _figure_maker_phase(ptsave,phase,x,z,t_cur,tstep,Pic,xd,zd):
    
    ptsave_b=os.path.join(ptsave,'phase')
    time_sim = "{:.3f}".format(t_cur)
    tick='Time = %s Myrs' %time_sim
    
    phase=phase[:,0,:]
    if not os.path.isdir(ptsave_b):
        os.mkdir(ptsave_b)
    
    fna='Fig'+str(tstep)+'.png'
    fn=os.path.join(ptsave_b,fna)
    fg = figure()
    ax0 = fg.gca()
    cmap= colors.ListedColormap(Pic.colorbar)
    cf=ax0.pcolormesh(x, z, phase,cmap=cmap,vmin=0,vmax=Pic.numPhase, shading='gouraud')
    cbar = fg.colorbar(cf, ax=ax0,orientation='horizontal',ticks = Pic.ticks)
    cbar.ax.set_xticklabels(Pic.labels, size = 3)  # horizontal colorbar
    ############################################################
    if (Pic.contour.yes == 'yes') | (Pic.quiver.yes == 'yes'): 
        if (Pic.contour.yes == 'yes'):
            field=Pic.contour.Scalar
            level=Pic.contour.vector
            CS = ax0.contour(xd, zd, field, levels = level, colors=Pic.contour.color,linewidths =Pic.contour.width)
            ax0.clabel(CS, fontsize=5, inline=1)

        
        if Pic.quiver.yes == 'yes':
            stp = Pic.quiver.nodes
            ax0.quiver(xd[0::stp,0::stp], zd[0::stp,0::stp], Pic.quiver.cx[0::stp,0::stp], Pic.quiver.cz[0::stp,0::stp],colors=Pic.quiver.color)

    if Pic.stream_line.yes == 'yes':
        if(Pic.stream_line.linewidth[0] == 'scale'):
            speed = np.sqrt(Pic.stream_line.cx**2+Pic.stream_line.cz**2)
            lw = Pic.stream_line.linewidth[1]*speed /np.max(speed)
        else:
            lw = Pic.stream_line.linewidth[1]
            
        size_x = np.int((np.size(xd[0,0:]))) 
        size_z = np.int((np.size(zd[0:,0])))
        XG = np.linspace(np.min(xd),np.max(xd),size_x)
        ZG = np.linspace(np.min(zd),np.max(zd),size_z)
           
            
        ax0.streamplot(XG, ZG,Pic.stream_line.cx, Pic.stream_line.cz,  color=Pic.stream_line.color, linewidth=lw)
     
            
            
           
    if(Pic.PTT_show.yes == 'yes'):
        ax1 = fg.gca()
        ax1.scatter(Pic.PTT_show.PTT_structure.x,Pic.PTT_show.PTT_structure.z,marker = 'P',c = 'black',s=0.1,alpha = 0.1)
    
    ############################################################### 
    ax0.set_aspect((np.max(x)-np.min(x)) / float(np.max(z)-np.min(z)) *Pic.aspectratio)
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    ax0.set_title(tick)
    ax0.set_ylim(np.min(z),np.max(z))
    ax0.set_xlim(np.min(x),np.max(x))
    ###############################################################    
    plt.draw()    # necessary to render figure before saving
    fg.savefig(fn,dpi=300,transparent=False)
    if(Pic.zoom.yes == 'yes'):
        ax0.set_ylim(Pic.zoom.zoom_coord[2],Pic.zoom.zoom_coord[3])
        ax0.set_xlim(Pic.zoom.zoom_coord[0],Pic.zoom.zoom_coord[1])
        fna = 'Zoom_Fig'+str(tstep)+'.png'
        fn=os.path.join(ptsave_b,fna)
        fg.savefig(fn,dpi=300,transparent=False)
    fg.clear()
    plt.close()
    
