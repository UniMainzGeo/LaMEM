#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 15:12:41 2020

@author: piccolo
"""


import sys,os,fnmatch
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import shutil
import time as perf
import glob
import re
from numpy.ma import masked_array
from  classes_visualization import * 
from vtk.numpy_interface import dataset_adapter as dsa



"""
Import useful function 
"""
from AUX_FUN import *
from AUX_FUN import _parse_line
from AUX_FUN import _file_list
from AUX_FUN import _figure_maker
from AUX_FUN import _figure_maker_phase
from AUX_FUN import _parse_grid


from initial_input import * 

reader=vtk.vtkXMLGenericDataObjectReader()
VTK_SET=reader.SetFileName
VTK_UPDATE=reader.Update
VTK_RELEASE=reader.UnRegister
VTK_POINTS=vtk.vtkPoints
VTK_OUT=reader.GetOutput
AllPoints=VTK_POINTS()
GetData=AllPoints.GetData
nodes_vtk_array= GetData()

VTK_SET_dyn=reader.SetFileName
VTK_UPDATE_dyn=reader.Update
VTK_RELEASE_dyn=reader.UnRegister
VTK_POINTS_dyn=vtk.vtkPoints
VTK_OUT_dyn=reader.GetOutput
AllPoints_dyn=VTK_POINTS()
GetData_dyn=AllPoints_dyn.GetData
nodes_vtk_array_dyn= GetData_dyn()

reader_ptt = vtk.vtkXMLGenericDataObjectReader()
VTK_SET_ptt=reader_ptt.SetFileName
VTK_UPDATE_ptt=reader_ptt.Update
VTK_RELEASE_ptt=reader_ptt.UnRegister
VTK_POINTS_ptt=vtk.vtkPoints
VTK_OUT_ptt=reader_ptt.GetOutput
AllPoints_ptt=VTK_POINTS()
GetData_ptt=AllPoints_ptt.GetData
nodes_vtk_array_ptt= GetData_ptt()
    


# If you follow the suggestion, change the name here.
# Do not touch it 
fname2 = fname+'.pvd'
fname=os.path.join(Folder,Test_Name,fname2)

time, Flist, n_tim_step =_file_list(fname)    # Retrive the file list and the time step information

max_time_step = 0
if ((PTT == 1) & (PTT_vis.criteria_ptt == 'number_melting_event')):
    for istp in Flist:
        fn=istp[1:istp.find('/')]
        t_cur=time[max_time_step]
        tstep=n_tim_step[max_time_step]
        print(tstep)
  
        Filename_ptt=os.path.join(Folder,Test_Name,fn,passive_tracers)
        VTK_SET_ptt(Filename_ptt)
        VTK_UPDATE_ptt()
        VTK_OUT_ptt().GetPoints()
        number_marker = VTK_OUT_ptt().GetNumberOfPoints()
        PTT_Hist = PTT_max_hist(number_marker)
        for i in range(number_marker):
            
            
            PTT_Hist.x[i],PTT_Hist.y[i],PTT_Hist.z[i] = VTK_OUT_ptt().GetPoint(i)
        
        ######################################################################
        # Retrieve ptt vector
        ######################################################################
        for name in PTT_vis.field_ptt:
            buf_ptt=vtk_to_numpy(VTK_OUT_ptt().GetPointData().GetArray(name))
            if name == 'Temperature [C]':
                PTT_Hist.T=buf_ptt
            if name == 'Pressure [MPa]' :
                PTT_Hist.p = buf_ptt
            if name == 'Mf [ ]'          :
                PTT_Hist.Mf = buf_ptt
            if name == 'Active'         :
                PTT_Hist.Active = buf_ptt
            if name == 'ID'             : 
                PTT_Hist.ID =buf_ptt
            if name == 'Phase'          :
                PTT_Hist.phase = buf_ptt
                
            #################################################################
            # Compute the relevant historical information 
            ################################################################
        if(max_time_step >=102):
            print('ahhaha')
        for n_m in range(number_marker):
            if((PTT_Hist.x[n_m] == -np.inf) | (PTT_Hist.z[n_m] == -np.inf ) | (PTT_Hist.z[n_m] == -np.inf ) ):
                PTT_Hist.Active[n_m] == 0
                PTT_Hist.x           == 0.0
                PTT_Hist.y           == 0.0
                PTT_Hist.z           == 0.0
                PTT_Hist.Mf          == 0.0
                PTT_Hist.p           == 0.0
                PTT_Hist.T           == 0.0
                PTT_Hist.phase       == 0 
                PTT_Hist.maxP        == 0.0
                PTT_Hist.num_melt_ev == 0
                PTT_Hist.activaction_time = -2.0 
                print('caught a rougue')
            
            if((PTT_Hist.Active[n_m] == 1) & (PTT_Hist.activation_time[n_m] == 0.0)):
                PTT_Hist.activation_time[n_m] = t_cur 
            if(PTT_Hist.T[n_m]>PTT_Hist.maxT[n_m]):
                PTT_Hist.maxT[n_m] = PTT_Hist.T[n_m]
                PTT_Hist.maxtimT[n_m] = t_cur
            if(PTT_Hist.p[n_m]>PTT_Hist.maxP[n_m]):
                PTT_Hist.maxP[n_m] = PTT_Hist.p[n_m]
                PTT_Hist.maxtimP[n_m] = t_cur
            if((PTT_Hist.Mf[n_m] > 0.0) & (PTT_Hist.still_melt[n_m] == 0)):
                PTT_Hist.num_melt_ev[n_m]+=1
            if(PTT_Hist.Mf[n_m] > 0.0):
                PTT_Hist.still_melt[n_m] = 1
            else:
                PTT_Hist.still_melt[n_m] = 0
                
        
        ##############################################################                
        # Add Field to vtu file
        ##############################################################


# Add data set and write VTK file
        old_data = VTK_OUT_ptt()
        New = dsa.WrapDataObject(old_data)
        New.PointData.append(PTT_Hist.num_melt_ev, "num_melt_ev")
        New.PointData.append(PTT_Hist.activation_time, "activation_time [Myrs]")
        New.PointData.append(PTT_Hist.maxP, "maxP [MPa]")
        
        Filename_ptt=os.path.join(Folder,Test_Name,fn,'update_passive.vtk')
        
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(Filename_ptt)
        writer.SetInputData(New.VTKObject)
        writer.Write()

    
                
                
                
                
                
        max_time_step +=1 
                    