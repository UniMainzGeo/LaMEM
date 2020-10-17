# Simple Python code that retrive P-T conditions of continental crust generation per each timestep 

import sys,os,fnmatch
import numpy as np
import matplotlib.pyplot as plt 
plt.switch_backend('agg')
import matplotlib
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import time
import glob
from matplotlib import rcParams
import time
import re



def Plot_Analytics_vs_Numerics(dir):

    ######################################################################################
    # Some useful function
    ######################################################################################
    reader=vtk.vtkXMLGenericDataObjectReader()
    VTK_SET=reader.SetFileName
    VTK_UPDATE=reader.Update
    VTK_POINTS=vtk.vtkPoints
    VTK_OUT=reader.GetOutput
    AllPoints=VTK_POINTS()
    GetData=AllPoints.GetData
    nodes_vtk_array= GetData()
    ROOT=np.sqrt
    QUADR=np.square
    OPEN=os.path.join
    perf=time.time
    ######################################################################################


    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    # Define a line search parser
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
        }




    Time=[]
    Flist=[]


    fname=("t13.pvd")


    F=open(fname,'r')
    d=F.readlines()
    F.seek(0)
    # Create a parser. Read the pvd file and generate a list with all the path to 
    #the timestep folder.
    
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
    # Set the dimension of the domain     
    dim1=129
    dim2=2
    dim3=5

    i=0
    #Initialize the variable to append 
    time1=[]
    max_lamem=[]
    max_anal=[]
    appendtime=time1.append
    append_ML=max_lamem.append
    append_MA=max_anal.append

    # Loop over the list of timestep

    for itime in Flist:
    
        start=perf()
        #####################################################
        itime1=itime[1:-1]      #Folder
        Filename=OPEN(itime1)   #Open the folder
        tm=Time[i]              #time
        buf=float(tm[1:-1])     #strtonumber
        print(itime)
        #####################################################
        s1=perf()
        VTK_SET(Filename)
        s2=perf()
        print(1,s2-s1)
        VTK_UPDATE()
        s3=perf()
        print(2,s3-s2)
        VTK_OUT().GetPoints(AllPoints)
        s4=perf()
        print(3,s4-s3)        
        end_init=perf()
        print('Initial stage',end_init-start)
        # Create the vtk Array
        start_array=perf()
        #Get the coordinates of the nodes
        nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
        x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:
        ,1] , nodes_nummpy_array[:,2]
        #Get Scalar Arrays
        temperature_vtk_array = reader.GetOutput().GetPointData().GetArray('temperature [C]')
        #Transfer to real variable and delete the other variable
        temp=vtk_to_numpy(temperature_vtk_array)
        end_array=perf()
        print('Create and delete Array took', end_array-start_array)
        
        Z_GRID=z.reshape(dim1,dim2,dim3)
        Z=Z_GRID[:,1,1]
        T_lamem=temp.reshape(dim1,dim2,dim3)
        T_lamem=T_lamem[:,1,1]
    ############################################################################
    # Analytical solution
    ############################################################################
    #T(x,t)=Tmax/()[...]  
        rho=3000                     #[kg/m^3]
        sigma=50*1e3                 #[m]
        k=3                          #[W/m2]
        Cp=1050                      #[J/kg]
        TMax=400                     #[C]
        kappa= float(k)/(Cp*rho)     #[m2/s]
        t=buf*365.25*24*60*60*1e6    #[s]
    
        T_anal  = float(TMax)/np.sqrt(1 + float(4*t*kappa)/sigma**2)*np.exp(-(Z*1e3)**2/(sigma**2 + 4*t*kappa)) #Taken from Anton's Class 

        ###########################################################################
        # Figure (1) Comparison between analytical solution and numerical one
        ##########################################################################
        filename=dir+'Temperature'+str(i)
        plt.figure(1)   
        plt.clf()
        plt.plot(Z, T_lamem,'r',Z, T_anal,'k--')
        plt.title('T at time '+ str(buf) + ' Myrs')
        plt.xlabel('Z,[km]')
        plt.ylabel('T,[C]')
        plt.axis((-50, 50, 0.0, TMax))
        plt.show()
        plt.pause(1e-1)
        plt.savefig(filename,dpi=600)

        ##########################################################################
        # Figure (2) Error between analytical solution and numerical one
        ##########################################################################

        err=T_anal-T_lamem

        filename= dir+'Error'+str(i) 
        plt.figure(2)   
        plt.clf()
        plt.plot(Z, err,'r--')
        plt.title('T at time '+ str(buf) + ' Myrs')
        plt.xlabel('Z,[km]')
        plt.ylabel('err T,[C]')
        plt.xlim((-50, 50))
        plt.show()
        plt.pause(1e-1)   
        plt.savefig(filename,dpi=600)


        # Is the error symmetric (check if the code is symmetric). If it is not, check
        # if there is random noises in your setup
        err2=err-np.flip(err)
        print('average err-flip(err) is',"{:.{}f}".format( np.mean(err2), 5 ), '[C]')
        
        #Append the maximum temperature 
        appendtime(float(t)/(365.25*24*60*60*1e6))
        append_ML(np.max(T_lamem))
        append_MA(np.max(T_anal))
        
        i+=1
        end=perf()
        print(end-start)

    max_anal=np.array(max_anal)

    max_lamem=np.array(max_lamem)
    filename= 'MaxError'
    plt.figure(2)   
    plt.clf()
    plt.plot(time1, max_anal-max_lamem,'r--')
    plt.title('Difference between LaMEM maximum temperature and analytical maximum')
    plt.xlabel('time,[Myrs]')
    plt.ylabel('err T,[C]')
    plt.show()
    plt.pause(1e-1)   
    plt.savefig(filename,dpi=600)
        

# for debugging:
#Plot_Analytics_vs_Numerics('./t12_Temperature_diffusion/')
