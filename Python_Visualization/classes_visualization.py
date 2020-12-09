######################################################################
# Classes Visualization
######################################################################
import numpy as np 
aspectratio=0.25


class contour():
    yes = 'no'
    field="temperature [C]"
    vector=[]
    color='k'
    text=0
    width= 0
    scalar = 0
        
class quiver():
        yes = 'no'
        nodes=0
        color='k'
        width= 0
        field="velocity [cm/yr]"
        cx = 0
        cy = 0
        cz = 0

class zoom():
        yes = 'no'
        zoom_int=0
        zoom_coord=(0,0.2)
        zoom_aspectr=(0.25)

class stream_line():
        yes = 'no'
        field = "velocity [cm/yr]"
        color = 'k'
        linewidth = ['scale',0.5]
        density   = 0.6
        cx = 0
        cy = 0
        cz = 0



     #   PTT_show = PTT_show()
        

class P_SURF():   
        name = 0
        clim = ()
        aspectratio = 0.25
        axis_lim = ['none', 0,0,0,0]
        n_component = 1
        component_show = ['mag','x','y','z']
        label       =  [r"v $cm\cdot yrs^(-1)$"]
 
class PTT_show():
        yes = 'no'
        PTT_Field_Sho = "num_melt_ev"
        timestep = 0 
        PTT_show_marker = '+'
        PTT_show_wake   = 10
        colormap   = 'jet'
        PTT_structure   = 'none'       
        wake  = 10 
        
class P_dym():   
        name = 0
        colorbar = 'Spectral_r'
        log = 'no'
        clim = ()
        shading= 'flat' 
        aspectratio = 0.25
        contour=contour()
        quiver = quiver()
        zoom   = zoom()
        stream_line = stream_line()
        axis_lim = ['none', 0,0,0,0]
        n_component = 1
        filter_air  = ['yes_show']
        component_show = ['mag','x','y','z']
        label       =  [r"v $cm\cdot yrs^(-1)$"]
        masked_value = ['no','cmin']
        PTT_show = PTT_show()

       
class P_phase():   
        name     = 'phase'
        numPhase = 0
        colorbar = 'Spectral_r'
        shading= 'flat' 
        aspectratio = 0.25
        zoom   = zoom()
        contour=contour()
        quiver=quiver()
        stream_line= stream_line()
        PTT_show = PTT_show()
        labels = ['air']
        ticks   = [0.5]
        

    
    
    
#Dynamic Properties [properties that are plotted name.pvd]
#####################################        
strain_rate_ii=P_dym()
strain_rate_ii.name='j2_strain_rate [1/s]'       
strain_rate_ii.colorbar='seismic'
strain_rate_ii.log= 'yes'
strain_rate_ii.clim =(-18,-12)
strain_rate_ii.shading='gouraud'
strain_rate_ii.aspectratio=aspectratio
strain_rate_ii.label   =  [r"$log_{10}(\dot{\eta}_{II}  s^(-1)$"]
########################################
Temperature=P_dym()
Temperature.name='temperature [C]'
Temperature.colorbar='Spectral_r'
Temperature.log= 'no'
Temperature.clim =(20,1800)
Temperature.shading='gouraud'
Temperature.aspectratio=aspectratio
Temperature.label   =  [r"$Temperature [C]$"]
############################################################
stress_ii=P_dym()
stress_ii.name="j2_dev_stress [MPa]"
stress_ii.colorbar='plasma'
stress_ii.log= 'yes'
stress_ii.clim =(-3,+3)
stress_ii.shading='gouraud'
stress_ii.aspectratio=aspectratio
stress_ii.label   =  [r"$log_{10}({\tau}_{II} [MPa]$"]
############################################################
Pressure=P_dym()
Pressure.name="pressure [MPa]"
Pressure.colorbar='Spectral_r'
Pressure.log= 'no'
Pressure.clim =(0,40000)
Pressure.shading='gouraud'
Pressure.aspectratio=aspectratio
Pressure.label   =  [r"$\bf{P} [MPa]$"]
############################################################
Pressure_Lith = P_dym()
Pressure_Lith.name = "litho_press [MPa]"
Pressure_Lith.colorbar='Spectral_r'
Pressure_Lith.log= 'no'
Pressure_Lith.clim =(0,40000)
Pressure_Lith.shading='gouraud'
Pressure_Lith.aspectratio=aspectratio
Pressure_Lith.label   =  [r"$P_{lithos} [MPa]$"]
############################################################
Pressure_Over = P_dym()
Pressure_Over.name = "litho_press [MPa]"
Pressure_Over.colorbar='Spectral_r'
Pressure_Over.log= 'no'
Pressure_Over.clim ='none'
Pressure_Over.shading='gouraud'
Pressure_Over.aspectratio=aspectratio
Pressure_Over.label   =  [r"$P_{lith}-\textbf{P} [MPa]$"]
############################################################
Velocity = P_dym()
Velocity.name =  "velocity [cm/yr]"
Velocity.colorbar='Spectral_r'
Velocity.log= 'no'
Velocity.clim =(-100,100)
Velocity.shading='gouraud'
Velocity.aspectratio=aspectratio
Velocity.n_component = 3
Velocity.component_show = ['mag','x','z']
Velocity.label =   [r"$Vel_{mag} cm \cdot yrs^{-1}$",r"$Vel_{x} cm \cdot yrs^{-1}$",r"$Vel_{z} cm \cdot yrs^{-1}$"] 
############################################################
Density = P_dym()
Density.name = "density [kg/m^3]"
Density.colorbar='Spectral_r'
Density.log= 'no'
Density.clim =(2700,4000)
Density.shading='gouraud'
Density.aspectratio=aspectratio
Density.label =   [r"$Density (\rho) [kg \cdot m^{-3}]$"]
############################################################
Visc_Creep = P_dym()
Visc_Creep.name = "visc_creep [Pa*s]"
Visc_Creep.colorbar='Spectral_r'
Visc_Creep.log= 'no'
Visc_Creep.clim =(18,24)
Visc_Creep.shading='gouraud'
Visc_Creep.aspectratio=aspectratio
Visc_Creep.label =   [r"$Viscosity_{visc} (\eta_{visc}) [Pa \cdot s]$"]
############################################################
Visc_Total = P_dym()
Visc_Total.name = "visc_creep [Pa*s]"
Visc_Total.colorbar='Spectral_r'
Visc_Total.log= 'no'
Visc_Total.clim =(18,24)
Visc_Total.shading='gouraud'
Visc_Total.aspectratio=aspectratio
Visc_Total.label =   [r"$Viscosity_{vep} (\eta_{vep}) [Pa \cdot s]$"]
############################################################
Displacement = P_dym()
Displacement.name = "tot_displ [km]"
Displacement.colorbar='Spectral_r'
Displacement.log= 'no'
Displacement.clim = 'none'
Displacement.shading='gouraud'
Displacement.aspectratio=aspectratio
Displacement.n_component = 3
Displacement.component_show = ['mag','x','z']
Displacement.label =   [r"$Displacement_{mag} (\textbf{u}) [km]$",r"$Displacement_{x} (\textbf{u}) [km]$",r"$Displacement_{z} (\textbf{u}) [km]$"]

#melt_fraction [ ]

Melt_Fraction = P_dym()
Melt_Fraction.name = "melt_fraction [ ]"
Melt_Fraction.colorbar='Blues'
Melt_Fraction.log= 'no'
Melt_Fraction.clim = (0.01,0.05)
Melt_Fraction.shading='gouraud'
Melt_Fraction.aspectratio=aspectratio
Melt_Fraction.label =   [r"$Melt Fraction ({Mf}) [n.d.]$"]
Melt_Fraction.masked_value = ['yes','cmin']


Melt_Extracted = P_dym()
Melt_Extracted.name = "Tot_Mext [ ]"
Melt_Extracted.colorbar='BuPu'
Melt_Extracted.log= 'no'
Melt_Extracted.clim = (0.01,0.5)
Melt_Extracted.shading='gouraud'
Melt_Extracted.aspectratio=aspectratio
Melt_Extracted.label =   [r"$Melt Extracted ({Mf_{Ext}}) [n.d.]$"]
Melt_Extracted.masked_value = ['yes','cmin']





##############################################################
#          Free Surf classes                                 #
##############################################################
Topography= P_SURF()
Amplitude = P_SURF()
#############################################################
# PT_Path classes                                           #
#############################################################
##
#Istantaneous PTT-grid
##


##
# Historical information 
##
class PTT_max_hist: 
        def __init__(self,number_marker):
            self.number_marker = number_marker
            self.ID            = np.zeros(self.number_marker,dtype=np.int)
            self.x             = np.zeros(self.number_marker,dtype=np.float)
            self.y             = np.zeros(self.number_marker,dtype=np.float)
            self.z             = np.zeros(self.number_marker,dtype=np.float)
            self.p             = np.zeros(self.number_marker,dtype=np.float)
            self.T             = np.zeros(self.number_marker,dtype=np.float)
            self.phase         = np.zeros(self.number_marker,dtype=np.int)
            self.Mf            = np.zeros(self.number_marker,dtype=np.float)
            self.Active        = np.zeros(self.number_marker,dtype=np.int)
            self.maxT          = np.zeros(self.number_marker,dtype=np.float)
            self.maxP          = np.zeros(self.number_marker,dtype=np.float)
            self.maxtimT       = np.zeros(self.number_marker,dtype=np.float)
            self.maxtimP       = np.zeros(self.number_marker,dtype=np.float)
            self.still_melt    = np.zeros(self.number_marker,dtype=np.int)
            self.num_melt_ev   = np.zeros(self.number_marker,dtype=np.int)
            self.activation_time = np.zeros(self.number_marker,dtype=np.float)


##
# Visualization routine commands
##

class PTT_Visualization ():
        lagrangian_grid = 'yes'
        number_of_nodes = (0,0,0)
        coloring        = ['monocolor','k']
        size_marker     =  0.5 
        criteria_ptt    = 'number_melting_event' # 'geometrical box', 'melt_fraction', 'number_melting_event' 
        geometrical_box = (0,0,0,0)
        value_criteria  = 0.0
        field_ptt       = ['Temperature [C]','Pressure [MPa]','Melt Fraction [ ]','ID','Active','Phase']
        


##
# Data Base for PTT  
##

class PTT_path : 
    def __init__(self,num_time, num_chosen):
        self.num_time = num_time
        self.num_chosen = num_chosen
        self.x_ptt      = np.zeros((self.num_time,self.num_chosen),dtype = np.float)
        self.y_ptt      = np.zeros((self.num_time,self.num_chosen),dtype = np.float)
        self.z_ptt      = np.zeros((self.num_time,self.num_chosen),dtype = np.float)
        self.p_ptt      = np.zeros((self.num_time,self.num_chosen),dtype = np.float)
        self.T_ptt      = np.zeros((self.num_time,self.num_chosen),dtype = np.float)
        self.mf_ptt     = np.zeros((self.num_time,self.num_chosen),dtype = np.float)
        self.active     = np.zeros((self.num_time,self.num_chosen),dtype = np.int)
        self.num_melt_ev = np.zeros((self.num_time,self.num_chosen),dtype=np.int)
        self.no_melt_ev  = np.zeros((self.num_time,self.num_chosen),dtype=np.int)
        
        self.time_ptt   = np.zeros((self.num_time),dtype = np.float)



