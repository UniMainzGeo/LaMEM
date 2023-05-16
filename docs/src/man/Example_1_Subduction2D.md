# 4.1 2D subduction setup with nonlinear rheologies

This is a worked-out example of how to create a 2D subduction setup with temperature dependent rheology.

The example input file can be found under:
```
/input_models/SubductionWithParticles/Subduction2D_FreeSurface_MATLABParticles_Nonlinear_DirectSolver.dat
```
whereas the MATLAB/Octave file is here:
```
/input_models/SubductionWithParticles/CreateMarkers_Subduction_Tdependent_FreeSurface_parallel.m
```

### 4.1.1 Input model geometry 
We start with having a look at the input file. With Visual Studio Code this can be done with
```
$ code Subduction2D_FreeSurface_MATLABParticles_Nonlinear_DirectSolver.dat
```

***Units***
At the beginning of the script, we specify the units
```
#===============================================================================
# Scaling
#===============================================================================

	units            = geo		# geological units 
	
	unit_temperature = 1000
	unit_length      = 1e3
	unit_viscosity   = 1e20
	unit_stress      = 1e9
``` 
Internally, LaMEM computes everything in non-dimensional units in order to preent round-off (numerical) errors. In order to transfer things from dimensional to non-dimensional units you need to specify these values.
LaMEM has 3 modes for units
```
units = none
units = SI
units = geo
```
as you can imagine, the first implies that all input is in non-dimensional units, the second that all is given in SI units. The third, and the one that is most typically used, employs geodynamically sensible units, such as temperature in celcius, length scales in kilometers, times in million years, etc.

The units that you define here are in SI-units, which implies that the typical lengthscale is 1000 meters, and the typical stress is $10^9$ Pa.

***Model domain***
The model domain and numerical reoslution are specified here:
```
#===============================================================================
# Grid & discretization parameters
#===============================================================================

# Number of cells for all segments
	nel_x = 256
	nel_y = 2
	nel_z = 64

# Coordinates of all segments (including start and end points)

	coord_x = -1500 1500
	coord_y = -10   10
	coord_z = -660  40
```

Note that LaMEM is a fully 3D code; there is no pure-2D part. Yet you can nevertheless do 2D simulations by selecting only 1 or 2 cells in the y-direction as done above. If you select multigrid as a solver, you need to have at least 2 cells here. Note that for 2D, the 'thin' direction should always be the y-direction.

The remaining part of this block defines the number of elements in each direction and the left/right, front/back and bottom/top coordinates (in km, as we use geo units).

Also note that you can always specify parameters on the command-line, that will *overrule* whatever is written in the input file. 
For example: 
```
../../bin/LaMEM -ParamFile Subduction2D_FreeSurface_MATLABParticles_Nonlinear_DirectSolver.dat -nel_x 64 -nel_z 16
```
will perform a simulation with 64 by 2 by 16 elements instead.




