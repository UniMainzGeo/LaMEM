# 3. Initial model setup

There are several ways in which you can define the initial geometry of your model setup:

[3.1 Build-in geometries](BuildinGeometries.md)

[3.2 Matlab/Octave geometries](GenerateModelSetup_MATLAB.md) 

[3.3 Geometry from geomIO-file](PolygonGeometry.md)

The build-in geometries and Matlab/Octave files allow you to specify the initial temperature structure. 
In addition, you can also independently set a temperarature profile with: 

[3.4 Temperature input file](TempInputFile.md)

In some cases you might be interested to start with a different temperature profile, such as a steady-state temperature profile for your setup. 
See:

[3.5 Temperature diffusion](TempDiffusion.md)

By default, the internal free surface is set to a constant value. 
In some cases that is insufficient and you can also set the initial topography with a file: 

[3.6 Topography input file](TopoFile.md)