# 3. Initial model setup

There are several ways in which you can define the initial geometry of your model setup. Click on the blue titles to get more extended info:

[3.1 Build-in geometries](BuildinGeometries.md)
This is for simpler geometries and can directly be set from the LaMEM `*.dat` file. The advantage is that it will work for any resolution, so if this is possible in your case go for it!

3.2 Julia based setup
More complicated setups can be constructed in julia by using the [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) package. This is the recommended way to do this.

[3.3 Matlab/Octave geometries](GenerateModelSetup_MATLAB.md) 
Previously, we used matlab scripts to create model setups. We keep this here for historical reasons (and because the geomIO workflow has not yet been fully ported to julia). 

[3.4 Geometry from geomIO-file](PolygonGeometry.md)

The build-in geometries and Matlab/Octave files allows you to specify the initial temperature setup as well. 
In addition, you can also, independently, set a temperarature profile with: 

[3.5 Temperature input file](TempInputFile.md)
The initial thermal structure can be set in julia or the matlab scripts (as described above) but also seperately if you want, as described here. 

[3.6 Temperature diffusion](TempDiffusion.md)
In some cases you might be interested to start with a different temperature profile, such as a steady-state temperature profile for your setup as described here.

[3.7 Topography input file](TopoFile.md)
By default, the internal free surface is set to a constant elevation. Instead, you can also set the initial topography with a file, as explained here.