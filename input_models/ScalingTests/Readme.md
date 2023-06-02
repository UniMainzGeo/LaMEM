# Scalability Tests

This summarizes scalability tests of LaMEM on different HPC machines.

In general, we consider 3 different setups (or miniapps):

- *ThermalDiffusion* - pure diffusion, no mechanics involved. Tests the effect of the diffusion equation itself. Whereas thermal diffusion is only part of the computational time of a typical simulation, it has the advantage that the code is much easier to port to GPUs, and will therefore  
- *FallingSpheres_Multigrid.dat* - variable (linear) viscous Stokes solver with 10 high-density spheres. In our experience this is a good test that represents more complicated viscous structures well. 
- *to be determined* - test case that includes nonlinear viscelastoplastic rheologies and a (simplified) topography, for example a landslide. This is still to be generated
  

### LUMI - StokesFallingSpheres

All cases are performed on LUMI-C with a viscosity contrast of 1000. As a coarse grid solver we employ SuperLU_DIST & choose the levels such that the coarse grid level remains the same in size

```
mpiexec -n 64 ./LaMEM -ParamFile FallingSpheres_Multigrid.dat -gmg_pc_mg_levels 5 -nel_x 128 -nel_y 128 -nel_z 128 -log_view
```

| resolution    | cores  | nodes | # MG levels | # iter | total time [s] | MG Apply [s] | time coarse [s] | logfile
|---------------|--------|-------| ------------|------- |----------------|--------------| --------------- | ------
| 32x32x32      |  1     |  1    |      3      |   37   |      7.6       |     4.2      |                 | slurm-3628882.out
| 64x64x64      |  8     |  1    |      4      |   58   |      43.0      |     38.9     |                 | slurm-3628768.out
| 128x128x128   |  64    |  1    |      5      |   71   |      97.7      |     89.0     |                 | slurm-3628803.out
| 128x128x128   |  128   |  1    |      5      |   70   |      44.0      |     39.5     |                 | slurm-3628709.out
| 128x128x256   |  128   |  1    |      5      |   80   |      109.0     |    100.0     |                 | slurm-3628755.out
| 256x256x256   |  512   |  4    |      6      |   75   |      110.0     |     99.5     |                 | slurm-3628873.out
| 512x512x512   |  512   |  4    |      7      |   ??   |      ??        |     ??       |                 | slurm-3628971.out



### MOGON2 - StokesFallingSpheres

These tests are performed on the skylake partition of MOGON2 (Mainz)

| resolution    | cores  | nodes | # MG levels | # iter | total time [s] | MG Apply [s] | time coarse [s] | logfile
|---------------|--------|-------| ------------|------- |----------------|--------------| --------------- | ------
| 128x128x128   |  64    |  2    |      3      |   -    |      -         |     --       |       OuM       | slurm-13547787.out
| 128x128x128   |  64    |  2    |      4      |   73   |      78.9      |     70.5     |                 | slurm-13547775.out
| 128x128x128   |  64    |  2    |      5      |   71   |      75.5      |     67.1     |                 | slurm-13547723.out
| 256x256x256   |  512   |  16   |      6      |   75   |      81.8      |     72.4     |                 | slurm-13548462.out
| 256x256x256   |  512   |  16   |      5      |   78   |      87.3      |     77.6     |                 | slurm-13548495.out
| 256x256x512   |  1024  |  16   |      7      |   --   |      ----      |     ----     |                 | slurm-13548914.out
| 512x512x512   |  4096  |  16   |      7      |   --   |      ----      |     ----     |                 | slurm-13548764.out
