# Scalability Tests

This summarizes scalability tests of LaMEM on different HPC machines.

In general, we consider 3 different setups (or miniapps):

- *Diffusion.dat* - pure diffusion, no mechanics involved. Tests the effect of the diffusion equation itself. Whereas thermal diffusion is only part of the computational time of a typical LaMEM simulation, it has the advantage that the code is much easier to port to GPUs, and will therefore serve as a good first benchmark.
- *FallingSpheres.dat* - variable (linear) viscous Stokes solver with 10 high-density spheres. In our experience this is a good test that represents more complicated viscous structures well. 
- *Volcano.dat* - test case that includes viscoelastoplastic rheologies and a topography, to simulate hangslope instabilities on volcanoes.
  
### LUMI - StokesFallingSpheres

All cases are performed on LUMI-C with a viscosity contrast of 1000. As a coarse grid solver we employ SuperLU_DIST & choose the levels such that the coarse grid level remains the same in size

```
mpiexec -n 64 ./LaMEM -ParamFile FallingSpheres.dat -gmg_pc_mg_levels 5 -nel_x 128 -nel_y 128 -nel_z 128 -log_view
```

| resolution    | cores  | nodes | # MG levels | # iter | total time [s] | MG Apply [s] | time coarse [s] | logfile
|---------------|--------|-------| ------------|------- |----------------|--------------| --------------- | ------
| 32x32x32      |  1     |  1    |      3      |   37   |      7.6       |     4.2      |                 | slurm-3628882.out
| 64x64x64      |  8     |  1    |      4      |   58   |      43.0      |     38.9     |                 | slurm-3628768.out
| 128x128x128   |  64    |  1    |      5      |   71   |      97.7      |     89.0     |                 | slurm-3628803.out
| 128x128x128   |  128   |  1    |      5      |   70   |      44.0      |     39.5     |                 | slurm-3628709.out
| 128x128x256   |  128   |  1    |      5      |   80   |      109.0     |    100.0     |                 | slurm-3628755.out
| 256x256x256   |  512   |  4    |      6      |   75   |      110.0     |     99.5     |                 | slurm-3628873.out
| 256x256x256   |  512   |  4    |      5      |   78   |      110.3     |     99.6     |                 | slurm-3630017.out



*GAMG* coarse grid solver of 32^3 solved on 64 cores with PC TELESCOPE:
example command on 512 cores:
```
srun ./LaMEM -ParamFile FallingSpheres.dat -gmg_pc_mg_levels 4 -nel_x 256 -nel_y 256 -nel_z 256 -log_view -crs_pc_type telescope -crs_pc_telescope_reduction_factor 8 -crs_telescope_pc_type gamg  -crs_pc_view

```

| resolution    | cores  | nodes | # MG levels | Telescope reduction | # iter | total time [s] | MG Apply [s] | time coarse [s] | logfile
|---------------|--------|-------| ------------| --------------------|------- |----------------|--------------| --------------- | ------
| 256x256x256   |  512   |  4    |      4      |   8                 | 124    |     208        |    123.1     |    1.9          | slurm-3631942.out
| 512x512x512   |  512   |  32   |      6      |   8                 | --     |     --         |    --    |    OoM          | slurm-3631954.out - 
| 512x512x512   |  512   |  43   |      6      |   8                 | --     |     --         |    --    |    --          | slurm-3635671.out


### MOGON2 - StokesFallingSpheres

These tests are performed on the skylake partition of MOGON2 (Mainz)

| resolution    | cores  | nodes | # MG levels | # iter | total time [s] | MG Apply [s] | time coarse [s] | logfile
|---------------|--------|-------| ------------|------- |----------------|--------------| --------------- | ------
| 128x128x128   |  64    |  2    |      3      |   -    |      -         |     --       |       OuM       | slurm-13547787.out
| 128x128x128   |  64    |  2    |      4      |   73   |      78.9      |     70.5     |                 | slurm-13547775.out
| 128x128x128   |  64    |  2    |      5      |   71   |      75.5      |     67.1     |                 | slurm-13547723.out
| 256x256x256   |  512   |  16   |      6      |   75   |      81.8      |     72.4     |                 | slurm-13548462.out
| 256x256x256   |  512   |  16   |      5      |   78   |      87.3      |     77.6     |                 | slurm-13548495.out
| 256x256x512   |  1024  |  16   |      7      |   85   |      93.8      |     82.9     |                 | slurm-13548914.out
| 512x512x512   |  4096  |  16   |      7      |   --   |      ----      |     ----     |                 | slurm-13548764.out



*GAMG* coarse grid solver of 32^3 solved on 64 cores with PC TELESCOPE:

| resolution    | cores  | nodes | # MG levels | Telescope reduction | # iter | total time [s] | MG Apply [s] | time coarse [s] | logfile
|---------------|--------|-------| ------------| --------------------|------- |----------------|--------------| --------------- | ------
| 128x128x128   |  64    |  2    |      3      |   1                 | 118    |     121.6      |    111.3     |    2.5          | slurm-13550518.out
| 256x256x256   |  512   |  16   |      4      |   8                 | 124    |     135.3      |    123.1     |    1.9          | slurm-13550507.out
| 512x512x512   |  4096  |  128  |      5      |   64                | OoM    |     --         |    --        |    --           | slurm-13550522.out
