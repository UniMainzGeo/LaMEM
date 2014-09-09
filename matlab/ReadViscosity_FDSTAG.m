%ReadViscosity_FDSTAG
%
% Reads various viscosity fields from disk, for an FDSTAG setup.

clear


addpath('./LaMEM_MATLAB')

% Read DA data
Data        =   PetscReadBinaryMatlab('Viscosity_Center_FDSTAG.dat');
Mu_cen      =   Data.ViscosityCenter.data;
Rho_cen     =   Data.DensityCenter.data;
Phase0_cen  =   Data.Phase_0.data;
Phase1_cen  =   Data.Phase_1.data;



% Data    = PetscReadBinaryMatlab('Density_Vz_FDSTAG.dat');
% Rho_Vz  = Data.VZ_DENSITY.data;
 
Data    = PetscReadBinaryMatlab('Viscosity_XY_FDSTAG.dat');
Eta_XY  = Data.Viscosity_XY.data;
Phase0_XY  =   Data.Phase_0.data;
Phase1_XY  =   Data.Phase_1.data;


Data    = PetscReadBinaryMatlab('Viscosity_XZ_FDSTAG.dat');
Eta_XZ  = Data.Viscosity_XZ.data;
Phase0_XZ  =   Data.Phase_0.data;
Phase1_XZ  =   Data.Phase_1.data;

Data    = PetscReadBinaryMatlab('Viscosity_YZ_FDSTAG.dat');
Eta_YZ  = Data.Viscosity_YZ.data;
Phase0_YZ  =   Data.Phase_0.data;
Phase1_YZ  =   Data.Phase_1.data;
