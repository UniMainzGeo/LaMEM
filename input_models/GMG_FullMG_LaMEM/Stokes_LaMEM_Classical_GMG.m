% Calling routine to use LaMEM in combination with a  Classical Geometric
% Multigrid Solver that employs a coupled multigrid (for both velocity and
% pressure), but classical

clear

addpath ../../matlab/

% Specify resolution of coarsest grid
res_coarse      =   4;      % # of elements in x/y/z direction
levels          =   3;   	% # of multigrid levels. We coarsen by a factor 2 every time.

% Specity viscosity structure
ViscosityBlock  =   1;      % Viscosity of the falling block

% Multigrid parameters
omega_1         =   0.5;        % Smoothening factor for downward steps
omega_2         =   0.5;        % Same during upwards steps
nu1             =   20;         % number of smoothening steps going down
nu2             =   nu1;

% Run routine to create stiffness matrixes at every multigrid level
CreateStiffnessMatrixes_LaMEM

% load input matrixes for all multigrid levels [created from LaMEM with the routine above]
load StiffnessMatrixes

% print information:
disp(' '); disp(' '); disp(' ')
disp(['Coupled Geometric Multigrid solver ...'])
for i=1:levels
    disp(['Level ',num2str(i),' : ',num2str(size(MG_Info.Numbering{i}.Number_P)),' elements'])
end
disp(' ');


Rhs_vec                                     =   MG_Info.Matrixes{levels}.Rhs_vec;
Sol_vec                                     =   zeros(size(Rhs_vec));
Sol_vec(MG_Info.BC{levels}.DirichletNodes)  =   MG_Info.BC{levels}.DirichletValues;

% perform V-cycles
iter        =   1;
max_error   =   realmax;
while max_error>1e-6
    
    %     %---------------------------------------------
    %     % Method 1: Use a Jacobi iterative solver:
    %     Sol_vec     = jacobi(Sol_vec, Rhs_vec, MG_Info, levels, omega_1, omega_2, nu1);
    %     %---------------------------------------------
    
    %     ---------------------------------------------
    %     Method 2:  Use MG:
    %     (1) Perform one V-cycle
    Sol_vec         =   V_cycle(MG_Info,Sol_vec, Rhs_vec,levels,nu1, nu2, omega_1, omega_2,'Interpolation');
    %     ---------------------------------------------
    
    % (2) Compute residual
    Res             =   Residual(MG_Info.Matrixes{levels}.A, Sol_vec, Rhs_vec, MG_Info, levels);  % compute residual
    Res(MG_Info.BC{levels}.DirichletNodes) = 0;
    
    %  Compute divergence
    Div             =   Res(MG_Info.Numbering{levels}.numV:end);            % P-part of residual vector
    
    max_div         =   max(abs(Div));
    max_error       =   max(abs(Res));
    max_error       =   max_div;
    
    res_iter(iter)  =   max_error;
    
    disp(['Iteration ',num2str(iter), '  Max Res.=',num2str( max(abs(Res))),'; norm(Residual)=',num2str(norm(Res)),' max. div=' ,num2str(max_div)])
    iter            =   iter + 1;
    
end

