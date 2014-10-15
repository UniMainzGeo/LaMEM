% Creates the stiffness matrixes at several levels for a falling  block
% setup, by calling LaMEM several times.


addpath ../../../matlab/
addpath ./MULTIGRID_FULL/

% res_coarse  =   4;      % resolution at coarsest grid
% levels      =   3;  

ViscosityBlock = 1;
% Create stiffness matrixes for the same setup, with different resolutions
for ilevel=1:levels
    
    nel 	=   res_coarse*2^(ilevel-1);
    Nx      =   nel+1;
    Ny      =   nel+1;
    Nz      =   nel+1;
    
    
    str     =   ['!../../bin/opt/LaMEM -restart 0 -save_timesteps 0 -save_breakpoints -1 -mumax ',num2str(ViscosityBlock),' -StokesSolver 4 -DumpStiffnessMatrixes -numProd 1 -nel ',num2str(Nx-1),',',num2str(Ny-1),',',num2str(Nz-1)];
    
    % run falling block setup
    eval(str);
    
    % read stiffness matrixes
    [VV,VP, PV, PP, approx_S,InfoVec]       =   PetscBinaryRead('stiffness.dat');
    [f,g]                                   =   PetscBinaryRead('rhs_vector.dat');

    % Create coordinate grids
    nel_x       =   InfoVec(1);
    nel_y    	=   InfoVec(2);
    nel_z   	=   InfoVec(3);
    
    [X,Y,Z] 	=   meshgrid(linspace(0,1,nel_x+1), linspace(0,1,nel_y+1), linspace(0,1,nel_z+1)); % corners
    XVy         =  (X(:,2:end,2:end) + X(:,1:end-1,1:end-1) + X(:,2:end,1:end-1) + X(:,1:end-1,2:end))/4;
    YVy         =  (Y(:,2:end,2:end) + Y(:,1:end-1,1:end-1) + Y(:,2:end,1:end-1) + Y(:,1:end-1,2:end))/4;
    ZVy         =  (Z(:,2:end,2:end) + Z(:,1:end-1,1:end-1) + Z(:,2:end,1:end-1) + Z(:,1:end-1,2:end))/4;
    
    XVx         =  (X(2:end,:,2:end) + X(1:end-1,:,1:end-1) + X(2:end,:,1:end-1) + X(1:end-1,:,2:end))/4;
    YVx         =  (Y(2:end,:,2:end) + Y(1:end-1,:,1:end-1) + Y(2:end,:,1:end-1) + Y(1:end-1,:,2:end))/4;
    ZVx         =  (Z(2:end,:,2:end) + Z(1:end-1,:,1:end-1) + Z(2:end,:,1:end-1) + Z(1:end-1,:,2:end))/4;
    
    XVz         =  (X(2:end,2:end,:) + X(1:end-1,1:end-1,:) + X(2:end,1:end-1,:) + X(1:end-1,2:end,:))/4;
    YVz         =  (Y(2:end,2:end,:) + Y(1:end-1,1:end-1,:) + Y(2:end,1:end-1,:) + Y(1:end-1,2:end,:))/4;
    ZVz         =  (Z(2:end,2:end,:) + Z(1:end-1,1:end-1,:) + Z(2:end,1:end-1,:) + Z(1:end-1,2:end,:))/4;
    
    XP          =  (X(2:end,2:end,2:end  )+X(2:end,1:end-1,2:end  ) + X(1:end-1,1:end-1,2:end  )+X(1:end-1,2:end,2:end) + ...
        X(2:end,2:end,1:end-1)+X(2:end,1:end-1,1:end-1) + X(1:end-1,1:end-1,1:end-1)+X(1:end-1,2:end,1:end-1))/8;
    YP          =  (Y(2:end,2:end,2:end  )+Y(2:end,1:end-1,2:end  ) + Y(1:end-1,1:end-1,2:end  )+Y(1:end-1,2:end,2:end) + ...
        Y(2:end,2:end,1:end-1)+Y(2:end,1:end-1,1:end-1) + Y(1:end-1,1:end-1,1:end-1)+Y(1:end-1,2:end,1:end-1))/8;
    ZP          =  (Z(2:end,2:end,2:end  )+Z(2:end,1:end-1,2:end  ) + Z(1:end-1,1:end-1,2:end  )+Z(1:end-1,2:end,2:end) + ...
        Z(2:end,2:end,1:end-1)+Z(2:end,1:end-1,1:end-1) + Z(1:end-1,1:end-1,1:end-1)+Z(1:end-1,2:end,1:end-1))/8;
    
    
    
    
    %----------------------------------------------------------------------
    %
    % Form element numbering etc.
    
    nnode_x                     =   InfoVec(4);
    nnode_y                     =   InfoVec(5);
    nnode_z                     =   InfoVec(6);
    
    node        =   1;
    dof_number  =   1;
    for iz=1:nnode_z;
        for iy =1:nnode_y;
            for ix=1:nnode_x;
                
                
                NodeNumber(iy,ix,iz)   =   node;
                
                for dof=1:3
                    DOF_Number(iy,ix,iz,dof)    =    dof_number;
                    dof_number                  =    dof_number+1;
                end
                
                node                   =   node+1;
            end
        end
    end
    
    % Set boundary conditions
    Number_Vx           = squeeze(DOF_Number(:,:,:,1));
    Number_Vy           = squeeze(DOF_Number(:,:,:,2));
    Number_Vz           = squeeze(DOF_Number(:,:,:,3));
    
    
    numV = max(DOF_Number(:));
    numP = nel_x*nel_y*nel_z;
    node = numV+1;
    for iz=1:nel_z;
        for iy=1:nel_y;
            for ix=1:nel_x;
                Number_P(ix,iy,iz)   =   node;
                node                   =   node+1;
            end
        end
    end
    %
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Set BC's
    DirichletNodes  = [];
    DirichletValues = [];
    
    % Vx
    num1            =   Number_Vx(:  ,1  ,:);   % Vx=0, left boundary
    num2            =   Number_Vx(:  ,end,:);   % Vx=0, right boundary
    num3            =   Number_Vx(end,:  ,:);   % dummy plane
    num4            =   Number_Vx(:  ,:  ,end); % dummy plane
    
    DirichletNodes  =   [DirichletNodes;  num1(:); num2(:); num3(:); num4(:)];
    DirichletValues =   [DirichletValues; num1(:)*0; num2(:)*0; num3(:)*0; num4(:)*0];
    
    % Vy
    num1            =   Number_Vy(1  ,:  ,:);   % Vy=0, front boundary
    num2            =   Number_Vy(end,:  ,:);   % Vy=0, back boundary
    num3            =   Number_Vy(:   ,end,:);  % dummy plane
    num4            =   Number_Vy(:  ,:  ,end); % dummy plane
    
    DirichletNodes  =   [DirichletNodes;  num1(:); num2(:); num3(:); num4(:)];
    DirichletValues =   [DirichletValues; num1(:)*0; num2(:)*0; num3(:)*0; num4(:)*0];
    
    
    % Vz
    num1            =   Number_Vz(:  ,:  ,1  );   % Vz=0, bottom boundary
    num2            =   Number_Vz(:  ,:  ,end);   % Vz=0, top boundary
    num3            =   Number_Vz(:  ,end,:  );   % dummy plane
    num4            =   Number_Vz(end,:  ,:  );   % dummy plane
    
    DirichletNodes  =   [DirichletNodes;  num1(:); num2(:); num3(:); num4(:)];
    DirichletValues =   [DirichletValues; num1(:)*0; num2(:)*0; num3(:)*0; num4(:)*0];
    
   
    Free = 1:(numV+numP);
    Free(DirichletNodes) = [];  % non-boundary nodes
    %
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Build matrixes and preconditioning 
    
    A                       =   [VV VP; PV PP*0];
    Rhs_vec         =   [f; g];
    
    diagS                               = -approx_S;               % 1/eta diagonal matrix
    A_taras                             = [diag(diag(VV)),    VP*0;   PV*0, diagS ];    % used as a jacobi smoothener later
    A_precondition                      = A_taras;

    %----------------------------------------------------------------------
    %
    % Store data for the current MG level
    MG_Info.Grids{ilevel}.XVx               =   XVx;
    MG_Info.Grids{ilevel}.YVx               =   YVx;
    MG_Info.Grids{ilevel}.ZVx               =   ZVx;
    MG_Info.Grids{ilevel}.XVy               =   XVy;
    MG_Info.Grids{ilevel}.YVy               =   YVy;
    MG_Info.Grids{ilevel}.ZVy               =   ZVy;
    MG_Info.Grids{ilevel}.XVz               =   XVz;
    MG_Info.Grids{ilevel}.YVz               =   YVz;
    MG_Info.Grids{ilevel}.ZVz               =   ZVz;
    MG_Info.Grids{ilevel}.XP                =   XP;
    MG_Info.Grids{ilevel}.YP                =   YP;
    MG_Info.Grids{ilevel}.ZP                =   ZP;
    
    MG_Info.Numbering{ilevel}.Number_Vx     =   Number_Vx;
    MG_Info.Numbering{ilevel}.Number_Vy     =   Number_Vy;
    MG_Info.Numbering{ilevel}.Number_Vz     =   Number_Vz;
    MG_Info.Numbering{ilevel}.Number_P      =   Number_P;
    MG_Info.Numbering{ilevel}.numV          =   numV;
    MG_Info.Numbering{ilevel}.numP          =   numP;
    
    MG_Info.Numerics{ilevel}.Nx             =   Nx;
    MG_Info.Numerics{ilevel}.Nz             =   Ny;
    MG_Info.Numerics{ilevel}.Nz             =   Nz;
    
    %     MG_Info.Numerics{ilevel}.penalty_fac    =   penalty_fac;
    
    
    MG_Info.Matrixes{ilevel}.A              = A;
    %     MG_Info.Matrixes{ilevel}.A_free         = A_free;
    MG_Info.Matrixes{ilevel}.Rhs_vec        = Rhs_vec;
    MG_Info.Matrixes{ilevel}.A_precondition = A_precondition;
    %     MG_Info.Matrixes{ilevel}.A_precondition_free = A_precondition_free;
    
    MG_Info.Matrixes{ilevel}.F              = VV;
    MG_Info.Matrixes{ilevel}.Bt             = VP;
    MG_Info.Matrixes{ilevel}.B              = PV;
    MG_Info.Matrixes{ilevel}.C              = PP;
    
    MG_Info.BC{ilevel}.DirichletNodes    	= DirichletNodes;
    MG_Info.BC{ilevel}.DirichletValues    	= DirichletValues;
    MG_Info.BC{ilevel}.Free                 = Free;
    
    %
    %----------------------------------------------------------------------
    
    
end

% cleanup directory
!rm *.dat
!rm *.info

% Store data into mat-file
save StiffnessMatrixes -v7.3


