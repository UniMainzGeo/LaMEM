% Creates the stiffness matrixes at the finest level and than use a
% galerkin interpolation to derive the other levels


%clear

addpath ../../../matlab/
addpath ./MULTIGRID_FULL/


%% Call LaMEM at the finest level
str     =   ['!../../bin/opt/LaMEM -restart 0 -save_timesteps 0 -save_breakpoints -1 -mumax ',num2str(ViscosityBlock),'-BC.Left 0 -StokesSolver 4 -DumpStiffnessMatrixes -numProd 1 -nel ',num2str(Nx),',',num2str(Ny),',',num2str(Nz)]

% str     =   ['!../../../bin/LaMEM -ParamFile Detachment_Folding_6a.dat -save_timesteps 0 -restart 0 -save_breakpoints -1 -StokesSolver 4 -DumpStiffnessMatrixes -numProd 1 -nel ',num2str(Nx),',',num2str(Ny),',',num2str(Nz)]

% str     =   ['!../../../bin/LaMEM  -ParamFile Subduction3D_FDSTAG.dat -restart 0 -save_timesteps 0 -save_breakpoints -1 -StokesSolver 4 -DumpStiffnessMatrixes -numProd 1 -nel ',num2str(Nx),',',num2str(Ny),',',num2str(Nz),' -NumPartX ',num2str(2),' -NumPartY ',num2str(2),' -NumPartZ ',num2str(2)];
        

% run falling block setup
eval(str);


% A first step is to create the grids numbering and restriction/prolongation operators for all levels
for ilevel=levels:-1:1
    
    
    if ilevel==levels
        % we are at the finest level
        nel_x = Nx;
        nel_y = Ny;
        nel_z = Nz;
    
    else
        nel_x = nel_x/2;
        nel_y = nel_y/2;
        nel_z = nel_z/2;
    
    end
    nnode_x = nel_x+1;
    nnode_y = nel_y+1;
    nnode_z = nel_z+1;
    
    % Node numbering
    DOF_Number  =   [];
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
    
    Number_P= [];
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
    DirichletNodes  =   [];
    DirichletValues =   [];
    
    % Vx
    num1            =   Number_Vx(:  ,1  ,:);   % Vx=0, left boundary
    num2            =   Number_Vx(:  ,end,:);   % Vx=0, right boundary
    num3            =   Number_Vx(end,:  ,:);   % dummy plane
    num4            =   Number_Vx(:  ,:  ,end); % dummy plane
    
    DirichletNodes  =   [DirichletNodes;  num1(:); num2(:); num3(:); num4(:)];
  %  DirichletValues =   [DirichletValues; num1(:)*0; num2(:)*0; num3(:)*0; num4(:)*0];
    
    % Vy
    num1            =   Number_Vy(1  ,:  ,:);   % Vy=0, front boundary
    num2            =   Number_Vy(end,:  ,:);   % Vy=0, back boundary
    num3            =   Number_Vy(:   ,end,:);  % dummy plane
    num4            =   Number_Vy(:  ,:  ,end); % dummy plane
    
    DirichletNodes  =   [DirichletNodes;  num1(:); num2(:); num3(:); num4(:)];
   % DirichletValues =   [DirichletValues; num1(:)*0; num2(:)*0; num3(:)*0; num4(:)*0];
    
    
    % Vz
    num1            =   Number_Vz(:  ,:  ,1  );   % Vz=0, bottom boundary
    num2            =   Number_Vz(:  ,:  ,end);   % Vz=0, top boundary
    num3            =   Number_Vz(:  ,end,:  );   % dummy plane
    num4            =   Number_Vz(end,:  ,:  );   % dummy plane
    
    DirichletNodes  =   [DirichletNodes;  num1(:); num2(:); num3(:); num4(:)];
    %DirichletValues =   [DirichletValues; num1(:)*0; num2(:)*0; num3(:)*0; num4(:)*0];
    
    
%    DirichletValues = Rhs_vec(DirichletNodes);
    
    Free                    = 1:(numV+numP);
    Free(DirichletNodes)    = [];  % non-boundary nodes
    %
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    

    %----------------------------------------------------------------------
    %
    % Store data for the current MG level
    MG_Info.Numbering{ilevel}.Number_Vx     =   Number_Vx;
    MG_Info.Numbering{ilevel}.Number_Vy     =   Number_Vy;
    MG_Info.Numbering{ilevel}.Number_Vz     =   Number_Vz;
    MG_Info.Numbering{ilevel}.Number_P      =   Number_P;
    MG_Info.Numbering{ilevel}.numV          =   numV;
    MG_Info.Numbering{ilevel}.numP          =   numP;
    
    MG_Info.Numerics{ilevel}.Nx             =   Nx;
    MG_Info.Numerics{ilevel}.Nz             =   Ny;
    MG_Info.Numerics{ilevel}.Nz             =   Nz;
   
    MG_Info.BC{ilevel}.DirichletNodes    	= DirichletNodes;
  %  MG_Info.BC{ilevel}.DirichletValues    	= DirichletValues;
    MG_Info.BC{ilevel}.Free                 = Free;
   
    
end





for ilevel=levels:-1:1
    
    
    if ilevel==levels
        % we are at the finest level
        
        
        
        % read stiffness matrixes
        [VV,VP, PV, PP, approx_S,InfoVec]       =       PetscBinaryRead('stiffness.dat');
        [f,g]                                   =       PetscBinaryRead('rhs_vector.dat');
        
        % Create coordinate grids
        nel_x       =   InfoVec(1);     nel_y    	=   InfoVec(2); nel_z   	=   InfoVec(3);
        
        %----------------------------------------------------------------------
        %
        % Build matrixes and preconditioner
        A                                       =       [VV VP; PV PP*0];
        Rhs_vec                                 =       [f; g];
        
        diagS                                   =       -approx_S;               % 1/eta diagonal matrix
        A_precondition                       	=       [diag(diag(VV)),    VP*0;   PV*0, diagS ];    % used as a jacobi smoothener later
        
    else
        
%         R=P';
%         %P = R';
%         
        
        
        A_coarse        =   R*A*P;
       
        
        A_p_coarse      =   R*A_precondition*P;
        Rhs_coarse      =   R*Rhs_vec(:);
        
        
        A               =   A_coarse;
        A_precondition  =   A_p_coarse;
        Rhs_vec         =   Rhs_coarse;
        
    end
    
    if ilevel>1
        % Compute restriction and prolongation operators for the current level
        [Pv]                =   ProlongationMatrixesVelocity(ilevel, MG_Info);
        [Pp]                =   ProlongationMatrixPressure(ilevel, MG_Info);
        P                   =   Pv+Pp;
        
        [R]                 =   RestrictionMatrixesVelocityPressure(ilevel, MG_Info);
        
        MG_Info.Numerics{ilevel}.P              =   P;  % prolongation matrix
        MG_Info.Numerics{ilevel}.R              =   R;  % restriction matrix
    end
    
        
    MG_Info.Matrixes{ilevel}.A              = A;
    MG_Info.Matrixes{ilevel}.Rhs_vec        = Rhs_vec;
    MG_Info.Matrixes{ilevel}.A_precondition = A_precondition;
    
    % Set values @ Dirichlet nodes
    MG_Info.BC{ilevel}.DirichletValues    	= Rhs_vec(MG_Info.BC{ilevel}.DirichletNodes);
    
end


% cleanup directory
!rm rhs_vector.dat
!rm stiffness.dat
!rm StiffnessMatrixesRhs_Binary.dat
!rm *.info

% Store data into mat-file
save StiffnessMatrixes -v7.3


