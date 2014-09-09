%ModifyInitialMesh_DetachmentFolding
%
% (1) Reads in the initial mesh (in parallel)
% (2) Modifies the geometry
% (3) Writes it to file
%
%

% NOTE:  The mesh that is read from disk is parallel. The mesh that is
% written out to disk is NOT parallel. This is not scalable. If this
% becomes an issue we'll change it.


clear

NewInitialMeshName = 'Test_InitialMesh.input';

% Read in the initial mesh 
% [X,Y,Z, info]	=   ReadInitialMesh('InitialMesh_DetFolding');

[X,Y,Z, info]	=   ReadInitialMesh('InitialMesh_DetFolding');

% Visualize the initial mesh
figure(1), clf
VisualizeInitialMesh(X,Y,Z)

%==========================================================================
% MODIFY THE MESH
ind_saltInterface = 13;

Step_Location   =   (min(X(:)) + max(X(:)))/2;     % x-location [m]
StepAmplitude    =   1000;                          %            [m]
StepWidth       =   1;
x               =   squeeze(X(:,:,1)) ;

% Create an ATAN based step (smooth) function
StepAddition    =   (atan((x-Step_Location)/StepWidth)/pi + 0.5)*StepAmplitude;
% StepAddition = x./max(x(:))*StepAmplitude



for iz=1:ind_saltInterface
    fac = (ind_saltInterface-iz)/(ind_saltInterface-1);      % 0 at salt interface and 1 at bottom of computational domain

    
    Z(:,:,iz) = Z(:,:,iz) + StepAddition*fac;
end


%
%==========================================================================





% Write the new initial mesh to file
WriteInitialMesh(X,Y,Z,info, NewInitialMeshName);

figure(2), clf
VisualizeInitialMesh(X,Y,Z)




% NOTE:
% In your LaMEM input file, you'll have to make sure that 
% the following lines are present: 
%
% InitialMeshFromFile 	=   1					
% InitialMeshFileName 	=	./InitialMesh/Test_InitialMesh.input
%