function Sol_vec_coarse = Restrict(Sol_vec_fine,level_fine, MG_Info, Method)
% Interpolates the error from a finer to a coarser grid
%
%

% Method      =   {'Interpolation','MatrixBased'};
% Method      =   Method{2};
level_coarse=   level_fine - 1;
        
switch Method
    case 'Interpolation'        % use 3D interpolation with the various grids 
        
        %--------------------------------------------------------------------------
        % Interpolate Vx-error  [1th force balance equation ]
        XVx_fine    = 	MG_Info.Grids{level_fine  }.XVx;
        YVx_fine    = 	MG_Info.Grids{level_fine  }.YVx;
        ZVx_fine    = 	MG_Info.Grids{level_fine  }.ZVx;
        Vx_fine     =   Sol_vec_fine(  MG_Info.Numbering{level_fine}.Number_Vx );
        Vx_fine     =   Vx_fine(1:end-1,:,1:end-1);     % we have 'dummy' planes because of the way the current version of LaMEM is written.
        
        XVx_coarse  = 	MG_Info.Grids{level_coarse}.XVx;
        YVx_coarse  = 	MG_Info.Grids{level_coarse}.YVx;
        ZVx_coarse  = 	MG_Info.Grids{level_coarse}.ZVx;
        Vx_coarse   =   interp3(XVx_fine, YVx_fine, ZVx_fine,Vx_fine, XVx_coarse, YVx_coarse, ZVx_coarse,'linear');
        
        Vx_coarse(:,:,end+1) = 0;   Vx_coarse(end+1,:,:) = 0; % add 'dummy' planes
        
        
        
        %--------------------------------------------------------------------------
        % Interpolate Vy-error  [2nd force balance equation ]
        XVy_fine    = 	MG_Info.Grids{level_fine  }.XVy;
        YVy_fine    = 	MG_Info.Grids{level_fine  }.YVy;
        ZVy_fine    = 	MG_Info.Grids{level_fine  }.ZVy;
        Vy_fine     =   Sol_vec_fine(  MG_Info.Numbering{level_fine}.Number_Vy );
        Vy_fine     =   Vy_fine(:,1:end-1,1:end-1);     % we have 'dummy' planes because of the way the current version of LaMEM is written.
        
        XVy_coarse  = 	MG_Info.Grids{level_coarse}.XVy;
        YVy_coarse  = 	MG_Info.Grids{level_coarse}.YVy;
        ZVy_coarse  = 	MG_Info.Grids{level_coarse}.ZVy;
        
        Vy_coarse   =   interp3(XVy_fine, YVy_fine, ZVy_fine,Vy_fine, XVy_coarse, YVy_coarse, ZVy_coarse,'linear');
        Vy_coarse(:,:,end+1) = 0;   Vy_coarse(:,end+1,:) = 0; % add 'dummy' planes
        
        
        %--------------------------------------------------------------------------
        % Interpolate Vz-error  [3rd force balance equation ]
        XVz_fine    = 	MG_Info.Grids{level_fine  }.XVz;
        YVz_fine    = 	MG_Info.Grids{level_fine  }.YVz;
        ZVz_fine    = 	MG_Info.Grids{level_fine  }.ZVz;
        Vz_fine     =   Sol_vec_fine(  MG_Info.Numbering{level_fine}.Number_Vz );
        Vz_fine     =   Vz_fine(1:end-1,1:end-1,:);     % we have 'dummy' planes because of the way the current version of LaMEM is written.
        
        XVz_coarse  = 	MG_Info.Grids{level_coarse}.XVz;
        YVz_coarse  = 	MG_Info.Grids{level_coarse}.YVz;
        ZVz_coarse  = 	MG_Info.Grids{level_coarse}.ZVz;
        
        Vz_coarse   =   interp3(XVz_fine, YVz_fine, ZVz_fine,Vz_fine, XVz_coarse, YVz_coarse, ZVz_coarse,'linear');
        
        
        Vz_coarse(end+1,:,:) = 0;   Vz_coarse(:,end+1,:) = 0; % add 'dummy' planes
        
        
        
        %--------------------------------------------------------------------------
        % Interpolate P-error  (continuity equation)
        %
        % Note: we will have to weigh this equation in the future (for faster
        % convergence if using variable viscosity.
        %
        P_fine      =   Sol_vec_fine(  MG_Info.Numbering{level_fine}.Number_P );
        XP_fine     = 	MG_Info.Grids{level_fine  }.XP;
        YP_fine     = 	MG_Info.Grids{level_fine  }.YP;
        ZP_fine     = 	MG_Info.Grids{level_fine  }.ZP;
        XP_coarse   = 	MG_Info.Grids{level_coarse}.XP;
        YP_coarse   = 	MG_Info.Grids{level_coarse}.YP;
        ZP_coarse   = 	MG_Info.Grids{level_coarse}.ZP;
        
        P_coarse    =   interp3(XP_fine,YP_fine, ZP_fine,P_fine, XP_coarse, YP_coarse, ZP_coarse,'nearest');
        
        
        % Fill matrix
        Sol_vec_coarse=   zeros(size(MG_Info.Matrixes{level_coarse}.Rhs_vec));
        Sol_vec_coarse( MG_Info.Numbering{level_coarse}.Number_Vx(:) )  =   Vx_coarse(:);
        Sol_vec_coarse( MG_Info.Numbering{level_coarse}.Number_Vy(:) )  =   Vy_coarse(:);
        Sol_vec_coarse( MG_Info.Numbering{level_coarse}.Number_Vz(:) )  =   Vz_coarse(:);
        Sol_vec_coarse( MG_Info.Numbering{level_coarse}.Number_P(:) )   =   P_coarse(:);
        
        
    case 'MatrixBased'
        
        
        % Do the same as above but with matrixes. In doing this, we assume
        % a factor coarsening between the two grids
        if isfield(MG_Info.Numerics{level_fine},'R')
            R 	=   MG_Info.Numerics{level_fine}.R;
        else
            R 	=   RestrictionMatrixesVelocityPressure(level_fine, MG_Info);
        end
        
        % Restrict to coarse grid
        Sol_vec_coarse      =   R*Sol_vec_fine;

    otherwise
        error('unknown method')
        
end

% Sol_vec_coarse( MG_Info.BC{level_coarse}.DirichletNodes ) = 0;
