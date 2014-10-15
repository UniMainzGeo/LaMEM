function Sol_vec_fine = Prolongate(Sol_vec_coarse,level_coarse, MG_Info, Method)
% Interpolates the error from a coarser to a finer grid
%
%

% Method      =       {'Interpolation','MatrixBased'};
% Method      =       Method{2};

level_fine  =       level_coarse+1;

switch Method
    case 'Interpolation'        % use 3D interpolation with the various grids
        
        
        %% Interpolate Vx-error  [1th force balance equation ]
        Vx_coarse   =   Sol_vec_coarse(  MG_Info.Numbering{level_coarse}.Number_Vx );
        Vx_coarse   =   Vx_coarse(1:end-1,:,1:end-1);     % we have 'dummy' planes because of the way the current version of LaMEM is written.
        
        XVx_coarse  = 	MG_Info.Grids{level_coarse}.XVx;
        YVx_coarse  = 	MG_Info.Grids{level_coarse}.YVx;
        ZVx_coarse  = 	MG_Info.Grids{level_coarse}.ZVx;
        
        % Extend grid at the boundaries, to prevent NaN's from occuring.
        % This extension should actually depend on the type of boundary condition
        % with free slip requiring different conditions than no slip.
        dz                            =   ZVx_coarse(2,2,2)-ZVx_coarse(1,1,1);
        dy                            =   YVx_coarse(2,2,2)-YVx_coarse(1,1,1);
        
        XVx_coarse(2:end+1,:,2:end+1) =   XVx_coarse;
        YVx_coarse(2:end+1,:,2:end+1) =   YVx_coarse;
        ZVx_coarse(2:end+1,:,2:end+1) =   ZVx_coarse;
        Vx_coarse(2:end+1,:,2:end+1) =   Vx_coarse;
        
        % Extend grids in y-direction -----------
        XVx_coarse(1,:,:)             =     XVx_coarse(2,:,:);
        XVx_coarse(end+1,:,:)         =     XVx_coarse(end,:,:);
        
        YVx_coarse(1,:,:)             =     YVx_coarse(2,:,:)   -dy;
        YVx_coarse(end+1,:,:)         =     YVx_coarse(end,:,:) +dy;
        
        ZVx_coarse(1,:,:)             =     ZVx_coarse(2,:,:);
        ZVx_coarse(end+1,:,:)         =     ZVx_coarse(end,:,:);
        
        Vx_coarse(1,:,:)              =     Vx_coarse(2,:,:);       % for free slip
        Vx_coarse(end+1,:,:)          =     Vx_coarse(end,:,:);     % for free slip
        %---------------------------------------
        
        % Extend grids in vertical direction ---
        XVx_coarse(:,:,1)             =     XVx_coarse(:,:,2);
        XVx_coarse(:,:,end+1)         =     XVx_coarse(:,:,end);
        
        YVx_coarse(:,:,1)             =     YVx_coarse(:,:,2);
        YVx_coarse(:,:,end+1)         =     YVx_coarse(:,:,end);
        
        ZVx_coarse(:,:,1)             =     ZVx_coarse(:,:,2)   -dz;
        ZVx_coarse(:,:,end+1)         =     ZVx_coarse(:,:,end) +dz;
        
        Vx_coarse(:,:,1)              =     Vx_coarse(:,:,2);       % for free slip
        Vx_coarse(:,:,end+1)          =     Vx_coarse(:,:,end);     % for free slip
        %---------------------------------------
        
        XVx_fine    = 	MG_Info.Grids{level_fine  }.XVx;
        YVx_fine    = 	MG_Info.Grids{level_fine  }.YVx;
        ZVx_fine    = 	MG_Info.Grids{level_fine  }.ZVx;
        
        Vx_fine     =   interp3(XVx_coarse, YVx_coarse, ZVx_coarse,Vx_coarse, XVx_fine, YVx_fine, ZVx_fine,'linear');
        
        Vx_fine(:,:,end+1) = 0;   Vx_fine(end+1,:,:) = 0; % add 'dummy' planes
        
        
        
        % [Pp]        =   ProlongationMatrixesVelocity(level_fine, MG_Info);
        % sol_fine    =   Pp*Sol_vec_coarse;
        % Vx_fine     =   sol_fine(  MG_Info.Numbering{level_fine}.Number_Vx );
        
        
        
        
        %% Interpolate Vy-error  [2nd force balance equation ]
        Vy_coarse   =   Sol_vec_coarse(  MG_Info.Numbering{level_coarse}.Number_Vy );
        Vy_coarse   =   Vy_coarse(:,1:end-1,1:end-1);     % we have 'dummy' planes because of the way the current version of LaMEM is written.
        
        XVy_coarse  = 	MG_Info.Grids{level_coarse}.XVy;
        YVy_coarse  = 	MG_Info.Grids{level_coarse}.YVy;
        ZVy_coarse  = 	MG_Info.Grids{level_coarse}.ZVy;
        
        % Extend grid at the boundaries, to prevent NaN's from occuring.
        % This extension should actually depend on the type of boundary condition
        % with free slip requiring different conditions than no slip.
        dz                              =   ZVy_coarse(2,2,2)-ZVy_coarse(1,1,1);
        dx                              =   XVy_coarse(2,2,2)-XVy_coarse(1,1,1);
        
        XVy_coarse(:,2:end+1,2:end+1)   =   XVy_coarse;
        YVy_coarse(:,2:end+1,2:end+1)   =   YVy_coarse;
        ZVy_coarse(:,2:end+1,2:end+1)   =   ZVy_coarse;
        Vy_coarse(:,2:end+1,2:end+1)   =   Vy_coarse;
        
        % Extend grids in x-direction -----------
        XVy_coarse(:,1,:)             =     XVy_coarse(:,2,:)-dx;
        XVy_coarse(:,end+1,:)         =     XVy_coarse(:,end,:)+dx;
        
        YVy_coarse(:,1,:)             =     YVy_coarse(:,2,:)   ;
        YVy_coarse(:,end+1,:)         =     YVy_coarse(:,end,:) ;
        
        ZVy_coarse(:,1,:)             =     ZVy_coarse(:,2,:);
        ZVy_coarse(:,end+1,:)         =     ZVy_coarse(:,end,:);
        
        Vy_coarse(:,1,:)              =     Vy_coarse(:,2,:);       % for free slip
        Vy_coarse(:,end+1,:)          =     Vy_coarse(:,end,:);     % for free slip
        %---------------------------------------
        
        % Extend grids in vertical direction ---
        XVy_coarse(:,:,1)             =     XVy_coarse(:,:,2);
        XVy_coarse(:,:,end+1)         =     XVy_coarse(:,:,end);
        
        YVy_coarse(:,:,1)             =     YVy_coarse(:,:,2);
        YVy_coarse(:,:,end+1)         =     YVy_coarse(:,:,end);
        
        ZVy_coarse(:,:,1)             =     ZVy_coarse(:,:,2)   -dz;
        ZVy_coarse(:,:,end+1)         =     ZVy_coarse(:,:,end) +dz;
        
        Vy_coarse(:,:,1)              =     Vy_coarse(:,:,2);       % for free slip
        Vy_coarse(:,:,end+1)          =     Vy_coarse(:,:,end);     % for free slip
        %---------------------------------------
        
        XVy_fine    = 	MG_Info.Grids{level_fine  }.XVy;
        YVy_fine    = 	MG_Info.Grids{level_fine  }.YVy;
        ZVy_fine    = 	MG_Info.Grids{level_fine  }.ZVy;
        
        Vy_fine     =   interp3(XVy_coarse, YVy_coarse, ZVy_coarse,Vy_coarse, XVy_fine, YVy_fine, ZVy_fine,'linear');
        
        
        Vy_fine(:,:,end+1) = 0;   Vy_fine(:,end+1,:) = 0; % add 'dummy' planes
        
        
        
        % [Pp]        =   ProlongationMatrixesVelocity(level_fine, MG_Info);
        % sol_fine    =   Pp*Sol_vec_coarse;
        % Vy_fine     =   sol_fine(  MG_Info.Numbering{level_fine}.Number_Vy );
        %
        
        
        
        %% Interpolate Vz-error  [3rd force balance equation ]
        Vz_coarse   =   Sol_vec_coarse(  MG_Info.Numbering{level_coarse}.Number_Vz );
        Vz_coarse   =   Vz_coarse(1:end-1,1:end-1,:);     % we have 'dummy' planes because of the way the current version of LaMEM is written.
        
        XVz_coarse  = 	MG_Info.Grids{level_coarse}.XVz;
        YVz_coarse  = 	MG_Info.Grids{level_coarse}.YVz;
        ZVz_coarse  = 	MG_Info.Grids{level_coarse}.ZVz;
        
        % Extend grid at the boundaries, to prevent NaN's from occuring.
        % This extension should actually depend on the type of boundary condition
        % with free slip requiring different conditions than no slip.
        dy                              =   YVz_coarse(2,2,2)-YVz_coarse(1,1,1);
        dx                              =   XVz_coarse(2,2,2)-XVz_coarse(1,1,1);
        
        XVz_coarse(2:end+1,2:end+1,:)   =   XVz_coarse;
        YVz_coarse(2:end+1,2:end+1,:)   =   YVz_coarse;
        ZVz_coarse(2:end+1,2:end+1,:)   =   ZVz_coarse;
        Vz_coarse(2:end+1,2:end+1,:)   =   Vz_coarse;
        
        % Extend grids in x-direction -----------
        XVz_coarse(:,1,:)             =     XVz_coarse(:,2,:)-dx;
        XVz_coarse(:,end+1,:)         =     XVz_coarse(:,end,:)+dx;
        
        YVz_coarse(:,1,:)             =     YVz_coarse(:,2,:)   ;
        YVz_coarse(:,end+1,:)         =     YVz_coarse(:,end,:) ;
        
        ZVz_coarse(:,1,:)             =     ZVz_coarse(:,2,:);
        ZVz_coarse(:,end+1,:)         =     ZVz_coarse(:,end,:);
        
        Vz_coarse(:,1,:)              =     Vz_coarse(:,2,:);       % for free slip
        Vz_coarse(:,end+1,:)          =     Vz_coarse(:,end,:);     % for free slip
        %---------------------------------------
        
        % Extend grids in y-direction -----------
        XVz_coarse(1,:,:)             =     XVz_coarse(2,:,:);
        XVz_coarse(end+1,:,:)         =     XVz_coarse(end,:,:);
        
        YVz_coarse(1,:,:)             =     YVz_coarse(2,:,:)   -dy;
        YVz_coarse(end+1,:,:)         =     YVz_coarse(end,:,:) +dy;
        
        ZVz_coarse(1,:,:)             =     ZVz_coarse(2,:,:);
        ZVz_coarse(end+1,:,:)         =     ZVz_coarse(end,:,:);
        
        Vz_coarse(1,:,:)              =     Vz_coarse(2,:,:);       % for free slip
        Vz_coarse(end+1,:,:)          =     Vz_coarse(end,:,:);     % for free slip
        %---------------------------------------
        
        XVz_fine    = 	MG_Info.Grids{level_fine  }.XVz;
        YVz_fine    = 	MG_Info.Grids{level_fine  }.YVz;
        ZVz_fine    = 	MG_Info.Grids{level_fine  }.ZVz;
        
        Vz_fine     =   interp3(XVz_coarse, YVz_coarse, ZVz_coarse,Vz_coarse, XVz_fine, YVz_fine, ZVz_fine,'linear');
        
        
        Vz_fine(end+1,:,:) = 0;   Vz_fine(:,end+1,:) = 0; % add 'dummy' planes
        
        
        [Pp]        =   ProlongationMatrixesVelocity(level_fine, MG_Info);
        sol_fine    =   Pp*Sol_vec_coarse;
        Vz_fine     =   sol_fine(  MG_Info.Numbering{level_fine}.Number_Vz );
        
        
        
        %% Interpolate P-error  [incompressibility equation ]
        %We use injection, rather than linear interpolation for this one
        P_coarse        =   Sol_vec_coarse(  MG_Info.Numbering{level_coarse}.Number_P );
        
        XP_coarse       = 	MG_Info.Grids{level_coarse}.XP;
        YP_coarse       = 	MG_Info.Grids{level_coarse}.YP;
        ZP_coarse       = 	MG_Info.Grids{level_coarse}.ZP;
        
        XP_fine         = 	MG_Info.Grids{level_fine  }.XP;
        YP_fine         = 	MG_Info.Grids{level_fine  }.YP;
        ZP_fine         = 	MG_Info.Grids{level_fine  }.ZP;
        
        P_fine          =   interp3(XP_coarse, YP_coarse, ZP_coarse,P_coarse, XP_fine, YP_fine, ZP_fine,'nearest');
        
        
        % The boundaries are outside the validity region of interp3, and is set
        % manually
        P_fine(:,1,:)   =   P_fine(:,2,:);
        P_fine(:,end,:) =   P_fine(:,end-1,:);
        
        P_fine(1,:,:)   =   P_fine(2,:,:);
        P_fine(end,:,:) =   P_fine(end-1,:,:);
        
        P_fine(:,:,1)   =   P_fine(:,:,2);
        P_fine(:,:,end) =   P_fine(:,:,end-1);
        
        % 'Inject pressure' points
        % P_coarse = P_coarse1;
        P_fine                                  = zeros(size(MG_Info.Grids{level_fine  }.XP));
        P_fine(1:2:end,1:2:end,1:2:end)         = P_coarse;
        P_fine(end,1:2:end,1:2:end)           	= P_coarse(end,:,:);
        P_fine(1:2:end,end,1:2:end)           	= P_coarse(:,end,:);
        P_fine(1:2:end,1:2:end,end)           	= P_coarse(:,:,end);
        
        
        P_fine(2:2:end-1,:,:)                   = P_fine(3:2:end,:,:);
        P_fine(:,2:2:end-1,:)                   = P_fine(:,3:2:end,:);
        P_fine(:,:,2:2:end-1)                   = P_fine(:,:,3:2:end);
        P_fine(end,end,:)                       = P_fine(end-1,end-1,:);
        P_fine(:,:,end)                         = P_fine(:,:,end-1);
        
        
        
        %% Fill matrix
        Sol_vec_fine( MG_Info.Numbering{level_fine}.Number_Vx(:) )  =   Vx_fine(:);
        Sol_vec_fine( MG_Info.Numbering{level_fine}.Number_Vy(:) )  =   Vy_fine(:);
        Sol_vec_fine( MG_Info.Numbering{level_fine}.Number_Vz(:) )  =   Vz_fine(:);
        Sol_vec_fine( MG_Info.Numbering{level_fine}.Number_P(:) )   =   P_fine(:);
        
        Sol_vec_fine = Sol_vec_fine(:);
        
        
        
        
    case 'MatrixBased'
        
        if isfield(MG_Info.Numerics{level_fine},'P')
            P = MG_Info.Numerics{level_fine}.P;
        else
            % Do it using a prolongation matrix for pressure points:
            [Pv]            =   ProlongationMatrixesVelocity(level_fine, MG_Info);
            [Pp]            =   ProlongationMatrixPressure(level_fine, MG_Info);
            P               =   Pv+Pp;
        end
        
       
        Sol_vec_fine    =   P*Sol_vec_coarse;
        
        
        
end


%% Set BC's to zero here (as we typically deal with errors that are interpolated)
% Sol_vec_fine( MG_Info.BC{level_fine}.DirichletNodes ) = 0;

