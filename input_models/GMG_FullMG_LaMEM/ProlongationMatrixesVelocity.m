function [Pp] = ProlongationMatrixesVelocity(level_fine, MG_Info);
% Computes prolongation matrixes for the velocity equations.

% level_fine      = 2;

level_coarse    = level_fine-1;


ind_i = []; ind_j = []; val   = [];

if 1==1
    %% Vx-velocity
    NumVx_coarse     =  MG_Info.Numbering{level_coarse}.Number_Vx;
    NumVx_fine       =  MG_Info.Numbering{level_fine}.Number_Vx;
    
    % Take care of boundary conditions
    NumVx_coarse(2:end+1,:,2:end+1) = NumVx_coarse;
    NumVx_coarse(1,:,:)             = NumVx_coarse(2,:,:);
    NumVx_coarse(end+1,:,:)     	= NumVx_coarse(end,:,:);
    NumVx_coarse(end-1,:,:)     	= NumVx_coarse(end-2,:,:);
    
    NumVx_coarse(:,:,1)             = NumVx_coarse(:,:,2);
    NumVx_coarse(:,:,end+1)     	= NumVx_coarse(:,:,end);
    NumVx_coarse(:,:,end-1)     	= NumVx_coarse(:,:,end-2);
    
    NumVx_sign                      =   ones(size(NumVx_coarse));
    
    % % in case pf no-slip, Vx=0 at boundary so sign should be negative
    % NumVx_sign(1,:,:)               =   -1;
    % NumVx_sign(end,:,:)             =   -1;
    % NumVx_sign(:,:,1)               =   -1;
    % NumVx_sign(:,:,end)             =   -1;
    
    
   
    for iiy=1:size(NumVx_fine,1)-1     % last plane is dummy one
        for iix=1:size(NumVx_fine,2)
            for iiz=1:size(NumVx_fine,3)-1 % last plane is dummy one
                
                ix = iix;
                iy = iiy;
                iz = iiz;
                
                % Vx  is averaged from 12 surrounding cells
                ind_fine    =   NumVx_fine(iy,ix,iz);
                ind_coarse    = [];               ind_val = [];
                
                
                if mod(ix-1,2)==0
                    % even - bilinear interpolation along the x-plane
                    
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)  ,fix(ix/2+1),fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)+1,fix(ix/2+1),fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)+1,fix(ix/2+1),fix(iz/2+1)+1)];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)  ,fix(ix/2+1),fix(iz/2+1)+1)];
                    
                    sign            = [     NumVx_sign(fix(iy/2+1)  ,fix(ix/2+1),fix(iz/2+1)  ); ...
                        NumVx_sign(fix(iy/2+1)+1,fix(ix/2+1),fix(iz/2+1)  ); ...
                        NumVx_sign(fix(iy/2+1)+1,fix(ix/2+1),fix(iz/2+1)+1); ...
                        NumVx_sign(fix(iy/2+1)  ,fix(ix/2+1),fix(iz/2+1)+1)];
                    
                    
                    
                    if      mod(iz,2)==0 & mod(iy,2)~=0
                        ind_val     =   [ 3/16 9/16  3/16 1/16];
                        
                    elseif  mod(iz,2)~=0 & mod(iy,2)~=0
                        ind_val     =   [ 1/16 3/16  9/16 3/16 ];
                        
                    elseif  mod(iz,2)~=0 & mod(iy,2)==0
                        ind_val     =   [ 3/16 1/16 3/16 9/16 ];
                        
                    elseif  mod(iz,2)==0 & mod(iy,2)==0
                        ind_val     =   [ 9/16 3/16 1/16 3/16 ];
                    end
                    
                    ind_val = ind_val(:).*sign;
                    
                    
                else
                    
                    % odd - trilinear interpolation in x,y,z direction
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)  ,fix(ix/2  ),fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)+1,fix(ix/2  ),fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)+1,fix(ix/2  ),fix(iz/2+1)+1)];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)  ,fix(ix/2  ),fix(iz/2+1)+1)];
                    
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)  ,fix(ix/2+1),fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)+1,fix(ix/2+1),fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)+1,fix(ix/2+1),fix(iz/2+1)+1)];
                    ind_coarse      =   [ind_coarse, NumVx_coarse(fix(iy/2+1)  ,fix(ix/2+1),fix(iz/2+1)+1)];
                    
                    sign = [    NumVx_sign(fix(iy/2+1)  ,fix(ix/2  ),fix(iz/2+1)  ); ...
                        NumVx_sign(fix(iy/2+1)+1,fix(ix/2  ),fix(iz/2+1)  ); ...
                        NumVx_sign(fix(iy/2+1)+1,fix(ix/2  ),fix(iz/2+1)+1); ...
                        NumVx_sign(fix(iy/2+1)  ,fix(ix/2  ),fix(iz/2+1)+1); ...
                        NumVx_sign(fix(iy/2+1)  ,fix(ix/2+1),fix(iz/2+1)  ); ...
                        NumVx_sign(fix(iy/2+1)+1,fix(ix/2+1),fix(iz/2+1)  ); ...
                        NumVx_sign(fix(iy/2+1)+1,fix(ix/2+1),fix(iz/2+1)+1); ...
                        NumVx_sign(fix(iy/2+1)  ,fix(ix/2+1),fix(iz/2+1)+1) ];
                    
                    
                    
                    % in the x-direction it is exactly in the middle
                    if mod(ix,2)==0
                        x_fac=0.5;
                    else
                        x_fac=0.5;
                    end
                    
                    if      mod(iz,2)==0 & mod(iy,2)~=0
                        ind_val     =   [[3/16 9/16  3/16 1/16]*x_fac    [ 3/16 9/16  3/16 1/16]*(1-x_fac)] ;
                        
                    elseif  mod(iz,2)~=0 & mod(iy,2)~=0
                        ind_val     =   [[1/16 3/16  9/16 3/16]*x_fac    [ 1/16 3/16  9/16 3/16 ]*(1-x_fac)];
                        
                        
                    elseif  mod(iz,2)~=0 & mod(iy,2)==0
                        ind_val     =   [[3/16 1/16 3/16 9/16 ]*x_fac    [ 3/16 1/16 3/16 9/16 ]*(1-x_fac)];
                        
                    elseif  mod(iz,2)==0 & mod(iy,2)==0
                        ind_val     =   [[9/16 3/16 1/16 3/16 ]*x_fac    [ 9/16 3/16 1/16 3/16 ]*(1-x_fac)];
                    end
                    
                    
                    ind_val = ind_val(:).*sign; % take BC's into account
                    
                end
                
                
                % Store all info in vectors, which are later added to a sparse
                % matrix
                ind_i = [ind_i; ind_fine*ones(size(ind_coarse(:)))];
                ind_j = [ind_j; ind_coarse(:)];
                val   = [val; ind_val(:)];
                
                
                
                
                
            end
        end
    end
    
    
    
    iiy=size(NumVx_fine,1);
    for iix=1:size(NumVx_fine,2)
        for iiz=1:size(NumVx_fine,3)-1 % last plane is dummy one
            
            ix = iix;
            iy = iiy;
            iz = iiz;
            
            % Vx  is averaged from 12 surrounding cells
            ind_fine    =   NumVx_fine(iy,ix,iz);
            ind_coarse  =   NumVx_coarse(end  ,fix(ix/2+1),fix(iz/2+1)  );
            
            ind_i = [ind_i; ind_fine];
            ind_j = [ind_j; ind_coarse];
            val   = [val; 1];
            
        end
    end
    
    iiz = size(NumVx_fine,3);
    for iiy=1:size(NumVx_fine,1)
        for iix=1:size(NumVx_fine,2)
            
            ix = iix;
            iy = iiy;
            iz = iiz;
            
            % Vx  is averaged from 12 surrounding cells
            ind_fine    =   NumVx_fine(iy,ix,iz);
            ind_coarse  =   NumVx_coarse(fix(iy/2+1)  ,fix(ix/2+1),end  );
            
            ind_i = [ind_i; ind_fine];
            ind_j = [ind_j; ind_coarse];
            val   = [val; 1];
            
        end
    end
    
    
end



if 1==1
    %% Vy-velocity
    NumVy_coarse     =  MG_Info.Numbering{level_coarse}.Number_Vy;
    NumVy_fine       =  MG_Info.Numbering{level_fine}.Number_Vy;
    
    % Take care of boundary conditions
    NumVy_coarse(:,2:end+1,2:end+1) = NumVy_coarse;
    NumVy_coarse(:,1,:)             = NumVy_coarse(:,2,:);
    NumVy_coarse(:,end+1,:)     	= NumVy_coarse(:,end,:);
    NumVy_coarse(:,end-1,:)     	= NumVy_coarse(:,end-2,:);
    
    NumVy_coarse(:,:,1)             = NumVy_coarse(:,:,2);
    NumVy_coarse(:,:,end+1)     	= NumVy_coarse(:,:,end);
    NumVy_coarse(:,:,end-1)     	= NumVy_coarse(:,:,end-2);
    
    NumVy_sign                      =   ones(size(NumVy_coarse));
    
    % % in case pf no-slip, Vx=0 at boundary so sign should be negative
    % NumVy_sign(:,1,:)               =   -1;
    % NumVy_sign(:,end,:)             =   -1;
    % NumVy_sign(:,:,1)               =   -1;
    % NumVy_sign(:,:,end)             =   -1;

    for iiy=1:size(NumVy_fine,1)     
        for iix=1:size(NumVy_fine,2)-1      % last plane is dummy one
            for iiz=1:size(NumVy_fine,3)-1  % last plane is dummy one
                
                ix = iix;
                iy = iiy;
                iz = iiz;
                
                % Vy  is averaged from 12 surrounding cells
                ind_fine    =   NumVy_fine(iy,ix,iz);
                ind_coarse 	=   [];               ind_val = [];
                
                
                if mod(iy-1,2)==0
                    % even - bilinear interpolation along the x-plane
                    
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2+1), fix(ix/2+1)  , fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2+1), fix(ix/2+1)+1, fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2+1), fix(ix/2+1)+1, fix(iz/2+1)+1)];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2+1), fix(ix/2+1)  , fix(iz/2+1)+1)];
                    
                    sign            = [     NumVy_sign(fix(iy/2+1), fix(ix/2+1)  , fix(iz/2+1)  ); ...
                                            NumVy_sign(fix(iy/2+1), fix(ix/2+1)+1, fix(iz/2+1)  ); ...
                                            NumVy_sign(fix(iy/2+1), fix(ix/2+1)+1, fix(iz/2+1)+1); ...
                                            NumVy_sign(fix(iy/2+1), fix(ix/2+1)  , fix(iz/2+1)+1)];
                    
                    
                    
                    if      mod(iz,2)==0 & mod(ix,2)~=0
                        ind_val     =   [ 3/16 9/16  3/16 1/16];
                        
                    elseif  mod(iz,2)~=0 & mod(ix,2)~=0
                        ind_val     =   [ 1/16 3/16  9/16 3/16 ];
                        
                    elseif  mod(iz,2)~=0 & mod(ix,2)==0
                        ind_val     =   [ 3/16 1/16 3/16 9/16 ];
                        
                    elseif  mod(iz,2)==0 & mod(ix,2)==0
                        ind_val     =   [ 9/16 3/16 1/16 3/16 ];
                    end
                    
                    ind_val = ind_val(:).*sign;
                    
                    
                else
                    
                    % odd - trilinear interpolation in x,y,z direction
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2  ), fix(ix/2+1)  , fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2  ), fix(ix/2+1)+1, fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2  ), fix(ix/2+1)+1, fix(iz/2+1)+1)];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2  ), fix(ix/2+1)  , fix(iz/2+1)+1)];
                    
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2+1), fix(ix/2+1)  , fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2+1), fix(ix/2+1)+1, fix(iz/2+1)  )];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2+1), fix(ix/2+1)+1, fix(iz/2+1)+1)];
                    ind_coarse      =   [ind_coarse, NumVy_coarse(fix(iy/2+1), fix(ix/2+1)  , fix(iz/2+1)+1)];
                    
                    
                    
                     sign           =   [   NumVy_sign(fix(iy/2  ), fix(ix/2+1)  , fix(iz/2+1)  );      ...
                                            NumVy_sign(fix(iy/2  ), fix(ix/2+1)+1, fix(iz/2+1)  );      ...
                                            NumVy_sign(fix(iy/2  ), fix(ix/2+1)+1, fix(iz/2+1)+1);      ...
                                            NumVy_sign(fix(iy/2  ), fix(ix/2+1)  , fix(iz/2+1)+1);      ...
                                            NumVy_sign(fix(iy/2+1), fix(ix/2+1)  , fix(iz/2+1)  );      ...
                                            NumVy_sign(fix(iy/2+1), fix(ix/2+1)+1, fix(iz/2+1)  );      ...
                                            NumVy_sign(fix(iy/2+1), fix(ix/2+1)+1, fix(iz/2+1)+1);      ...
                                            NumVy_sign(fix(iy/2+1), fix(ix/2+1)  , fix(iz/2+1)+1)   ];

                    % in the y-direction it is exactly in the middle
                    if mod(iy,2)==0
                        y_fac=0.5;
                    else
                        y_fac=0.5;
                    end
                    
                    if      mod(iz,2)==0 & mod(ix,2)~=0
                        ind_val     =   [[3/16 9/16  3/16 1/16]*y_fac    [ 3/16 9/16  3/16 1/16]*(1-y_fac)] ;
                        
                    elseif  mod(iz,2)~=0 & mod(ix,2)~=0
                        ind_val     =   [[1/16 3/16  9/16 3/16]*y_fac    [ 1/16 3/16  9/16 3/16 ]*(1-y_fac)];
                        
                        
                    elseif  mod(iz,2)~=0 & mod(ix,2)==0
                        ind_val     =   [[3/16 1/16 3/16 9/16 ]*y_fac    [ 3/16 1/16 3/16 9/16 ]*(1-y_fac)];
                        
                    elseif  mod(iz,2)==0 & mod(ix,2)==0
                        ind_val     =   [[9/16 3/16 1/16 3/16 ]*y_fac    [ 9/16 3/16 1/16 3/16 ]*(1-y_fac)];
                    end
                    
                    
                    ind_val = ind_val(:).*sign; % take BC's into account
                    
                end
                
                
                % Store all info in vectors, which are later added to a sparse
                % matrix
                ind_i = [ind_i; ind_fine*ones(size(ind_coarse(:)))];
                ind_j = [ind_j; ind_coarse(:)];
                val   = [val; ind_val(:)];
                
                
                
                
                
            end
        end
    end
    
    
    iix=size(NumVy_fine,2);
    for iiy=1:size(NumVy_fine,1)
        for iiz=1:size(NumVy_fine,3)-1 % last plane is dummy one
            
            ix = iix;
            iy = iiy;
            iz = iiz;
            
            % Vx  is averaged from 12 surrounding cells
            ind_fine    =   NumVy_fine(iy,ix,iz);
            ind_coarse  =   NumVy_coarse(fix(iy/2+1)  ,end,fix(iz/2+1)  );
            
            ind_i = [ind_i; ind_fine];
            ind_j = [ind_j; ind_coarse];
            val   = [val; 1];
            
        end
    end
    
    
    iiz = size(NumVy_fine,3);
    for iiy=1:size(NumVy_fine,1)
        for iix=1:size(NumVy_fine,2)
            
            ix = iix;
            iy = iiy;
            iz = iiz;
            
            % Vx  is averaged from 12 surrounding cells
            ind_fine    =   NumVy_fine(iy,ix,iz);
            ind_coarse  =   NumVy_coarse(fix(iy/2+1)  ,fix(ix/2+1),end  );
            
            ind_i = [ind_i; ind_fine];
            ind_j = [ind_j; ind_coarse];
            val   = [val; 1];
            
        end
    end
    
    
end


if 1==1
    %% Vz-velocity
    NumVz_coarse     =  MG_Info.Numbering{level_coarse}.Number_Vz;
    NumVz_fine       =  MG_Info.Numbering{level_fine}.Number_Vz;
    
    % Take care of boundary conditions
    NumVz_coarse(2:end+1,2:end+1,:) = NumVz_coarse;
    NumVz_coarse(:,1,:)             = NumVz_coarse(:,2,:);
    NumVz_coarse(:,end+1,:)     	= NumVz_coarse(:,end,:);
    NumVz_coarse(:,end-1,:)     	= NumVz_coarse(:,end-2,:);
    
    NumVz_coarse(1,:,:)             = NumVz_coarse(2,:,:);
    NumVz_coarse(end+1,:,:)     	= NumVz_coarse(end,:,:);
    NumVz_coarse(end-1,:,:)     	= NumVz_coarse(end-2,:,:);
    
    NumVz_sign                      =   ones(size(NumVz_coarse));
    
    % % in case pf no-slip, Vz=0 at boundary so sign should be negative
    % NumVz_sign(:,1,:)               =   -1;
    % NumVz_sign(:,end,:)             =   -1;
    % NumVz_sign(1,:,:)               =   -1;
    % NumVz_sign(end,:,:)             =   -1;
    
    for iiy=1:size(NumVz_fine,1)-1          % last plane is dummy one
        for iix=1:size(NumVz_fine,2)-1      % last plane is dummy one
            for iiz=1:size(NumVz_fine,3)  
                
                ix = iix;
                iy = iiy;
                iz = iiz;
                
                % Vz
                ind_fine    =   NumVz_fine(iy,ix,iz);
                ind_coarse 	=   [];               ind_val = [];
                
                
                if mod(iz-1,2)==0
                    % even - bilinear interpolation along the x-plane
                    
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1),   fix(ix/2+1)  ,  fix(iz/2+1))];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1),   fix(ix/2+1)+1,  fix(iz/2+1))];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1)+1, fix(ix/2+1)+1,  fix(iz/2+1))];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1)+1, fix(ix/2+1)  ,  fix(iz/2+1))];
                    
                    sign            = [     NumVz_sign(fix(iy/2+1),   fix(ix/2+1)  , fix(iz/2+1)); ...
                                            NumVz_sign(fix(iy/2+1),   fix(ix/2+1)+1, fix(iz/2+1)); ...
                                            NumVz_sign(fix(iy/2+1)+1, fix(ix/2+1)+1, fix(iz/2+1)); ...
                                            NumVz_sign(fix(iy/2+1)+1, fix(ix/2+1)  , fix(iz/2+1))];
                    
                    if      mod(iy,2)==0 & mod(ix,2)~=0
                        ind_val     =   [ 3/16 9/16  3/16 1/16];
                        
                    elseif  mod(iy,2)~=0 & mod(ix,2)~=0
                        ind_val     =   [ 1/16 3/16  9/16 3/16 ];
                        
                    elseif  mod(iy,2)~=0 & mod(ix,2)==0
                        ind_val     =   [ 3/16 1/16 3/16 9/16 ];
                        
                    elseif  mod(iy,2)==0 & mod(ix,2)==0
                        ind_val     =   [ 9/16 3/16 1/16 3/16 ];
                    end
                    
                    ind_val = ind_val(:).*sign;
                    
                    
                else
                    
                    % odd - trilinear interpolation in x,y,z direction
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1),   fix(ix/2+1)  , fix(iz/2) ) ];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1),   fix(ix/2+1)+1, fix(iz/2) ) ];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1)+1, fix(ix/2+1)+1, fix(iz/2) ) ];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1)+1, fix(ix/2+1)  , fix(iz/2) ) ];
                    
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1),   fix(ix/2+1)  , fix(iz/2)+1) ];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1),   fix(ix/2+1)+1, fix(iz/2)+1) ];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1)+1, fix(ix/2+1)+1, fix(iz/2)+1) ];
                    ind_coarse      =   [ind_coarse, NumVz_coarse(fix(iy/2+1)+1, fix(ix/2+1)  , fix(iz/2)+1) ];
                    
                     sign           =   [   NumVz_sign(fix(iy/2+1)  , fix(ix/2+1)  , fix(iz/2)  ); 	...
                                            NumVz_sign(fix(iy/2+1)  , fix(ix/2+1)+1, fix(iz/2)  ); 	...
                                            NumVz_sign(fix(iy/2+1)+1, fix(ix/2+1)+1, fix(iz/2)  ); 	...
                                            NumVz_sign(fix(iy/2+1)+1, fix(ix/2+1)  , fix(iz/2)  ); 	...
                                            NumVz_sign(fix(iy/2+1)  , fix(ix/2+1)  , fix(iz/2)+1); 	...
                                            NumVz_sign(fix(iy/2+1)  , fix(ix/2+1)+1, fix(iz/2)+1); 	...
                                            NumVz_sign(fix(iy/2+1)+1, fix(ix/2+1)+1, fix(iz/2)+1); 	...
                                            NumVz_sign(fix(iy/2+1)+1, fix(ix/2+1)  , fix(iz/2)+1)   ];

                    % in the y-direction it is exactly in the middle
                    if mod(iz,2)==0
                        z_fac=0.5;
                    else
                        z_fac=0.5;
                    end
                    
                    if      mod(iy,2)==0 & mod(ix,2)~=0
                        ind_val     =   [[3/16 9/16  3/16 1/16]*z_fac    [ 3/16 9/16  3/16 1/16]*(1-z_fac)] ;
                        
                    elseif  mod(iy,2)~=0 & mod(ix,2)~=0
                        ind_val     =   [[1/16 3/16  9/16 3/16]*z_fac    [ 1/16 3/16  9/16 3/16 ]*(1-z_fac)];
                        
                        
                    elseif  mod(iy,2)~=0 & mod(ix,2)==0
                        ind_val     =   [[3/16 1/16 3/16 9/16 ]*z_fac    [ 3/16 1/16 3/16 9/16 ]*(1-z_fac)];
                        
                    elseif  mod(iy,2)==0 & mod(ix,2)==0
                        ind_val     =   [[9/16 3/16 1/16 3/16 ]*z_fac    [ 9/16 3/16 1/16 3/16 ]*(1-z_fac)];
                    end
                    
                    
                    ind_val = ind_val(:).*sign; % take BC's into account
                    
                end
                
                
                % Store all info in vectors, which are later added to a sparse
                % matrix
                ind_i = [ind_i; ind_fine*ones(size(ind_coarse(:)))];
                ind_j = [ind_j; ind_coarse(:)];
                val   = [val; ind_val(:)];
                
                
                
                
                
            end
        end
    end
    
    
    iix=size(NumVz_fine,2);
    for iiy=1:size(NumVz_fine,1)
        for iiz=1:size(NumVz_fine,3) 
            
            ix = iix;
            iy = iiy;
            iz = iiz;
            
            % Vx  is averaged from 12 surrounding cells
            ind_fine    =   NumVz_fine(iy,ix,iz);
            ind_coarse  =   NumVz_coarse(fix(iy/2+1)  ,end,fix(iz/2+1)  );
            
            ind_i = [ind_i; ind_fine];
            ind_j = [ind_j; ind_coarse];
            val   = [val; 1];
            
        end
    end
    
    iiy=size(NumVz_fine,1);
    for iix=1:size(NumVz_fine,2)
        for iiz=1:size(NumVz_fine,3) % last plane is dummy one
            
            ix = iix;
            iy = iiy;
            iz = iiz;
            
            ind_fine    =   NumVz_fine(iy,ix,iz);
            ind_coarse  =   NumVz_coarse(end  ,fix(ix/2+1),fix(iz/2+1)  );
            
            ind_i = [ind_i; ind_fine];
            ind_j = [ind_j; ind_coarse];
            val   = [val; 1];
            
        end
    end
  
    
end


Pp              = sparse(ind_i, ind_j, val, max(MG_Info.Numbering{level_fine}.Number_P(:)),max(MG_Info.Numbering{level_coarse}.Number_P(:)));


% % Prolongation matrix is simply the transpose of the restriction one multiplied by a factor:
% 
% Sol_vec                     =   zeros(max(MG_Info.Numbering{level_coarse}.Number_P(:)),1);
% Sol_vec(NumVz_coarse(:))    =   1;
% 
% Vz_coarse                   =   MG_Info.Grids{level_coarse}.ZVz+10; %.*MG_Info.Grids{level_coarse}.YVx;
% Vz_coarse(:,end+1,:) = 3;
% Vz_coarse(end+1,:,:) = 1;
% 
% 
% % % extract 3D matrix
% Sol_vec(MG_Info.Numbering{level_coarse}.Number_Vz) = Vz_coarse;
% 
% 
% 
% %
% 
% % Testing the restriction operatator
% % P_fine                  = (MG_Info.Grids{level_fine}.ZP);        % a P-value at a fine grid
% 
% 
% 
% 
% % % Prolongate
% % sol_vec_coarse          =   ones(size(sol_vec_coarse));
% 
% %
% sol_fine_prolongated    =    Pp*Sol_vec;
% Vz_fine_prolongated      =   sol_fine_prolongated(MG_Info.Numbering{level_fine}.Number_Vz);
% 
% 
% Vz_fine_prolongated(:,:,2)

