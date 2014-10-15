function [Rp] = RestrictionProlongationMatrixVelocity(level_fine, MG_Info);
%
% It interpolates the Vx values, defined at Vx-points, from a
% finer to a coarser mesh, using a technique recommended by a number of
% publications that is some kind of weighted linear interpolation.
%
% In addition, we need to take care that we have additional 'dummy' planes
% in the velocity matrixes

level_coarse    = level_fine-1;
Rp              = sparse(max(MG_Info.Numbering{level_coarse}.Number_P(:)),max(MG_Info.Numbering{level_fine}.Number_P(:)));


%% Vx-VELOCITY
NumVx_coarse     =  MG_Info.Numbering{level_coarse}.Number_Vx;
NumVx_fine       =  MG_Info.Numbering{level_fine}.Number_Vx;

% averaging scheme - the 4 cells at the same level receive the largest
% weight; the other 8 less.
%
% We have to take care of boundaries, which is why this routine is a bit
% uggly, but in this manner we can do all in a single loop

for iy=1:size(NumVx_coarse,1)-1     % last plane is dummy one
    for ix=1:size(NumVx_coarse,2)
        for iz=1:size(NumVx_coarse,3)-1 % last plane is dummy one
            
            % Vx  is averaged from 12 surrounding cells 
            ind_fine    = [];               ind_val = [];
            
            if iy>1 & iz>1
                ind_fine    = [ind_fine, NumVx_fine(2*iy-1,2*ix-1,2*iz-1)];     ind_val = [ind_val, 1/8];
            end
            if iy>1
                ind_fine    = [ind_fine, NumVx_fine(2*iy-1,2*ix-1,2*iz  )];     ind_val = [ind_val, 1/8];
            end
            
            if iz>1
                ind_fine    = [ind_fine, NumVx_fine(2*iy  ,2*ix-1,2*iz-1)];     ind_val = [ind_val, 1/8];
            end
            
            ind_fine    = [ind_fine, NumVx_fine(2*iy  ,2*ix-1,2*iz  )];     ind_val = [ind_val, 1/8];
            
            if ix>1                     % takes care of left boundary
                if iy>1 & iz>1
                    ind_fine    = [ind_fine, NumVx_fine(2*iy-1,2*ix-2,2*iz-1)];     ind_val = [ind_val, 1/16];
                end
                
                if iz>1
                ind_fine    = [ind_fine, NumVx_fine(2*iy  ,2*ix-2,2*iz-1)];     ind_val = [ind_val, 1/16];
                end
                
                ind_fine    = [ind_fine, NumVx_fine(2*iy  ,2*ix-2,2*iz  )];     ind_val = [ind_val, 1/16];
                if iy>1
                    ind_fine    = [ind_fine, NumVx_fine(2*iy-1,2*ix-2,2*iz  )];     ind_val = [ind_val, 1/16];
                end
                
            end
            
            if ix<size(NumVx_coarse,2)  % takes care of right boundary
                if iy>1 & iz>1
                    ind_fine    = [ind_fine, NumVx_fine(2*iy-1,2*ix  ,2*iz-1)];     ind_val = [ind_val, 1/16];
                end
                
                if iz>1
                    ind_fine    = [ind_fine, NumVx_fine(2*iy  ,2*ix  ,2*iz-1)];     ind_val = [ind_val, 1/16];
                end
                
                ind_fine    = [ind_fine, NumVx_fine(2*iy  ,2*ix  ,2*iz  )];     ind_val = [ind_val, 1/16];
                if iy>1
                    
                    ind_fine    = [ind_fine, NumVx_fine(2*iy-1,2*ix  ,2*iz  )];     ind_val = [ind_val, 1/16];
                end
                
            end
            ind_coarse  =   NumVx_coarse(iy,ix,iz);
            ind_val     =   ind_val./sum(ind_val);      % normalize to 1
            
            
            Rp(ind_coarse, ind_fine(:)) = ind_val;
            
        end
    end
end

% Deal with the dummy planes
iy=size(NumVx_coarse,1);
for ix=1:size(NumVx_coarse,2)
    for iz=1:size(NumVx_coarse,3) % last plane is dummy one
        
        ind_fine    = NumVx_fine(2*iy-1,2*ix-1,2*iz-1);     ind_val = 1;
        
        ind_coarse  =   NumVx_coarse(iy,ix,iz);
        Rp(ind_coarse, ind_fine(:)) = ind_val;
        
    end
end

for iy=1:size(NumVx_coarse,1)     % last plane is dummy one
    for ix=1:size(NumVx_coarse,2)
        iz          =   size(NumVx_coarse,3); 
        
        ind_fine    =   NumVx_fine(2*iy-1,2*ix-1,2*iz-1);     ind_val = 1;
        ind_coarse  =   NumVx_coarse(iy,ix,iz);
        Rp(ind_coarse, ind_fine(:)) = ind_val;
         
    end
end

% full(sum(Rp(NumVx_coarse(:),:),2))

%% Vy-VELOCITY
NumVy_coarse     =  MG_Info.Numbering{level_coarse}.Number_Vy;
NumVy_fine       =  MG_Info.Numbering{level_fine}.Number_Vy;
for iy=1:size(NumVy_coarse,1)     
    for ix=1:size(NumVy_coarse,2)-1         % last plane is dummy one
        for iz=1:size(NumVy_coarse,3)-1     % last plane is dummy one
            
            % Vy  is averaged from 12 surrounding cells 
            ind_fine    = [];               ind_val = [];
            
            if iy>1
                if ix>1 & iz>1
                    ind_fine    = [ind_fine, NumVy_fine(2*iy-2,2*ix-1,2*iz-1)];     ind_val = [ind_val, 1/16];
                end
                
                ind_fine        = [ind_fine, NumVy_fine(2*iy-2,2*ix  ,2*iz  )];     ind_val = [ind_val, 1/16];
                if iz>1
                    ind_fine    = [ind_fine, NumVy_fine(2*iy-2,2*ix  ,2*iz-1)];     ind_val = [ind_val, 1/16];
                end
                if ix>1
                    ind_fine    = [ind_fine, NumVy_fine(2*iy-2,2*ix-1,2*iz  )];     ind_val = [ind_val, 1/16];
                end
                
            end
            
            if ix>1 & iz>1
                ind_fine    = [ind_fine, NumVy_fine(2*iy-1,2*ix-1,2*iz-1)];     ind_val = [ind_val, 1/8];
            end
                ind_fine    = [ind_fine, NumVy_fine(2*iy-1,2*ix  ,2*iz  )];     ind_val = [ind_val, 1/8];
            if iz>1
                ind_fine    = [ind_fine, NumVy_fine(2*iy-1,2*ix  ,2*iz-1)];     ind_val = [ind_val, 1/8];
            end
            if ix>1
                ind_fine    = [ind_fine, NumVy_fine(2*iy-1,2*ix-1,2*iz  )];     ind_val = [ind_val, 1/8];
            end
            
            if iy<size(NumVy_coarse,1)
                if ix>1 & iz>1
                    ind_fine    = [ind_fine, NumVy_fine(2*iy  ,2*ix-1,2*iz-1)];     ind_val = [ind_val, 1/16];
                end
                ind_fine    = [ind_fine, NumVy_fine(2*iy  ,2*ix  ,2*iz  )];     ind_val = [ind_val, 1/16];
                if iz>1
                    ind_fine    = [ind_fine, NumVy_fine(2*iy  ,2*ix  ,2*iz-1)];     ind_val = [ind_val, 1/16];
                end
                if ix>1
                    ind_fine    = [ind_fine, NumVy_fine(2*iy  ,2*ix-1,2*iz  )];     ind_val = [ind_val, 1/16];
                end
            end
            
            
            ind_coarse  =   NumVy_coarse(iy,ix,iz);
            ind_val     =   ind_val./sum(ind_val);      % normalize to 1
            
            
            Rp(ind_coarse, ind_fine(:)) = ind_val;
            
        end
    end
end


% Deal with the dummy planes
ix=size(NumVy_coarse,2);
for iy=1:size(NumVy_coarse,1)
    for iz=1:size(NumVy_coarse,3) % last plane is dummy one
        
        ind_fine    = NumVy_fine(2*iy-1,2*ix-1,2*iz-1);     ind_val = 1;
        
        ind_coarse  =   NumVy_coarse(iy,ix,iz);
        Rp(ind_coarse, ind_fine(:)) = ind_val;
    end
end


for iy=1:size(NumVy_coarse,1)     % last plane is dummy one
    for ix=1:size(NumVy_coarse,2)
        iz          =   size(NumVy_coarse,3); 
        
        ind_fine    =   NumVy_fine(2*iy-1,2*ix-1,2*iz-1);     ind_val = 1;
        ind_coarse  =   NumVy_coarse(iy,ix,iz);
        Rp(ind_coarse, ind_fine(:)) = ind_val;
         
    end
end


% full(sum(Rp(NumVy_coarse(:),:),2))




%% Vz-VELOCITY
NumVz_coarse     =  MG_Info.Numbering{level_coarse}.Number_Vz;
NumVz_fine       =  MG_Info.Numbering{level_fine}.Number_Vz;
for iy=1:size(NumVz_coarse,1)-1     % last plane is dummy one
    for ix=1:size(NumVz_coarse,2)-1         % last plane is dummy one
        for iz=1:size(NumVz_coarse,3)     
            
            % Vz  is averaged from 12 surrounding cells 
            ind_fine    = [];               ind_val = [];
            
            if iz>1
                if iy>1 & ix>1
                    ind_fine    = [ind_fine, NumVz_fine(2*iy-1,2*ix-1,2*iz-2)];     ind_val = [ind_val, 1/16];
                end
                if iy>1
                    ind_fine    = [ind_fine, NumVz_fine(2*iy-1,2*ix  ,2*iz-2)];     ind_val = [ind_val, 1/16];
                end
                ind_fine        = [ind_fine, NumVz_fine(2*iy  ,2*ix  ,2*iz-2)];     ind_val = [ind_val, 1/16];
                if ix>1
                    ind_fine    = [ind_fine, NumVz_fine(2*iy  ,2*ix-1,2*iz-2)];     ind_val = [ind_val, 1/16];
                end
            end
            
            if ix>1 & iy>1
                ind_fine    = [ind_fine, NumVz_fine(2*iy-1,2*ix-1,2*iz-1)];     ind_val = [ind_val, 1/8];
            end
            if iy>1
                ind_fine    = [ind_fine, NumVz_fine(2*iy-1,2*ix  ,2*iz-1)];     ind_val = [ind_val, 1/8];
            end
            ind_fine        = [ind_fine, NumVz_fine(2*iy  ,2*ix  ,2*iz-1)];     ind_val = [ind_val, 1/8];
            if ix>1
                ind_fine    = [ind_fine, NumVz_fine(2*iy  ,2*ix-1,2*iz-1)];     ind_val = [ind_val, 1/8];
            end
            
            if iz<size(NumVz_coarse,3)
                if ix>1 & iy>1
                    ind_fine    = [ind_fine, NumVz_fine(2*iy-1,2*ix-1,2*iz  )];     ind_val = [ind_val, 1/16];
                end
                if iy>1
                    ind_fine    = [ind_fine, NumVz_fine(2*iy-1,2*ix  ,2*iz  )];     ind_val = [ind_val, 1/16];
                end
                
                ind_fine        = [ind_fine, NumVz_fine(2*iy  ,2*ix  ,2*iz  )];     ind_val = [ind_val, 1/16];
                if ix>1
                    ind_fine    = [ind_fine, NumVz_fine(2*iy  ,2*ix-1,2*iz  )];     ind_val = [ind_val, 1/16];
                end
            end
            
            
            ind_coarse  =   NumVz_coarse(iy,ix,iz);
            ind_val     =   ind_val./sum(ind_val);      % normalize to 1
            
            
            Rp(ind_coarse, ind_fine(:)) = ind_val;
            
        end
    end
end


% Deal with the dummy planes
ix=size(NumVz_coarse,2);
for iy=1:size(NumVz_coarse,1)
    for iz=1:size(NumVz_coarse,3) % last plane is dummy one
        
        ind_fine    = NumVz_fine(2*iy-1,2*ix-1,2*iz-1);     ind_val = 1;
        
        ind_coarse  =   NumVz_coarse(iy,ix,iz);
        Rp(ind_coarse, ind_fine(:)) = ind_val;
    end
end

% Deal with the dummy planes
iy=size(NumVz_coarse,1);
for ix=1:size(NumVz_coarse,2)
    for iz=1:size(NumVz_coarse,3) % last plane is dummy one
        
        ind_fine    = NumVz_fine(2*iy-1,2*ix-1,2*iz-1);     ind_val = 1;
        
        ind_coarse  =   NumVz_coarse(iy,ix,iz);
        Rp(ind_coarse, ind_fine(:)) = ind_val;
        
    end
end


%% Pressure
% We use simple injection for this
NumP_coarse     =  MG_Info.Numbering{level_coarse}.Number_P;
NumP_fine       =  MG_Info.Numbering{level_fine}.Number_P;
for iy=1:size(NumP_coarse,1)
    for ix=1:size(NumP_coarse,2)
        for iz=1:size(NumP_coarse,3)
            
            % Pressure is averaged from the 8 surrounding cells
            ind_fine    = NumP_fine(2*iy-1:2*iy,2*ix-1:2*ix,2*iz-1:2*iz);   
            ind_coarse  = NumP_coarse(iy,ix,iz);
            
            Rp(ind_coarse, ind_fine(:)) = 1/8;
            
        end
    end
end






% full(sum(Rp(NumVz_coarse(:),:),2))





% Prolongation matrix is simply the transpose of the restriction one multiplied by a factor:
% Pp                      =   Rp'*8;
% 
% 
% 
% % Testing the restriction operatator
% Vx_fine                  = (MG_Info.Grids{level_fine}.ZVx);        % a P-value at a fine grid
% 
% % add dummy planes to coordinates (for debugging)
% Vx_fine(end+1,:,:) = Vx_fine(end,:,:)*1e100 ;
% Vx_fine(:,:,end+1) = Vx_fine(:,:,end)*1e100 ;
% 
% 
% % Set P-solution in solution vector
% Sol_vec                 =   zeros(size(MG_Info.Matrixes{level_fine}.Rhs_vec));
% Sol_vec(NumVx_fine(:))   =   Vx_fine(:);
% 
% 
% % Restrict
% sol_vec_coarse          =   Rp*Sol_vec;
% 
% 
% % extract 3D matrix
% Vx_coarse                =   sol_vec_coarse(MG_Info.Numbering{level_coarse}.Number_Vx);
% 


% % Prolongate 
% sol_vec_coarse          =   ones(size(sol_vec_coarse));

% 
% sol_fine_prolongated    =   Pp'*sol_vec_coarse;
% P_fine_prolongated      =   sol_fine_prolongated(MG_Info.Numbering{level_fine}.Number_P);
% 


% 
% Vx_fine_real        =   Vx_fine(1:end-1,:,1:end-1);
% Vx_coarse_real      =   Vx_coarse(1:end-1,:,1:end-1);



