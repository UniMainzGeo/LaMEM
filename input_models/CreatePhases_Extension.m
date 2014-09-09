%CreatePhases
%
% Create a 3D array with phase distributions

clear



W       =   40e3;
L       =   40e3;
H       =   10e3;

nump_x  =   101;
nump_y  =   101;
nump_z  =   31;





dx_reg  =   W/(nump_x-1);
dy_reg  =   L/(nump_y-1);
dz_reg  =   H/(nump_z-1);
x_left  =   -W/2;
y_front =   -L/2;
z_bot   =   0;

[X,Y,Z] =   meshgrid([x_left:dx_reg:x_left+W], [y_front:dy_reg:y_front+L], [z_bot:dz_reg:z_bot+H]);

Phase   =   zeros(size(X));     %   Contains phases
Temp    =   zeros(size(X));     %   Contains temperatures

if round((nump_y-1)/4)~=(nump_y-1)/4
    disp(['watch it: probably mistake with number of particles in y-direction'])
end
if round((nump_x-1)/4)~=(nump_x-1)/4
    disp(['watch it: probably mistake with number of particles in x-direction'])
end
if round((nump_z-1)/4)~=(nump_z-1)/4
    disp(['watch it: probably mistake with number of particles in z-direction'])
end


% Create one or two perturbations====================================
Setup = 2;      %   1-one pert. t2-two perturn

if Setup==1
    % one throughgoing perturbation
    ParticleOutput   =   'ExtensionSetup_OneThroughgoing.dat';

    
    % First perturbation:
    Xc  =   0;   Yc  =   -L/2;   Zc  =  0;
    
    W1          =   3e3;
    H1          =   3e3;
    L1          =   2*L;
    ind         =   find( abs((X-Xc))<W1/2 & abs((Y-Yc))<L1/2 & abs((Z-Zc))<H1/2 );
    Phase(ind)  =   1;         % weak perturbation
    
elseif Setup==2
    % one half throughgoing perturbation
    ParticleOutput   =   'ExtensionSetup_HalfThroughgoing.dat';
    
    % First perturbation:
      Xc  =   0;   Yc  =   -L/2;   Zc  =  0;
    
    W1          =   3e3;
    H1          =   3e3;
    L1          =   L;
    ind         =   find( abs((X-Xc))<W1/2 & abs((Y-Yc))<L1/2 & abs((Z-Zc))<H1/2 );
    Phase(ind)  =   1;         % weak perturbation
    
  
elseif Setup==3
    
    ParticleOutput   =   'ExtensionSetup_TwoQuarterThroughgoing.dat';
    
    % First perturbation:
    Xc  =   W/4;   Yc  =   -L/2;   Zc  =  0;
    W1          =   3e3;
    H1          =   3e3;
    L1          =   L/2;
    ind         =   find( abs((X-Xc))<W1/2 & abs((Y-Yc))<L1/2 & abs((Z-Zc))<H1/2 );
    Phase(ind)  =   1;         % weak perturbation
    
    
    Xc  =   -W/4;   Yc  =   L/2;   Zc  =  0;
    W1          =   3e3;
    H1          =   3e3;
    L1          =   L/2;
    ind         =   find( abs((X-Xc))<W1/2 & abs((Y-Yc))<L1/2 & abs((Z-Zc))<H1/2 );
    Phase(ind)  =   1;         % weak perturbation
    
    
    
end

%==========================================================================




%==========================================================================
% Set initial temperature distribution - in Celcius
%
Temp = (H-Z)./H*0 + 0.5 + (rand(size(Z))-0.5)*0.05 + 0*1000;

%
%==========================================================================



% Plot in matlab
figure(1), clf, hold on
marker = {'r','k','g','b','m','y'};
for phase=0:max(Phase(:))
    ind = find(Phase==phase);
    ind = find(abs(Phase-phase)<1e-9);
    ind_plot = phase+1;
    if ind_plot>length(marker); ind_plot = ind_plot-fix(ind_plot/length(marker))*length(marker)+1; end
    plot3(X(ind),Y(ind), Z(ind),[marker{ind_plot},'.']);
end
view(3)
axis equal





% Create a vector to be saved in PETSc format
PhaseOrig   = Phase;
PhaseVec(1) = nump_z;
PhaseVec(2) = nump_y;
PhaseVec(3) = nump_x;
Phase       = permute(Phase,[2 1 3]);
Temp        = permute(Temp, [2 1 3]);

PhaseVec    = [PhaseVec(:); Phase(:); Temp(:)];

% Save data to file
PetscBinaryWrite(ParticleOutput, PhaseVec);




% 
% 
% %==========================================================================
% % Generate a surface that describes the same data
% %==========================================================================
% 
% % Upper and lower boundaries
% Xtop    =   [];     Ytop    =   [];     Ztop    =   [];
% Xbot    =   [];     Ybot    =   [];     Zbot    =   [];
% InnerPhase  =     1;
% for ix=1:size(X,2);
%     for iy=1:size(X,1);
%         ind =   find(PhaseOrig(iy,ix,:)==InnerPhase);
% 
%         if ~isempty(ind)
%             deps = 1e-4;
%             Xtop    =   [Xtop; X(iy,ix,ind(end))];
%             Ytop    =   [Ytop; Y(iy,ix,ind(end))];
%             Ztop    =   [Ztop; Z(iy,ix,ind(end))-deps];
% 
%             Xbot    =   [Xbot; X(iy,ix,ind(1))];
%             Ybot    =   [Ybot; Y(iy,ix,ind(1))];
%             Zbot    =   [Zbot; Z(iy,ix,ind(1))];
%         end
% 
%     end
% end
% 
% 
% % Front and back boundaries
% Xfront   =   [];     Yfront   =   [];     Zfront   =   [];
% Xback    =   [];     Yback    =   [];     Zback    =   [];
% InnerPhase  =     1;
% for ix=1:size(X,2);
%     for iz=1:size(X,3);
%         ind =   find(PhaseOrig(:,ix,iz)==InnerPhase);
% 
%         if ~isempty(ind)
%             deps = 1e-4;
%             Xback    =   [Xback; X(ind(end),ix,iz)];
%             Yback    =   [Yback; Y(ind(end),ix,iz)-deps];
%             Zback    =   [Zback; Z(ind(end),ix,iz)];
% 
%             Xfront    =   [Xfront; X(ind(1),ix,iz)];
%             Yfront    =   [Yfront; Y(ind(1),ix,iz)+deps];
%             Zfront    =   [Zfront; Z(ind(1),ix,iz)];
%         end
% 
%     end
% end
% 
% % Left boundaries
% Xleft   =   [];     Yleft   =   [];     Zleft   =   [];
% InnerPhase  =     1;
% for iy=1:size(X,1);
%     for iz=1:size(X,3);
%         ind =   find(PhaseOrig(iy,:,iz)==InnerPhase);
% 
%         if ~isempty(ind)
%             if (ind(1))==1
%                 Xleft    =   [Xleft; X(iy,ind(1),iz)];
%                 Yleft    =   [Yleft; Y(iy,ind(1),iz)];
%                 Zleft    =   [Zleft; Z(iy,ind(1),iz)];
%             end
%         end
% 
%     end
% end
% 
% % 'Patch' the surfaces together
% TriTop          =   delaunay(Xtop,Ytop);
% TriBot          =   delaunay(Xbot,Ybot);
% 
% 
% % Create the front part----------------------------------------------------
% TriFront        =   delaunay(Xfront,Zfront);
% 
% % 'Inspect' the triangles. Only those that are inside are allowed
% TriReal = ones(length(TriFront),1);
% for itri = 1:length(TriFront);
% 
%     tri     =   TriFront(itri,:);
%     xmean   =   mean(Xfront(tri));
%     ymean   =   mean(Yfront(tri));
%     zmean   =   mean(Zfront(tri));
% 
%     x_int   =   round((xmean-x_left)/dx_reg)+1;
%     z_int   =   round((zmean-z_bot)/dz_reg)+1;
% 
%     if (x_int==0); x_int = 1; end
% 
%     if PhaseOrig(fix(size(PhaseOrig,1)/2),x_int,z_int)~=InnerPhase;
%         TriReal(itri) = 0;
%     end
% end
% ind = find(TriReal==0);
% TriFront(ind,:) = [];
% %--------------------------------------------------------------------------
% 
% % % Create the back part-----------------------------------------------------
% TriBack        =   delaunay(Xback,Zback);
% 
% % 'Inspect' the triangles. Only those that are inside are allowed
% TriReal = ones(length(TriBack),1);
% for itri = 1:length(TriBack);
% 
%     tri     =   TriBack(itri,:);
%     xmean   =   mean(Xback(tri));
%     ymean   =   mean(Yback(tri));
%     zmean   =   mean(Zback(tri));
% 
%     x_int   =   round((xmean-x_left)/dx_reg)+1;
%     z_int   =   round((zmean-z_bot)/dz_reg)+1;
% 
%     if PhaseOrig(fix(size(PhaseOrig,1)/2),x_int,z_int)~=InnerPhase;
%         TriReal(itri) = 0;
%     end
% end
% ind = find(TriReal==0);
% TriBack(ind,:) = [];
% %TriBack=[];
% % %--------------------------------------------------------------------------
% 
% TriLeft        =   delaunay(Yleft,Zleft);
% 
% TriCombined     =   [TriTop; TriBot+ length(Xtop); TriFront+length(Xtop)+length(Xbot); TriBack+length(Xtop)+length(Xbot)+length(Xfront); TriLeft+length(Xtop)+length(Xbot)+length(Xfront)+length(Xback)];
% 
% XCombined       =   [Xtop; Xbot; Xfront; Xback; Xleft];
% YCombined       =   [Ytop; Ybot; Yfront; Yback; Yleft];
% ZCombined       =   [Ztop; Zbot; Zfront; Zback; Zleft];
% ZCombined=ZCombined-0.005;
% 
% % Save data in Petsc format
% PassiveMat(:,1) =   XCombined;
% PassiveMat(:,2) =   YCombined;
% PassiveMat(:,3) =   ZCombined;
% PassiveMat      =   [length(XCombined);    XCombined(:);     YCombined(:);     ZCombined(:); ...
%     length(TriCombined);  TriCombined(:,1); TriCombined(:,2); TriCombined(:,3)];
% %PassiveMat      =   sparse(PassiveMat);
% PetscBinaryWrite(PassiveOutput, PassiveMat);
% 


figure(2), clf, hold on
%trisurf(TriTop,Xtop,Ytop,Ztop,'facecolor','r')
%trisurf(TriBot,Xbot,Ybot,Zbot,'facecolor','r')

%trisurf(TriCombined,XCombined,YCombined,ZCombined,'facecolor','r')

view(3), axis equal, box on
axis( [min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) min(Z(:)) max(Z(:)) ]);

view(0,0)

%save test TriCombined

