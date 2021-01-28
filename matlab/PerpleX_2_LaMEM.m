%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perple_X to LaMEM (January 2021)
% Credit to:
% Lisa Rummel for writing the first version
% Andrea Piccolo for modifications
% Arne Spang for the latest version
%
% Important: 
%   1)  You need to copy this file to the Perple_X directory, where you did
%       your work
%   2)  using werami, you need to make sure to create an output file that
%       includes, in the following order:
%           bulk density     | melt fraction   | melt density
%
%   Once it is generated, you can inspect the diagram with
%   Visualize_LaMEM_PhaseDiagram('YOUR_PHASE_DIAGRAM_NAME_HERE.in')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% input
% name of your PerpleX output file:
filename = 'sediments_1.tab';
% name of the .mat and .in file to be produced (without ending (no .in)):
outname  = 'sediments_1';
% plot the phase diagram or not
plotFlag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% don't change anything after here

% things that might change some day
fake_rho_melt = 2700;
numCol        = 5;

% open file
fID = fopen(filename);

% loop through top of the file and check for expected format
for i = 1 : 13
  line = fgetl(fID);
  
  if i == 4 || i == 8 || i == 13
    line = strsplit(line);
    if ~strcmp(line(1),'T(K)') && ~strcmp(line(1),'P(bar)')
      error('Unexpected format!')
    end
  end

  if i == 7
    line = strsplit(line);
    numT = str2double(line(2));
  end
  
  if i == 11
    line = strsplit(line);
    numP = str2double(line(2));
  end
end

% preallocate vectors
nLines   = numT*numP;
ooM      = floor(log10(nLines));
factor   = round(nLines/(10^ooM));
len      = factor * 10^ooM;
tenPer   = len/10;
T        = zeros(nLines,1);
P        = zeros(nLines,1);
rho      = zeros(nLines,1);
f_melt   = zeros(nLines,1);
rho_melt = zeros(nLines,1);

% read the matrix
iter = 0;
fprintf('Reading values...\n')
for i = 1 : nLines
  line = fgetl(fID);
  line = strsplit(line);
  T(i)        = str2double(line(2));
  P(i)        = str2double(line(3));
  rho(i)      = str2double(line(4));
  f_melt(i)   = str2double(line(5));
  rho_melt(i) = str2double(line(6));
  if mod(i,tenPer) == 0
    iter = iter + 10;
    fprintf('%d%% done \n',iter)
  end
end
fprintf('100%% done! \n')

% remove NaNs
ind            = find(isnan(f_melt));
f_melt(ind)    = 0;
rho_melt(ind)  = fake_rho_melt;

% melt fraction should not be in per cent
f_melt    = f_melt/100;

% reshape to matrices and vectors
rho_bulk  = reshape(rho,numT,numP);
Melt      = reshape(f_melt,numT,numP);
rho_melt  = reshape(rho_melt,numT,numP);
rho_solid = (rho_bulk - Melt.*rho_melt)./(1-Melt);
T_K       = T(1:numT);
P_bar     = P([1:numT:(numP-1)*numT+1]);

% save .mat file
save([outname '.mat'],'P_bar','T_K','rho_bulk','rho_solid','Melt','rho_melt');

% plot
if plotFlag
  figure()
  pcolor(T_K,P_bar/1e3,rho_bulk'); shading interp
  colorbar
  ylabel('Pressure [kbar]')
  xlabel('Temperature [K]')
  title('Bulk density [kg/m^3]')
  
  figure()
  pcolor(T_K,P_bar/1e3,rho_solid'); shading interp
  colorbar
  ylabel('Pressure [kbar]')
  xlabel('Temperature [K]')
  title('Solid density [kg/m^3]')
  
  figure()
  inds              = find(rho_melt == fake_rho_melt);
  rhoMeltPlot       = rho_melt;
  rhoMeltPlot(inds) = NaN;
  pcolor(T_K,P_bar/1e3,rhoMeltPlot'); shading interp
  colorbar
  ylabel('Pressure [kbar]')
  xlabel('Temperature [K]')
  title('Melt density [kg/m^3]')
  
  figure()
  pcolor(T_K,P_bar/1e3,Melt'); shading interp
  colorbar
  ylabel('Pressure [kbar]')
  xlabel('Temperature [K]')
  title('Melt Content []')
end

% setup grid that works for LaMEM (buffer is limited to 44100 lines)
nN        = 200;
nT        = nN;
nP        = nN;
T_interp  = linspace(T_K(1),T_K(end),nT);
P_interp  = linspace(P_bar(1),P_bar(end),nP);
[P,T]     = meshgrid(P_bar,T_K);
[P2D,T2D] = meshgrid(P_interp,T_interp);

% interpolate onto that grid
fprintf('Interpolating to LaMEM grid...\n')
rho_interp      = griddata(T,P,rho_solid,T2D,P2D,'linear');
melt_interp     = griddata(T,P,Melt,T2D,P2D,'linear');
rho_melt_interp = griddata(T,P,rho_melt,T2D,P2D,'linear');
fprintf('Done!\n')

% vectorize grids for LaMEM
rhovec       = rho_interp(:);
Meltvec      = melt_interp(:);
rhoMeltvec   = rho_melt_interp(:);
Tvec         = T2D(:);              % [K]
Pvec         = P2D(:);              % [bar]

% assemble output matrix
MATRIX       = [rhoMeltvec, Meltvec, rhovec, Tvec, Pvec];

% grid stepping
dT1 = diff(T_K); 
dT  = mean(dT1); 
dP1 = diff(P_bar);
dP  = mean(dP1);

% structure of output file:
% % Line 1:           Number of columns (5)
% % Line 2-49:        room to insert comments
% % Line 50:          lowest temperature
% % Line 51:          temperature steps
% % Line 52:          number of nodes along T
% % Line 53:          lowest pressure
% % Line 54:          pressure steps
% % Line 55:          number of nodes along P
% % Line 56-end:      information in 5 columns
% % rho_melt, melt fraction, rho_solid, T [K], P [bar]

%% create output file
fprintf('Writing LaMEM input to %s...\n',[outname '.in'])

fileID     = fopen([outname '.in'],'w');

fprintf(fileID, '%d\n',numCol);     % number of columns
fprintf(fileID, '%s\n',filename);   % name of PerpleX output

% block for comments
for i = 1 : 47
  fprintf(fileID, '\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All the information needed    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, '%5f\n',T_K(1));
fprintf(fileID, '%5f\n',dT);
fprintf(fileID, '%d\n',nT);
fprintf(fileID, '%5f\n',P_bar(1));
fprintf(fileID, '%5f\n',dP);
fprintf(fileID, '%d\n',nP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple Loop to put all the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : nP*nT
    fprintf(fileID, '%5f %5f %5f %5f %5f\n',MATRIX(i,:));
end
fclose(fileID);

fprintf('Done! \n')
