%% Convert CT-scans from .tif images to 3D raw image data

clc, clear all , close all

%% Add important paths

addpath('/home/anton/PROG/LaMEM/matlab')

%% Rock phase
rock = 0; % default rock phase ID

%% Define parameters

% Circular    = logical(1); % circular scans flag (square otherwise)
% SizeOfVoxel = 4.34e-6;    % micrometers
% Dir         = 'porous/';  % specify directory with your CT-data
% Ext         = '*.TIF';    % file exenstion in the directory

Circular    = logical(0);   % circular scans flag (square otherwise)
SizeOfVoxel = 3.528e-6;     % micrometers
Dir         = 'granular/';  % specify directory with your CT-data
Ext         = '*.tif';      % file exenstion in the directory

filename   = 'voxel3D';     % file-name to store voxel-based data set

%==========================================================================

%% open files

curdir = pwd;

cd(Dir);

files = dir(Ext);

NumInp = length(files);

if(Circular)

    %% read circular scans

    % compute diameters
    D = zeros(NumInp, 2);

    fprintf('Reading diameters ...\n')

    for i=1:NumInp
        
        % read image
        P = imread(files(i).name);
        
        % determine image crop bounds
        [row, col] = find(P == rock);

        nxmin = min(row);
        nxmax = max(row);
        nymin = min(col);
        nymax = max(col);

        D(i, 1) = nxmax-nxmin+1;
        D(i, 2) = nymax-nymin+1;

    end
    
    % compute mean diameter (take minimum, stay on a safe side)
    DA = min(D(:));
    
    % split square
    DAL = round(DA/2);
    DAR = DA - DAL; 
    
    % generate normalized pixel radius matix
    R = (DA-1)/2; 
    x = linspace(-R, R, DA);
    x = x/R;
    x = x.^2;
    
    RMAT = repmat(x, DA, 1) + repmat(x', 1, DA);

    % create corner index mask
    [corners] = find(RMAT >= 1.0);

    % read scans
    P3D = zeros(DA, DA, NumInp);
    
    fprintf('Reading scans ... \n')

    for i=1:NumInp
        
        % read image
        P = imread(files(i).name);

        % determine image crop bounds
        [row, col] = find(P == rock);

        nxmin = min(row);
        nxmax = max(row);
        nymin = min(col);
        nymax = max(col);

        % get center point
        CX = nxmin + round((nxmax-nxmin)/2);
        CY = nymin + round((nymax-nymin)/2);
        
        % correct bounds to crop exact square circumscribed around circle
        nxmin = CX - DAL + 1;
        nxmax = CX + DAR;
        nymin = CY - DAL + 1;
        nymax = CY + DAR;
              
        % crop image
        P = P(nxmin:nxmax, nymin:nymax);
        
        % apply corner mask
        P(corners) = rock;
        
        % store image scan
        P3D(:, :, i) = P;
        
    end

else
    
    %% read square scans

    P = imread(files(1).name);
    
    [DA, ~] = size(P);

    % read scans
    P3D = zeros(DA, DA, NumInp);
    
    fprintf('Reading scans ... \n')

    for i=1:NumInp
        
        % read image
        P = imread(files(i).name);

        % store image scan
        P3D(:, :, i) = P;
        
    end
end

% return to current directory
cd(curdir);

% set flags for fixed cells
if(rock == 0)
    P3D = ~P3D;
end

fprintf('Dumping raw image to disk ... \n')

fp = fopen(strcat(filename, '-', num2str(DA), 'x', num2str(DA), 'x', num2str(NumInp), '.raw'), 'w');

fwrite(fp, P3D, 'uint8'); 

fclose(fp);

% compute and print model size for creating partitioning file
R = SizeOfVoxel*DA/2; 
H = SizeOfVoxel*NumInp/2;

fprintf('Please create partitioning file with the following domain sizes: \n')

fprintf('coord_x = %g %g\n', -R, R); 
fprintf('coord_y = %g %g\n', -R, R);
fprintf('coord_z = %g %g\n', -H, H);
