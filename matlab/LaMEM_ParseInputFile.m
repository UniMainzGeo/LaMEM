function [npart_x,npart_y,npart_z, varargout] = LaMEM_ParseInputFile(File)
%
% This reads the LaMEM *.dat input file and retrieves 
%   # of particles in all directions
%   grid dimensions & left/front/bottom coordinates
%   # of elements in all directions
%
% Warning: this routine currently only works with a constant grid-spacing
% (to be modified!)


str             =   fileread(File);             % Read file as a full string
lines           =   regexp(str,'\n','split');   % split into lines

% We always need to know the # of markers in every cell
Grid.nmark_x    =   RetrieveValue(lines,'nmark_x');
Grid.nmark_y    =   RetrieveValue(lines,'nmark_y');
Grid.nmark_z    =   RetrieveValue(lines,'nmark_z');

if nargout>3
    % The grid info itself is only required for non-parallel cases
    
    Grid.nel_x      =   RetrieveValue(lines,'nel_x');
    Grid.nel_y      =   RetrieveValue(lines,'nel_y');
    Grid.nel_z      =   RetrieveValue(lines,'nel_z');
    
    Grid.coord_x    =   RetrieveValue(lines,'coord_x');
    Grid.coord_y    =   RetrieveValue(lines,'coord_y');
    Grid.coord_z    =   RetrieveValue(lines,'coord_z');
    
    Grid.W          =   diff(Grid.coord_x);
    Grid.L          =   diff(Grid.coord_y);
    Grid.H          =   diff(Grid.coord_z);
    
    % Create a 3D mesh, for each of the particles
    npart_x =   Grid.nmark_x;
    npart_y =   Grid.nmark_y;
    npart_z =   Grid.nmark_z;
    
    nump_x  =   Grid.nel_x*npart_x;
    nump_y  =   Grid.nel_y*npart_y;
    nump_z  =   Grid.nel_z*npart_z;
    
    W       =   Grid.W;
    L       =   Grid.L;
    H       =   Grid.H;
    
    dx      =   W/nump_x;
    dy      =   L/nump_y;
    dz      =   H/nump_z;
    x       =   [Grid.coord_x(1)    + dx*0.5 : dx : Grid.coord_x(2)  - dx*0.5 ];
    y       =   [Grid.coord_y(1)    + dy*0.5 : dy : Grid.coord_y(2)  - dy*0.5 ];
    z       =   [Grid.coord_z(1)    + dz*0.5 : dz : Grid.coord_z(2)  - dz*0.5 ];
    [X,Y,Z] =   meshgrid(x,y,z);
    
    
    varargout{1} = Grid;
    varargout{2} = X;
    varargout{3} = Y;
    varargout{4} = Z;
    varargout{5} = W;
    varargout{6} = L;
    varargout{7} = H;
    
else
    npart_x     =   Grid.nmark_x;
    npart_y     =   Grid.nmark_y;
    npart_z     =   Grid.nmark_z;
    varargout   =   {};
end



% -------------------------------------------------------------------------
function Value = RetrieveValue(lines, Keyword)
% find the line @ which the keyword occurs and parse it 

a       =   regexp(lines,Keyword);
ind     =   find(~cellfun(@isempty,a));     % line at which it occurs


str_1   =   lines{ind}; 

% 1) Filter out comments from string
id      =   findstr(str_1,'#');
if ~isempty(id)
    str_1(id(1):end) = [];
end

% 2) Determine where the '=' symbol is
id      =   findstr(str_1,'=');
str_1(1:id) = [];

% 3) Transfer to numerical value
Value   =   str2num(str_1);
