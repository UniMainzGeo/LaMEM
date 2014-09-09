function  SURFACE_DATA = ExtractInternalSurface(TimestepNumber,SimulationName)
% ExtractInternalSurface
%
% Converts the internal free surface from LaMEM from various VTK data sets into a
% MATLAB format so we can re-use it for further analysis
%
% Example usage:
%   TimestepNumber  =   10;
%   name            =   'DetFolding_4';
%   SURFACE_DATA    =   ExtractInternalSurface(TimestepNumber,name)
%


% reconstruct the directory number
num_str =  num2str(TimestepNumber+100000,'%i');
num_str(1)='0';
dirname =   ['Timestep_',num_str];

% reconstruct the full filename
name = [SimulationName,'_Surface_topography_step',num_str,'-mesh.pvts'];

% display information on screen
disp(['Reading the internal free surface of simulation: ',SimulationName])
disp(['From directory:                                  ',dirname])


dir_cur =   pwd;
cd(dirname)

% 1: Open the main file, that tells us where the subfiles are located
fid = fopen(name,'r');
if fid<0
    error(['Cannot open file ',name])
end
    
Piece = 0;
while ~feof(fid)
    line = fgetl(fid);
    
    % read line
    if strfind(line,'WholeExtent="')
        % size of full domain
        start_l       = strfind(line,'WholeExtent="');
        start_l       = start_l+13;
        SUBDOMAIN.WholeExtent   = sscanf(line(start_l:end),'%i');
    end
    
    
    if strfind(line,'<Piece Extent="')
        % size and filename of each of the subdomains
        Piece = Piece+1;
        
        % Interprete this line, with information about the subdomains.
        
        % Get the size of the subdomain:
        start_l       = strfind(line,'Piece Extent="');
        start_l       = start_l+14;
        SUBDOMAIN(Piece).PieceExtent   = sscanf(line(start_l:end),'%i'); % get the piece extend of this
        
        % Retrieve the name of the file we need to read
        start_l       =   strfind(line,'Source="');
        start_l       =   start_l+8;
        end_l         =   strfind(line,'"');
        end_l         =   end_l(end)-1;
        
        SUBDOMAIN(Piece).FileName      =   line(start_l:end_l);
    end
end
fclose(fid);


% Step 2: open each of the subfiles

for ifile=1:length(SUBDOMAIN)
    fid = fopen(SUBDOMAIN(ifile).FileName,'r');
    Piece = 0;
    num_PointDataArray = 1;
    while ~feof(fid)
        line = fgetl(fid);
        
        % Read the coordinates
        if strfind(line,'<Points>')
            
            % move two lines further
            line = fgetl(fid);
            num = 1;
            
            % read in the coordinates
            while isempty(strfind(line,'</Points>'))
                line = fgetl(fid);
                data = sscanf(line,'%f %f %f');
                if ~isempty(data)
                    SUBDOMAIN(ifile).Coordinates(num,:)   = data;
                end
                
                num=num+1;
            end
            
            % read in
            
        end
        
        % read the names of each of the fields that have point data
        if strfind(line,'<PointData Scalars="')
            start_l       =   strfind(line,'<PointData Scalars="');
            start_l       =   start_l+21;
            end_l         =   strfind(line,'"');
            end_l         =   end_l(end)-1;
            
            % Read in all the point data
            [a,count]     = sscanf(line(start_l:end_l),' %s ');
            a             = strsplit(line(start_l:end_l),' ');      % split string into substrings
            PointData     = a(1:count);
            numPointData  = length(PointData);
            SUBDOMAIN(ifile).PointData.Names = PointData;
            
        end
        
        
        if strfind(line,'Name="')
            % next line
            % line = fgetl(fid);
            
            % Read data
            num=1;
            while isempty(strfind(line,'</DataArray>'))
                line = fgetl(fid);
                data = sscanf(line,'%f');
                if ~isempty(data)
                    SUBDOMAIN(ifile).PointData.Data(num,num_PointDataArray)   = data;
                end
                num=num+1;
            end
            num_PointDataArray = num_PointDataArray+1;
        end
        
    end
    fclose(fid);
    
end

% At this stage, we have read in all the data of all subfiles, and added it
% to the SUBDOMAIN structure.

% Next we need to transfer these datasets of the subdomains into 2D arrays


% Initialize arrays
X2d = zeros((SUBDOMAIN(1).WholeExtent(2:2:end)'+1));
Y2d = zeros((SUBDOMAIN(1).WholeExtent(2:2:end)'+1));
Z2d = zeros((SUBDOMAIN(1).WholeExtent(2:2:end)'+1));
for num=1:numPointData
    DataArray{num} = zeros((SUBDOMAIN(1).WholeExtent(2:2:end)'+1));
    
end



for Piece = 1:length(SUBDOMAIN);
    
    
    
    % Figure out the number of elements in every direction:
    n   =  SUBDOMAIN(Piece).PieceExtent+1;
    Nel =  n(2:2:end)-n(1:2:end)+1;
    
    % Coordinates
    X2d_piece = reshape(SUBDOMAIN(Piece).Coordinates(:,1),Nel');
    X2d(n(1):n(2),n(3):n(4),n(5):n(6))  = X2d_piece;
    
    Y2d_piece = reshape(SUBDOMAIN(Piece).Coordinates(:,2),Nel');
    Y2d(n(1):n(2),n(3):n(4),n(5):n(6))  = Y2d_piece;
    
    Z2d_piece = reshape(SUBDOMAIN(Piece).Coordinates(:,3),Nel');
    Z2d(n(1):n(2),n(3):n(4),n(5):n(6))  = Z2d_piece;
    
    % Work on all Data arrays
    for num=1:numPointData
        Data_piece = reshape(SUBDOMAIN(Piece).PointData.Data(:,num),Nel');
        
        Data =  DataArray{num};
        
        Data(n(1):n(2),n(3):n(4),n(5):n(6))  = Data_piece;
        
        DataArray{num} = Data;
        
    end
    
    
end


% Finally, rename the data we want into a usefull format, and store them
% into a SURFACE_DATA structure

SURFACE_DATA.X2d = X2d;
SURFACE_DATA.Y2d = Y2d;
SURFACE_DATA.Z2d = Z2d;
for num=1:numPointData
    str = char(SUBDOMAIN(1).PointData.Names(num));
    str(strfind(str,'[')) = [];
    str(strfind(str,']')) = [];
    str(strfind(str,'''')) = [];
    
    str = ['SURFACE_DATA.',str,'=DataArray{num};'];
    eval(str);
    
end

% save as MATLAB file in the directory
save SURFACE_DATA SURFACE_DATA


cd(dir_cur);


